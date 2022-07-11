/*  File: modmap.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Feb  1 18:28 2019 (rd109)
 * Created: Sat Oct 27 20:37:44 2018 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"
#include <ctype.h>
#include "utils.h"
#include "seqhash.h"

#ifdef OMP
#include <omp.h>
#endif

struct {
  int k ;
  int w ;
  int s ;
  int B ;
} params ;

/*******************************************************************/

/* object to hold a reference sequence */

  
typedef struct {
  Seqhash *hasher ;
  int tableBits ;		/* max 34 so size < 2^32 so index is 32bit */
  U32 size ;			/* size of value array, and other arrays over this set */
  U64 tableSize ; 		/* = 1 << tableBits */
  U64 tableMask ;		/* = tableSize - 1 */
  U32 *index ;			/* this is the primary table - size tableSize */
  U64 *value ;			/* the hashed values */
  U16 *depth ;			/* depth at each index */
  U8  *info ;			/* bits for various things */
  U32 max ;			/* number of entries in the set - must be less than size */
} Modset ;

typedef struct {
  Modset *ms ;
  U32 size ;			/* size of index, offset, id */
  U32 max ;		      	/* number of indices in the reference */
  U32 *index ; 			/* consecutive Modset index values */
  U32 *offset ;			/* base offset within relevant sequence */
  U32 *id ;			/* index into refDict */
  U32 *depth ;		  	/* number of times this index is seen; size ms->size */
  U32 *rev ;			/* reverse index from mod to ref; size max */
  U32 *loc ;			/* offsets into rev for each mod ; size ms->max */
  DICT *dict ;			/* set of reference names */
  Array len ;			/* of U32, indexed over refId */
} Reference ;

int numThreads = 1 ;		/* default to serial - reset if multi-threaded */
FILE *outFile ;			/* initialise to stdout at start of main() */
BOOL isVerbose = FALSE ;
  

static inline void msSetCopy0 (Modset *ms, U32 i) { ms->info[i] &= 0xfc ; }
static inline void msSetCopy1 (Modset *ms, U32 i) { ms->info[i] = (ms->info[i] & 0xfc) | 1 ; }
static inline void msSetCopy2 (Modset *ms, U32 i) { ms->info[i] = (ms->info[i] & 0xfc) | 2 ; }
static inline void msSetCopyM (Modset *ms, U32 i) { ms->info[i] |= 3 ; }
static inline BOOL msIsCopy0  (Modset *ms, U32 i) { return ((ms->info[i] & 3) == 0) ; }
static inline BOOL msIsCopy1  (Modset *ms, U32 i) { return ((ms->info[i] & 3) == 1) ; }
static inline BOOL msIsCopy2  (Modset *ms, U32 i) { return ((ms->info[i] & 3) == 2) ; }
static inline BOOL msIsCopyM  (Modset *ms, U32 i) { return ((ms->info[i] & 3) == 3) ; }
static inline int  msCopy (Modset *ms, U32 i) { return (ms->info[i] & 3) ; }


Modset *modsetCreate (Seqhash *sh, int bits, U32 size)
{
  if (bits < 20 || bits > 34) die ("table bits %d must be between 20 and 34", bits) ;
  Modset *ms = new0 (1, Modset) ;
  ms->hasher = sh ;
  ms->tableBits = bits ;
  ms->tableSize = (U64)1 << ms->tableBits ;
  ms->tableMask = ms->tableSize - 1 ;
  ms->index = new0 (ms->tableSize, U32) ;
  if (size >= (ms->tableSize >> 2)) die ("Modset size %u is too big for %d bits", size, bits) ;
  else if (size) ms->size = size ;
  else ms->size = (ms->tableSize >> 2) - 1 ;
  ms->value = new (ms->size, U64) ;
  ms->depth = new0 (ms->size, U16) ;
  ms->info = new0 (ms->size, U8) ;
  return ms ;
}

void modsetDestroy (Modset *ms)
{ free (ms->index) ; free (ms->value) ; free (ms->info) ; free (ms) ; }

BOOL modsetPack (Modset *ms)	/* compress per-item arrays */
{ if (ms->size == ms->max+1) return FALSE ;
  resize (ms->value, ms->size, ms->max+1, U64) ;
  resize (ms->depth, ms->size, ms->max+1, U16) ;
  resize (ms->info, ms->size, ms->max+1, U8) ;
  ms->size = ms->max+1 ;
  return TRUE ;
}

U32 modsetIndexFind (Modset *ms, U64 hash, int isAdd)
{
  U64 diff = 0 ;
  U64 offset = hash & ms->tableMask ;
  U32 index = ms->index[offset] ;
  while (index && (ms->value[index] != hash))
    { if (!diff) diff = ((hash >> ms->tableBits) & ms->tableMask) | 1 ; /* odd so comprime */
      offset = (offset + diff) & ms->tableMask ;
      index = ms->index[offset] ;
    }
  if (!index && isAdd)
    { index = ms->index[offset] = ++ms->max ;
      if (ms->max >= ms->size) die ("hashTableSize %u is too small for %u", ms->size, ms->max) ;
      ms->value[index] = hash ;
    }
  return index ;
}

void modsetDepthPrune (Modset *ms, int min, int max)
{
  U32 i ;
  U32 N = ms->max ; ms->max = 0 ;
  memset (ms->index, 0, ms->tableSize*sizeof(U32)) ;
  for (i = 1 ; i <= N ; ++i)	/* NB index runs from 1..max */
    if (ms->depth[i] >= min && (!max || ms->depth[i] < max))
      { modsetIndexFind (ms, ms->value[i], TRUE) ;
	ms->info[ms->max] = ms->info[i] ;
	ms->depth[ms->max] = ms->depth[i] ;
      }
  fprintf (stderr, "  pruned Modset from %d to %d with min %d <= depth < max %d\n",
	   N, ms->max, min, max) ;
}

void modsetWrite (Modset *ms, FILE *f)
{ if (fwrite ("MSHSTv1",8,1,f) != 1) die ("failed to write modset header") ;
  if (fwrite (&ms->tableBits,sizeof(int),1,f) != 1) die ("failed to write bits") ;
  U32 size = ms->max+1 ; if (fwrite (&size,sizeof(U32),1,f) != 1) die ("failed to write size") ;
  seqhashWrite (ms->hasher, f) ;
  if (fwrite (ms->index,sizeof(U32),ms->tableSize,f) != ms->tableSize) die ("fail write index") ;
  if (fwrite (ms->value,sizeof(U64),ms->max+1,f) != ms->max+1) die ("failed to write value") ;
  if (fwrite (ms->depth,sizeof(U16),ms->max+1,f) != ms->max+1) die ("failed to write depth") ;
  if (fwrite (ms->info,sizeof(U8),ms->max+1,f) != ms->max+1) die ("failed to write info") ;
}

Modset *modsetRead (FILE *f)
{ char name[8] ;
  if (fread (name,8,1,f) != 1) die ("failed to read modset header") ;
  if (strcmp (name, "MSHSTv1")) die ("bad reference header") ;
  int bits ; if (fread (&bits,sizeof(int),1,f) != 1) die ("failed to read bits") ;
  U32 size ; if (fread (&size,sizeof(U32),1,f) != 1) die ("failed to read size") ;
  Seqhash *sh = seqhashRead (f) ;
  Modset *ms = modsetCreate (sh, bits, size) ;
  if (fread (ms->index,sizeof(U32),ms->tableSize,f) != ms->tableSize) die ("failed read index") ;
  if (fread (ms->value,sizeof(U64),size,f) != size) die ("failed to read value") ;
  if (fread (ms->depth,sizeof(U16),size,f) != size) die ("failed to read depth") ;
  if (fread (ms->info,sizeof(U8),size,f) != size) die ("failed to read info") ;
  ms->max = size - 1 ;
  return ms ;
}

BOOL modsetMerge (Modset *ms1, Modset *ms2)
{
  U32 i, index ;
  /* first check that the hashers are identical */
  Seqhash *sh1 = ms1->hasher, *sh2 = ms2->hasher ;
  if (sh1->w != sh2->w || sh1->k != sh2->k || sh1->factor1 != sh2->factor1) return FALSE ;
  /* need to expand size of ms1 to make space */
  U64 newSize = ms1->max + ms2->max + 1 ;
  if (newSize >= (ms1->tableSize >> 2)) newSize = (ms1->tableSize >> 2) - 1 ;
  resize (ms1->value, ms1->size, newSize, U64) ; 
  resize (ms1->depth, ms1->size, newSize, U16) ; 
  resize (ms1->info, ms1->size, newSize, U8) ;
  ms1->size = newSize ;
  /* then pass through ms2 adding into ms1 */
  for (i = 1 ; i <= ms2->max ; ++i)
    { index = modsetIndexFind (ms1, ms2->value[i], TRUE) ;
      U32 d = ms1->depth[index] + (U32) ms2->depth[i] ; if (d > U16MAX) d = U16MAX ;
      ms1->depth[index] = d ;
      int c = msCopy(ms1,index) + msCopy(ms2,i) ; if (c > 3) c = 3 ;
      ms1->info[index] &= 0x3 ; ms1->info[index] |= c ;
    }
  return TRUE ;
}

void modsetSummary (Modset *ms, FILE *f)
{
  seqhashReport (ms->hasher, f) ;
  fprintf (f, "MS table bits %d size %lu number of entries %u",
	   ms->tableBits, ms->tableSize, ms->max) ;
  if (!ms->max) { fputc ('\n', f) ; return ; }
  U32 i, copy[4] ; copy[0] = copy[1] = copy[2] = copy[3] = 0 ;
  Array h = arrayCreate (256, U32) ;
  for (i = 1 ; i <= ms->max ; ++i)
    { if (ms->depth[i] < arrayMax(h)) ++arr(h,ms->depth[i],U32) ; /* more efficient to check */
      else ++array(h,ms->depth[i],U32) ;
      ++copy[msCopy(ms,i)] ;
    }
  U64 sum = 0, tot = 0 ;
  for (int i = 0 ; i < arrayMax(h) ; ++i) { sum += arr(h,i,U32) ; tot += i * arr(h,i,U32) ; }
  I64 htot = tot / 2 ;
  for (int i = 0 ; i < arrayMax(h) ; ++i) { htot -= i*arr(h,i,U32) ; if (htot < 0) break ; }
  fprintf (f, " total count %lu\nMS average depth %.1f N50 depth %u",
	   tot, tot / (double)sum , i) ;
  if (copy[0] < ms->max)
    fprintf (f, " copy0 %u copy1 %u copy2 %u copyM %u", copy[0], copy[1], copy[2], copy[3]) ;
  fputc ('\n', f) ;
  arrayDestroy (h) ;
}

Reference *referenceCreate (Modset *ms, U32 size)
{
  if (!ms || !ms->size) die ("modset must be initialised before reference") ;
  if (!size) die ("refCreate must have size > 0") ;
  Reference *ref = new0 (1, Reference) ;
  ref->ms = ms ;
  ref->depth = ms->max ? new0 (ms->max, U32) : new0 (ms->size, U32) ;
  ref->size = size ;
  ref->index = new (size, U32) ;
  ref->offset = new (size, U32) ;
  ref->id = new (size, U32) ;
  ref->dict = dictCreate (1024) ;
  ref->len = arrayCreate (1024, U32) ;
  /* NB don't create loc and rev here */
  return ref ;
}
typedef struct { U32 index ; U32 pos ; } Seed ;

void referenceDestroy (Reference *ref)
{ free (ref->depth) ;
  if (ref->loc) free (ref->loc) ; 
  if (ref->rev) free (ref->rev) ;
  free (ref->index) ; free (ref->offset) ; free (ref->id) ;
  dictDestroy (ref->dict) ; arrayDestroy (ref->len) ;
  free (ref) ;
}

void referencePack (Reference *ref)
{
  resize (ref->depth, ref->ms->size, ref->ms->max+1, U32) ;
  resize (ref->index, ref->size, ref->max, U32) ;
  resize (ref->offset, ref->size, ref->max, U32) ;
  resize (ref->id, ref->size, ref->max, U32) ;
  ref->size = ref->max ;
  
  ref->rev = new (ref->size, U32) ;			/* now build rev and loc */
  ref->loc = new (ref->ms->max + 1, U32) ; ref->loc[0] = 0 ;
  // U32 i ;
  for (U32 i = 1 ; i <= ref->ms->max ; ++i)
    { ref->loc[i] = ref->loc[i-1] + ref->depth[i-1] ; }
  memset (ref->depth, 0, ref->ms->size*sizeof(U32)) ; /* build this up again in loop below */
  U32 *ri = ref->index ;
  for (U32 i = 0 ; i < ref->max ; ++i, ++ri)
    ref->rev[ref->loc[*ri] + ref->depth[*ri]++] = i ;
}

void referenceFastaRead (Reference *ref, char *filename, BOOL isAdd)
{
  U64 totLen = 0 ;
  
  dna2indexConv['N'] = dna2indexConv['n'] = 0 ; /* to get 2-bit encoding */
  SeqIO *si = seqIOopen (filename, dna2indexConv, FALSE) ;
  if (!si) die ("failed to read reference sequence file %s", filename) ;
  while (seqIOread (si)) 
    { int id ;
      if (!dictAdd (ref->dict, sqioId(si), &id))
	die ("duplicate ref sequence name %s", sqioId(si)) ;
      array (ref->len, id, int) = si->seqLen ;
      totLen += si->seqLen ;
      SeqhashRCiterator *mi = modRCiterator (ref->ms->hasher, sqioSeq(si), si->seqLen) ;
      U64 hash ; int pos ;
      while (modRCnext (mi, &hash, &pos, 0))
	{ U32 index = modsetIndexFind (ref->ms, hash, isAdd) ;
	  if (index)
	    { if (ref->max+1 >= ref->size) die ("reference size overflow") ;
	      ref->index[ref->max] = index ;
	      ++ref->depth[index] ;
	      ref->offset[ref->max] = pos ;
	      ref->id[ref->max] = id ;
	      ++ref->max ;
	    }
	}
      seqhashRCiteratorDestroy (mi) ;
    }
  seqIOclose (si) ;

  fprintf (outFile, "  %d hashes from %d reference sequences, total length %lu\n",
	   ref->max, dictMax(ref->dict), totLen) ;
  U32 i, *d = &ref->depth[1] ; U32 n1 = 0, n2 = 0, nM = 0 ;
  for (i = 1 ; i <= ref->ms->max ; ++i, ++d)
    if (*d == 1) { msSetCopy1 (ref->ms, i) ; ++n1 ; }
    else if (*d == 2) { msSetCopy2 (ref->ms, i) ; ++n2 ; }
    else { msSetCopyM (ref->ms, i) ; ++nM ; }
  fprintf (outFile, "  %d copy 1, %d copy 2, %d multiple\n", n1, n2, nM) ;

  if (isAdd) modsetPack (ref->ms) ;
  referencePack (ref) ;
}

void referenceWrite (Reference *ref, char *root)
{
  FILE *f ;
  if (!(f = fopenTag (root, "mod", "w"))) die ("failed to open %s.mod to write", root) ;
  modsetWrite (ref->ms, f) ;
  fclose (f) ;
  if (!(f = fopenTag (root, "ref", "w"))) die ("failed to open %s.ref to write", root) ;
  if (fwrite ("RFMSHv1",8,1,f) != 1) die ("failed to write reference header") ;
  U32 size = ref->max ; if (fwrite (&size,sizeof(U32),1,f) != 1) die ("failed to write size") ;
  if (fwrite (&ref->max,sizeof(U32),1,f) != 1) die ("failed to write max") ;
  if (fwrite (ref->index,sizeof(U32),size,f) != size) die ("failed write ref index") ;
  if (fwrite (ref->offset,sizeof(U32),size,f) != size) die ("failed write ref offset") ;
  if (fwrite (ref->id,sizeof(U32),size,f) != size) die ("failed write ref id") ;
  if (fwrite (ref->depth,sizeof(U32),ref->ms->max+1,f) != ref->ms->max+1) die ("fail depth") ;
  if (fwrite (ref->rev,sizeof(U32),size,f) != size) die ("fail rev") ;
  if (fwrite (ref->loc,sizeof(U32),ref->ms->max+1,f) != ref->ms->max+1) die ("fail loc") ;
  if (!arrayWrite (ref->len, f)) die ("failed write ref len") ;
  if (!dictWrite (ref->dict, f)) die ("failed write ref dict") ;
  /* write the dict last, because it can terminate in arbitrary byte alignment */
  fclose (f) ;
}  

Reference *referenceRead (char *root)
{
  FILE *f ;
  if (!(f = fopenTag (root, "mod", "r"))) die ("failed to open %s.mod to read", root) ;
  Modset *ms = modsetRead (f) ;
  fclose (f) ;
  if (!(f = fopenTag (root, "ref", "r"))) die ("failed to open %s.ref to read", root) ;
  char name[8] ;
  if (fread (name,8,1,f) != 1) die ("failed to read reference header") ;
  if (strcmp (name, "RFMSHv1")) die ("bad reference header") ;
  U32 size ; if (fread (&size,sizeof(U32),1,f) != 1) die ("failed to read size") ;
  Reference *ref = referenceCreate (ms, size) ;
  if (fread (&ref->max,sizeof(U32),1,f) != 1) die ("failed to read max") ;
  if (fread (ref->index,sizeof(U32),size,f) != size) die ("failed read ref index") ;
  if (fread (ref->offset,sizeof(U32),size,f) != size) die ("failed read ref offset") ;
  if (fread (ref->id,sizeof(U32),size,f) != size) die ("failed read ref id") ;
  if (fread (ref->depth,sizeof(U32),ref->ms->max+1,f) != ref->ms->max+1) die ("fail depth") ;
  ref->rev = new (size, U32) ; ref->loc = new (ref->ms->max+1, U32) ;
  if (fread (ref->rev,sizeof(U32),size,f) != size) die ("fail rev") ;
  if (fread (ref->loc,sizeof(U32),ref->ms->max+1,f) != ref->ms->max+1) die ("fail loc") ;
  arrayDestroy (ref->len) ; if (!(ref->len = arrayRead (f))) die ("failed read ref len") ;
  dictDestroy (ref->dict) ; if (!(ref->dict = dictRead (f))) die ("failed read ref dict") ;
  fclose (f) ;
  return ref ;
}

/************************************************************/


void queryProcess (Reference *ref, char *filename)
{
  int len ;
  Array seeds = 0 ;

  dna2indexConv['N'] = dna2indexConv['n'] = 0 ; /* to get 2-bit encoding */
  SeqIO *si = seqIOopen (filename, dna2indexConv, FALSE) ;
  if (!si) die ("failed to read query sequence file %s", filename) ;
  while (seqIOread (si)) 
    { SeqhashRCiterator *mi = modRCiterator (ref->ms->hasher, sqioSeq(si), si->seqLen) ;
      U64 hash ; int pos ;
      seeds = arrayReCreate (seeds, 1024, Seed) ;
      int missed = 0, copy[4] ; copy[1] = copy[2] = copy[3] = 0 ;
      while (modRCnext (mi, &hash, &pos, 0))
	{ U32 index = modsetIndexFind (ref->ms, hash, FALSE) ;
	  Seed *s = arrayp(seeds,arrayMax(seeds),Seed) ;
	  s->index = index ; s->pos = pos ;
	  if (index) ++copy[msCopy(ref->ms,index)] ;
	  else ++missed ;
	}
      seqhashRCiteratorDestroy (mi) ;
      fprintf (outFile, "Q\t%s\t%d\t%d miss, %d copy1, %d copy2, %d multi, %.2f hit\n",
	       sqioId(si), si->seqLen, missed, copy[1], copy[2], copy[3],
	       (arrayMax(seeds)-missed)/(double)(arrayMax(seeds))) ;

      int i ;
      U32 loc0 = 0, locN, i0, iN ;
      int n1 = 0, n2 = 0 ;
      for (i = 0 ; i < arrayMax(seeds) ; ++i)
	{ Seed *s = arrp (seeds, i, Seed) ;
	  if (s->index && !msIsCopyM(ref->ms,s->index)) /* ignore multi-hits for now */
	    { U32 loc = ref->rev[ref->loc[s->index]] ;
	      BOOL is1 = msIsCopy1(ref->ms,s->index) ;
	      if (isVerbose)
		{ if (is1)
		    printf ("  %6d\t%s %d\n", s->pos,
			    dictName(ref->dict,ref->id[loc]), ref->offset[loc]) ;
		  else
		    { U32 loc2 = ref->rev[ref->loc[s->index]+1] ;
		      printf ("  %6d\t%s %d\t%s %d\n", s->pos,
			      dictName(ref->dict,ref->id[loc]), ref->offset[loc],
			      dictName(ref->dict,ref->id[loc2]), ref->offset[loc2]) ;
		    }
		}
	      BOOL endBlock = FALSE ;
	      if (!loc0 || ref->id[loc] != ref->id[loc0]) endBlock = TRUE ;
	      else if (loc0 < locN)
		{ if (loc < locN) endBlock = TRUE ;
		  int d = locN - loc0 - iN + i0 ; if (d > 50 || d < -50) endBlock = TRUE ;
		}
	      else if (loc0 > locN)
		{ if (loc > locN) endBlock = TRUE ;
		  int d = loc0 - locN - iN + i0 ; if (d > 50 || d < -50) endBlock = TRUE ;
		}
	      if (endBlock && loc0 && !is1) /* try the second loc */
		{ loc = ref->rev[ref->loc[s->index]+1] ;
		  endBlock = FALSE ;
		  if (ref->id[loc] != ref->id[loc0]) endBlock = TRUE ;
		   else if (loc0 < locN)
		     { if (loc < locN) endBlock = TRUE ;
		       int d = locN - loc0 - iN + i0 ; if (d > 50 || d < -50) endBlock = TRUE ;
		     }
		   else if (loc0 > locN)
		     { if (loc > locN) endBlock = TRUE ;
		       int d = loc0 - locN - iN + i0 ; if (d > 50 || d < -50) endBlock = TRUE ;
		     }
		}
	      if (endBlock)
		{ if (n1 > 2)
		    fprintf (outFile, "M\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d %d\t%.2f\t%.2f\n",
			     sqioId(si), arrp(seeds,i0,Seed)->pos, arrp(seeds,iN,Seed)->pos, len,
			     dictName(ref->dict,ref->id[loc0]),
			     ref->offset[loc0], ref->offset[locN],
			     n1, n2, (n1+n2) / (double)((locN>loc0) ? (locN-loc0) : (loc0-locN)),
			     n1 / (double)copy[1]) ;
		  n1 = 0 ; n2 = 0 ; loc0 = loc ; i0 = i ;
		}
	      if (is1) ++n1 ; else ++n2 ; locN = loc ; iN = i ;
	    }
	}
      if (n2 > 2)
	fprintf (outFile, "M\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d %d\t%.2f\t%.2f\n",
		 sqioId(si), arrp(seeds,i0,Seed)->pos, arrp(seeds,iN,Seed)->pos, len,
		 dictName(ref->dict,ref->id[loc0]),
		 ref->offset[loc0], ref->offset[locN],
		 n1, n2, (n1+n2) / (double)((locN>loc0) ? (locN-loc0) : (loc0-locN)),
		 n1 / (double)copy[1]) ;
    }
  seqIOclose (si) ;
  
  arrayDestroy (seeds) ;
}

/************************************************************/


/************* end of file ************/


// void usage (void)
// { fprintf (stderr, "Usage: modmap <commands>\n") ;
//   fprintf (stderr, "Commands are executed in order - set parameters before using them!\n") ;
//   fprintf (stderr, "  -K | --kmer <kmer size> [%d]\n", params.k) ;
//   fprintf (stderr, "  -W | --window <window> [%d]\n", params.w) ;
//   fprintf (stderr, "  -S | --seed <random number seed> [%d]\n", params.s) ;
//   fprintf (stderr, "  -B | --tableBits <hash index table bitcount> [%d]\n", params.B) ;
//   fprintf (stderr, "  -v | --verbose : toggle verbose mode\n") ;
//   fprintf (stderr, "  -t | --threads <number of threads for parallel ops> [%d]\n", numThreads) ;
//   fprintf (stderr, "  -o | --output <output filename> : '-' for stdout\n") ;
//   fprintf (stderr, "  -f | --referenceFasta <reference fasta file>\n") ;
//   fprintf (stderr, "  -w | --referenceWrite <file stem> : writes reference hash files\n") ;
//   fprintf (stderr, "  -r | --referenceRead <file stem> : read reference hash files\n") ;
//   fprintf (stderr, "  -q | --query <query fasta file>\n") ;
// }

// int main (int argc, char *argv[])
// {
//   --argc ; ++argv ;		/* eat program name */

//   outFile = stdout ;
  
//   timeUpdate (stdout) ;		/* initialise timer */
// #ifdef OMP
//   numThreads = omp_get_max_threads () ;
//   omp_set_num_threads (numThreads) ;
// #endif

//   /* command line parameters */
//   params.k = 19 ;
//   params.w = 31 ;
//   params.s = 17 ;
//   params.B = 28 ;

//   if (!argc) usage () ;

//   int i ;			/* generically useful variables */
//   FILE *f ;

//   Reference *ref = 0 ;

//   while (argc) {
//     if (**argv != '-')
//       die ("option/command %s does not start with '-': run without arguments for usage", *argv) ;
//     fprintf (stderr, "COMMAND %s", *argv) ;
//     for (i = 1 ; i < argc && *argv[i] != '-' ; ++i) fprintf (stderr, " %s", argv[i]) ;
//     fputc ('\n', stderr) ;
    
// #define ARGMATCH(x,y,n)	((!strcmp (*argv, x) || (!strcmp (*argv,y))) && argc >= n && (argc -= n, argv += n))
//     if (ARGMATCH("-K","--kmer",2)) params.k = atoi(argv[-1]) ;
//     else if (ARGMATCH("-W","--window",2)) params.w = atoi(argv[-1]) ;
//     else if (ARGMATCH("-S","--seed",2)) params.s = atoi(argv[-1]) ;
//     else if (ARGMATCH("-B","--tableBits",2)) params.B = atoi(argv[-1]) ;
//     else if (ARGMATCH("-t","--threads",2))
//       {
// #ifdef OMP
// 	numThreads = atoi(argv[-1]) ;
// 	if (numThreads > omp_get_max_threads ()) numThreads = omp_get_max_threads () ;
// 	omp_set_num_threads (numThreads) ;
// #else
// 	fprintf (stderr, "  can't set thread number - not compiled with OMP\n") ;
// #endif
//       }
//     else if (ARGMATCH("-v","--verbose",1)) isVerbose = !isVerbose ;
//     else if (ARGMATCH("-o","--output",2))
//       { if (!strcmp (argv[-1], "-"))
// 	  outFile = stdout ;
// 	else if (!(outFile = fopen (argv[-1], "w")))
// 	  { fprintf (stderr, "can't open output file %s - resetting to stdout\n", argv[-1]) ;
// 	    outFile = stdout ;
// 	  }
//       }
//     else if (ARGMATCH("-f","--referenceFasta",2))
//       { if (params.k <= 0 || params.w <= 0) die ("k %d, w %d must be > 0", params.k, params.w) ;
// 	Seqhash *hasher = seqhashCreate (params.k, params.w, params.s) ;
// 	fprintf (outFile, "  modmap initialised with k = %d, w = %d, random seed = %d\n",
// 		 params.k, params.w, params.s) ;
// 	Modset *ms = modsetCreate (hasher, params.B, 0) ;
// 	ref = referenceCreate (ms, 1 << 26) ;
// 	referenceFastaRead (ref, argv[-1], TRUE) ;
//       }
//     else if (ARGMATCH("-q","--query",2))
//       { if (!ref) die ("need to read a reference before processing query sequences") ;
// 	if (!(f = fopen (argv[-1], "r"))) die ("failed to open query file %s", argv[-1]) ;
// 	queryProcess (ref, argv[-1]) ;
//       }
//     else if (ARGMATCH("-r","--referenceRead",2))
//       { if (ref) referenceDestroy (ref) ;
// 	ref = referenceRead (argv[-1]) ;
//       }
//     else if (ARGMATCH("-w","--referenceWrite",2))
//       referenceWrite (ref, argv[-1]) ;
//     else die ("unkown command %s - run without arguments for usage", *argv) ;

//     timeUpdate (outFile) ;
//   }

//   fprintf (outFile, "total resources used: ") ; timeTotal (outFile) ;
//   if (outFile != stdout) { printf ("total resources used: ") ; timeTotal (stdout) ; }
// }
