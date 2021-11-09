/*  File: seqio.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description: buffered package to read arbitrary sequence files - much faster than readseq
 * Exported functions:
 * HISTORY:
 * Last edited: Jan 30 15:41 2019 (rd109)
 * Created: Fri Nov  9 00:21:21 2018 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"

SeqIO *seqIOopen (char *filename, int* convert, BOOL isQual)
{
  SeqIO *si = new0 (1, SeqIO) ;
  if (!strcmp (filename, "-")) si->f = gzdopen (fileno (stdin), "r") ;
  else si->f = gzopen (filename, "r") ;
  if (!si->f) { free(si) ; return 0 ; }
  si->bufSize = 1<<24 ;
  si->b = si->buf = new (si->bufSize, char) ;
  si->convert = convert ;
  si->isQual = isQual ;
  si->nb = gzread (si->f, si->buf, si->bufSize) ;
  if (!si->nb)
    { fprintf (stderr, "sequence file %s unreadable or empty\n", filename) ;
      seqIOclose (si) ;
      return 0 ;
    }
  si->line = 1 ;
  if (*si->buf == '>')
    { si->type = FASTA ; si->isQual = FALSE ;
      if (!si->convert) si->convert = dna2textAmbigConv ; /* default: need to remove whitespace */
    }
  else if (*si->buf == '@') si->type = FASTQ ;
  else if (*si->buf == 'B')
    { si->type = BINARY ;
      if (si->nb <= 8) die ("binary file too short\n") ;
      si->qualThresh = si->buf[1] ;
      si->isQual = (si->qualThresh > 0) ;
      si->nSeq = *(U64*)si->buf & 0xffffffffffffLL ; /* 6 bytes of ff */
      si->b += 8 ; si->nb -= 8 ;
    }
  else
    { fprintf (stderr, "sequence file %s is unknown type\n", filename) ;
      seqIOclose (si) ;
      return 0 ;
    }
  return si ;
}

void seqIOclose (SeqIO *si)
{ free (si->buf) ;
  if (si->seqBuf) free (si->seqBuf) ; if (si->qualBuf) free (si->qualBuf) ;
  gzclose (si->f) ;
  free (si) ;
}

/********** local routines for seqIOread() ***********/
 
static int bufRefill (SeqIO *si)
{ si->b -= si->recStart ;		/* will be position after move */
  memmove (si->buf, si->buf + si->recStart, si->b - si->buf) ;
  si->idStart -= si->recStart ; si->descStart -= si->recStart ; /* adjust all the offsets */
  si->seqStart -= si->recStart ; si->qualStart -= si->recStart ;
  si->recStart = 0 ;
  si->nb = gzread (si->f, si->b, si->buf + si->bufSize - si->b) ;
  return si->nb ;
}

static int bufExtend (SeqIO *si)
{
  char *newbuf = new (si->bufSize*2, char) ;
  memcpy (newbuf, si->buf, si->bufSize) ;
  free (si->buf) ; si->buf = newbuf ;
  si->b = si->buf + si->bufSize ; si->nb = si->bufSize ;
  si->nb = gzread (si->f, si->b, si->bufSize) ;
  si->bufSize *= 2 ;
  return si->nb ;
}

#define bufAdvanceEndRecord(si) \
  { ++si->b ; \
    if (!--si->nb) \
      { if (si->recStart) bufRefill (si) ; \
        else if (si->b == si->buf + si->bufSize) bufExtend (si) ; \
      } \
  }

#define bufAdvanceInRecord(si) \
  { bufAdvanceEndRecord(si) ; \
    if (!si->nb) \
      { fprintf (stderr, "incomplete sequence record line %d\n", si->line) ; return FALSE ; } \
  } 

#define bufConfirmNbytes(si, n)	\
  { if (si->nb < n && si->recStart) bufRefill (si) ; \
    if (si->nb < n && si->b == si->buf + si->bufSize) bufExtend (si) ; \
    if (si->nb < n) \
      { fprintf (stderr, "incomplete sequence record line %d\n", si->line) ; return FALSE ; } \
  } 

#include <ctype.h>

BOOL seqIOread (SeqIO *si)
{
  if (!si->nb) return FALSE ;
  si->recStart = si->b - si->buf ;
  
  if (si->type == BINARY)
    { bufConfirmNbytes (si, 3*sizeof(int)) ;
      si->idLen = *((int*)si->b) ; si->b += sizeof(int) ;
      si->descLen = *((int*)si->b) ; si->b += sizeof(int) ;
      si->seqLen = *((int*)si->b) ; si->b += sizeof(int) ;
      si->nb -= 3*sizeof(int) ;
      int nBytes = si->idLen+1 + si->descLen+1
	+ (si->seqLen+3)/4 + si->isQual ? (si->seqLen+7)/8 : 0 ;
      nBytes = 4*((nBytes+3)/4) ;
      bufConfirmNbytes (si, nBytes) ;
      si->idStart = si->b - si->buf ;
      si->descStart = si->idStart + si->idLen+1 ;
      if (si->seqLen > si->seqBufSize)
	{ if (si->seqBufSize) { free (si->seqBuf) ; if (si->isQual) free (si->qualBuf) ; }
	  si->seqBufSize = 1 << 20 ; while (si->seqBufSize < si->seqLen) si->seqBufSize <<= 1 ;
	  si->seqBuf = new (si->seqBufSize, char) ;
	  if (si->isQual) si->qualBuf = new (si->seqBufSize, char) ;
	}
      si->seqStart = si->descStart + si->descLen+1 ;
      sqioSeqUnpack ((U8*)(si->buf+si->seqStart), si->seqBuf, si->seqLen) ;
      if (si->isQual)
	{ si->qualStart = si->seqStart + (si->seqLen+3)/4 ;
	  sqioQualUnpack ((U8*)(si->buf+si->qualStart), si->qualBuf, si->seqLen, si->qualThresh) ;
	}
      ++si->nSeq ;
      si->b += nBytes ; si->nb -= nBytes ;
      return TRUE ;
    }

  /* if get to here then this is a text file, FASTA or FASTQ */
  
  if (si->type == FASTA)
    { if (*si->b != '>') die ("no initial > for FASTA record line %d", si->line) ; }
  else if (si->type == FASTQ)
    { if (*si->b != '@') die ("no initial @ for FASTQ record line %d", si->line) ; }
  bufAdvanceInRecord(si) ; si->idStart = si->b - si->buf ;
  while (!isspace(*si->b)) bufAdvanceInRecord(si) ;
  si->idLen = si->b - sqioId(si) ;
  if (*si->b != '\n') /* a space or tab - whatever follows on this line is description */
    { *si->b = 0 ; bufAdvanceInRecord(si) ;
      si->descStart = si->b - si->buf ;
      while (*si->b != '\n') bufAdvanceInRecord(si) ;
      si->descLen = si->b - sqioDesc(si) ;
    }
  else { si->descLen = si->descStart = 0 ; }
  *si->b = 0 ;
  ++si->line ; bufAdvanceInRecord(si) ;	              /* line 2 */
  si->seqStart = si->b - si->buf ;
  if (si->type == FASTA)
    { while (si->nb && *si->b != '>')
	{ while (*si->b != '\n') bufAdvanceInRecord(si) ;
	  ++si->line ; bufAdvanceEndRecord(si) ;
	}
      char *s = sqioSeq(si), *t = s ;
      while (s < si->b) if ((*t++ = si->convert[*s++]) < 0) --t ;
      si->seqLen = t - sqioSeq(si) ;
    }
  else if (si->type == FASTQ)
    { while (*si->b != '\n') bufAdvanceInRecord(si) ;
      si->seqLen = si->b - sqioSeq(si) ;
      if (si->convert)
	{ char *s = sqioSeq(si) ;
	  while (s < si->b) { *s = si->convert[*s] ; ++s ; }
	}
      ++si->line ; bufAdvanceInRecord(si) ; 	      /* line 3 */
      if (*si->b != '+') die ("missing + FASTQ line %d", si->line) ;
      while (*si->b != '\n') bufAdvanceInRecord(si) ; /* ignore remainder of + line */
      ++si->line ; bufAdvanceInRecord(si) ;	      /* line 4 */
      si->qualStart = si->b - si->buf ;
      while (*si->b != '\n') bufAdvanceInRecord(si) ;
      if (si->b - si->buf - si->qualStart != si->seqLen)
	die ("qual not same length as seq line %d", si->line) ;
      if (si->isQual) { char *q = sqioQual(si), *e = q + si->seqLen ; while (q < e) *q++ -= 33 ; }
      ++si->line ; bufAdvanceEndRecord(si) ;
    }

  ++si->nSeq ;
  return TRUE ;
}

/* something for the future */

SeqIO *seqIObinaryOpen (char *filename, int nSeq) ; /* n = 0 if number unknown */
  
void seqIObinaryWrite (SeqIO *si) ;

/********** some routines to pack sequence and qualities **************/

void sqioSeqPack (char *s, U8 *u, int len) /* compress s into (len+3)/4 u  */
{
  int i ;
  while (len > 4)
    { *u = 0 ; for (i = 0 ; i < 4 ; ++i) *u = (*u << 2) | dna2indexConv[*s++] ;
      len -= 4 ; ++u ;
    }
  if (len) { *u = 0 ; for (i = 0 ; i < len ; ++i) *u = (*u << 2) | dna2indexConv[*s++] ; }
}

void sqioSeqUnpack (U8 *u, char *s, int len) /* uncompress (len+3)/4 u into s */
{
  int i ;
  while (len > 4)
    { for (i = 4 ; i-- ; ) { s[i] = (*u & 0x3) ; *u >>= 2 ; }
      len -= 4 ; ++u ; s += 4 ;
    }
  if (len) for (i = len ; --i ; ) { s[i] = (*u & 0x3) ; *u >>= 2 ; }
}

void sqioQualPack (char *q, U8 *u, int len, int thresh) /* compress q into (len+7)/8 u  */
{
  int i ;
  while (len > 8)
    { *u = 0 ; for (i = 0 ; i < 8 ; ++i) { if (*q++ >= thresh) *u |= 1 ; *u <<= 1 ; }
      len -= 8 ; ++u ;
    }
  if (len) { *u = 0 ; for (i = 0 ; i < len ; ++i) { if (*q++ >= thresh) *u |= 1 ; *u <<= 1 ; } }
}

void sqioQualUnpack (U8 *u, char *q, int len, int thresh) /* uncompress (len+7)/8 u into q */
{
  int i ;
  while (len > 8)
    { for (i = 8 ; i-- ; ) { q[i] = (*u & 0x1) ? thresh : 0 ; *u >>= 1 ; }
      len -= 8 ; ++u ; q += 8 ;
    }
  if (len) for (i = len ; --i ; ) { q[i] = (*u & 0x1) ? thresh : 0 ; *u >>= 1 ; }
}

/*********** standard conversion tables **************/

int dna2textConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2, 'A',  -2, 'C',  -2,  -2,  -2, 'G',  -2,  -2,  -2,  -2,  -2,  -2, 'N',  -2,
  -2,  -2,  -2,  -2, 'T',  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2, 'A',  -2, 'C',  -2,  -2,  -2, 'G',  -2,  -2,  -2,  -2,  -2,  -2, 'N',  -2,
  -2,  -2,  -2,  -2, 'T',  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2textAmbigConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, '-',  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2, 'A', 'B', 'C', 'D',  -2,  -2, 'G', 'H',  -2,  -2, 'K',  -2, 'M', 'N',  -2,
  -2,  -2, 'R', 'S', 'T',  -2, 'V', 'W',  -2, 'Y',  -2,  -2,  -2,  -2,  -2,  -2,
  -2, 'A', 'B', 'C', 'D',  -2,  -2, 'G', 'H',  -2,  -2, 'K',  -2, 'M', 'N',  -2,
  -2,  -2, 'R', 'S', 'T',  -2, 'V', 'W',  -2, 'Y',  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2textAmbig2NConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2, 'A', 'N', 'C', 'N',  -2,  -2, 'G', 'N',  -2,  -2, 'N',  -2, 'N', 'N',  -2,
  -2,  -2, 'N', 'N', 'T',  -2, 'N', 'N',  -2, 'N',  -2,  -2,  -2,  -2,  -2,  -2,
  -2, 'A', 'N', 'C', 'N',  -2,  -2, 'G', 'N',  -2,  -2, 'N',  -2, 'N', 'N',  -2,
  -2,  -2, 'N', 'N', 'T',  -2, 'N', 'N',  -2, 'N',  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2indexConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   0,  -2,   1,  -2,  -2,  -2,   2,  -2,  -2,  -2,  -2,  -2,  -2,   4,  -2,
  -2,  -2,  -2,  -2,   3,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2,   0,  -2,   1,  -2,  -2,  -2,   2,  -2,  -2,  -2,  -2,  -2,  -2,   4,  -2,
  -2,  -2,  -2,  -2,   3,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2binaryConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   1,  -2,   2,  -2,  -2,  -2,   4,  -2,  -2,  -2,  -2,  -2,  -2,  15,  -2,
  -2,  -2,  -2,  -2,   8,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2,   1,  -2,   2,  -2,  -2,  -2,   4,  -2,  -2,  -2,  -2,  -2,  -2,  15,  -2,
  -2,  -2,  -2,  -2,   8,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2binaryAmbigConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,   0,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   1,  14,   2,  13,  -2,  -2,   4,  11,  -2,  -2,  12,  -2,   3,  15,  -2,
  -2,  -2,   5,   6,   8,  -2,   7,   9,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2,   1,  14,   2,  13,  -2,  -2,   4,  11,  -2,  -2,  12,  -2,   3,  15,  -2,
  -2,  -2,   5,   6,   8,  -2,   7,   9,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int aa2textConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2, 'A', 'X', 'C', 'D', 'E', 'F', 'G', 'H', 'I',  -2, 'K', 'L', 'M', 'N',  -2,
 'P', 'Q', 'R', 'S', 'T',  -2, 'V', 'W', 'X', 'Y', 'X',  -2,  -2,  -2,  -2,  -2,
  -2, 'A', 'X', 'C', 'D', 'E', 'F', 'G', 'H', 'I',  -2, 'K', 'L', 'M', 'N',  -2,
 'P', 'Q', 'R', 'S', 'T',  -2, 'V', 'W', 'X', 'Y', 'X',  -2,  -2,  -2,  -2,  -2,
} ;

int aa2indexConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   0,  20,   1,   2,   3,   4,   5,   6,   7,  -2,   8,   9,  10,  11,  -2,
  12,  13,  14,  15,  16,  -2,  17,  18,  20,  19,  20,  -2,  -2,  -2,  -2,  -2,
  -2,   0,  20,   1,   2,   3,   4,   5,   6,   7,  -2,   8,   9,  10,  11,  -2,
  12,  13,  14,  15,  16,  -2,  17,  18,  20,  19,  20,  -2,  -2,  -2,  -2,  -2,
} ;

int noConv[] = {
   0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,
  96,  97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127
} ;

/******** end of file ********/
