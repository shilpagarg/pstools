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


Modset *modsetCreate (Seqhash *sh, int bits, U32 size) ;
void modsetDestroy (Modset *ms) ; 
void modsetWrite (Modset *ms, FILE *f) ;
Modset *modsetRead (FILE *f) ;

/* this is the key low level function, both to insert new hashes and find existing ones */
U32 modsetIndexFind (Modset *ms, U64 hash, int isAdd) ;

/* the following act on the whole set */
void modsetSummary (Modset *ms, FILE *f) ;
BOOL modsetPack (Modset *ms)	; /* reduce size to max+1 and compress value; TRUE if changes */
void modsetDepthPrune (Modset *ms, int min, int max) ;
BOOL modsetMerge (Modset *ms1, Modset *ms2) ;


Reference *referenceCreate (Modset *ms, U32 size);

void referenceDestroy (Reference *ref);

void referencePack (Reference *ref);

void referenceFastaRead (Reference *ref, char *filename, BOOL isAdd);

void referenceWrite (Reference *ref, char *root);

Reference *referenceRead (char *root);
/************************************************************/

typedef struct { U32 index ; U32 pos ; } Seed ;

void queryProcess (Reference *ref, char *filename);


static inline void msSetCopy0 (Modset *ms, U32 i) { ms->info[i] &= 0xfc ; }
static inline void msSetCopy1 (Modset *ms, U32 i) { ms->info[i] = (ms->info[i] & 0xfc) | 1 ; }
static inline void msSetCopy2 (Modset *ms, U32 i) { ms->info[i] = (ms->info[i] & 0xfc) | 2 ; }
static inline void msSetCopyM (Modset *ms, U32 i) { ms->info[i] |= 3 ; }
static inline BOOL msIsCopy0  (Modset *ms, U32 i) { return ((ms->info[i] & 3) == 0) ; }
static inline BOOL msIsCopy1  (Modset *ms, U32 i) { return ((ms->info[i] & 3) == 1) ; }
static inline BOOL msIsCopy2  (Modset *ms, U32 i) { return ((ms->info[i] & 3) == 2) ; }
static inline BOOL msIsCopyM  (Modset *ms, U32 i) { return ((ms->info[i] & 3) == 3) ; }
static inline int  msCopy (Modset *ms, U32 i) { return (ms->info[i] & 3) ; }


/************************************************************/

/************* end of file ************/
