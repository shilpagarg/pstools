/*  File: seqhash.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description: header file for seqhash package - minimizers and moshers
 * Exported functions: see below
 * HISTORY:
 * Last edited: Jan 26 22:27 2019 (rd109)
 * Created: Mon Mar  5 08:43:45 2018 (rd)
 *-------------------------------------------------------------------
 */

#include "utils.h"

typedef struct {
  int seed ;			/* seed ultimately */
  int k ;			/* kmer */
  int w ;			/* window */
  U64 mask ;			/* 2*k bits */
  int shift1, shift2 ;
  U64 factor1, factor2 ;
  U64 patternRC[4] ;		/* one per base */
} Seqhash ;

typedef struct {
  Seqhash *sh ;
  char *s, *sEnd ;     		/* sequence currently being hashed, end marker */
  U64 h, hRC ;			/* current hash values */
  U64 *hashBuf ;		/* buffer of length w holding hashes for current window */
  BOOL *fBuf ;			/* buffer of length w holding isForward for current window */
  int base ;			/* start of buf in sequence */
  int iStart, iMin ;		/* position in buf of start of current window, next min */
  BOOL isDone ;
} SeqhashRCiterator ;

Seqhash *seqhashCreate (int k, int w, int seed) ;
static void seqhashDestroy (Seqhash *sh) { free (sh) ; }

void seqhashWrite (Seqhash *sh, FILE *f) ;
Seqhash *seqhashRead (FILE *f) ;
void seqhashReport (Seqhash *sh, FILE *f) ;

// iterator to extract minimizers from a sequence
// NB sequence must continue to exist through the life of the iterator
SeqhashRCiterator *minimizerRCiterator (Seqhash *sh, char *s, int len) ;
BOOL minimizerRCnext (SeqhashRCiterator *si, U64 *u, int *pos, BOOL *isF) ; /* return u,pos,isF */

// modimizer extracts hashes that are divisible by m->w
// this is faster and more robust to errors - same mean density without evenness guarantees
SeqhashRCiterator *modRCiterator (Seqhash *sh, char *s, int len) ;
BOOL modRCnext (SeqhashRCiterator *si, U64 *u, int *pos, BOOL *isF) ; /* returns (u, pos, isF) */

static void seqhashRCiteratorDestroy (SeqhashRCiterator *si)
{ free (si->hashBuf) ; free (si->fBuf) ; free (si) ; }

/******* end of file ********/
