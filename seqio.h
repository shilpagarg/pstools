/*  File: seqio.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 27 11:21 2018 (rd109)
 * Created: Sat Nov 10 08:51:49 2018 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include <zlib.h>

typedef enum { FASTA, FASTQ, BINARY } SeqIOtype ;

typedef struct {
  SeqIOtype type ;
  U64 nSeq ;			/* number of sequences */
  int idLen, descLen, seqLen ;
  int idStart, descStart, seqStart, qualStart ;
  BOOL isQual ;			/* if set then convert qualities by subtracting 33 (FASTQ) */
  int qualThresh ;		/* used for binary representation of qualities */
  /* below here private */
  int bufSize ;
  int nb ;			/* nb is how many characters left to read in the buffer */
  int line, recStart ;		/* recStart is the offset for the current record */
  gzFile f ;
  char *buf, *b ;		/* b is current pointer in buf */
  int *convert ;
  int  seqBufSize ;
  char *seqBuf, *qualBuf ;	/* only used in binary mode, because buf then holds compressed forms */
} SeqIO ;

/* Reads FASTA or FASTQ, gzipped or not. */
/* Philosophy here is to read blocks of 8Mb and provide direct access into the buffer. */
/* So the user does not own the pointers. */
/* Add 0 terminators to ids.  Convert sequences in place if convert != 0, and quals if isQual. */
/* Potential for future storage of (packed?) sequences with counts up front. */ 

SeqIO *seqIOopen (char *filename, int* convert, BOOL isQual) ; /* can use "-" for stdin */
BOOL seqIOread (SeqIO *si) ;
void seqIOclose (SeqIO *si) ;
#define sqioId(si)   ((si)->buf+(si)->idStart)
#define sqioDesc(si) ((si)->buf+(si)->descStart)
#define sqioSeq(si)  ((si)->type == BINARY ? (si)->seqBuf : (si)->buf+(si)->seqStart)
#define sqioQual(si) ((si)->type == BINARY ? (si)->qualBuf : (si)->buf+(si)->qualStart)

void sqioSeqPack (char *s, U8 *u, int len) ; /* compress s into (len+3)/4 u */
void sqioSeqUnpack (U8 *u, char *s, int len) ; /* uncompress (len+3)/4 u into s */
void sqioQualPack (char *q, U8 *u, int len, int thresh) ; /* compress q into (len+7)/8 u */
void sqioQualUnpack (U8 *u, char *q, int len, int thresh) ; /* uncompress (len+7)/8 u into q */

extern int dna2textConv[] ;
extern int dna2textAmbigConv[] ;
extern int dna2textAmbig2NConv[] ;
extern int dna2indexConv[] ;
extern int dna2binaryConv[] ;
extern int dna2binaryAmbigConv[] ;
static const char index2char[] = "acgtn" ;
static const char binary2char[] = "-ACMGRSVTWYHKDBN" ;
extern int aa2textConv[] ;
extern int aa2indexConv[] ;
static const char index2aa[] = "ACDEFGHIKLMNPQRSTVWYX*" ;
extern int noConv[] ;

/******************************************************************/
