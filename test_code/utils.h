/*  File: utils.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2011
 *-------------------------------------------------------------------
 * Description: includes standard system headers and own headers
 * Exported functions:
 * HISTORY:
 * Last edited: Jan 12 13:09 2019 (rd109)
 * Created: Wed Jan  5 16:13:48 2011 (rd)
 *-------------------------------------------------------------------
 */

#include <stdio.h>		/* FILE etc. */
#include <stdlib.h>		/* malloc(), free(), ... notation */
#include <stdint.h>		/* for standard size int types */
#include <string.h>		/* memset() */
#include <limits.h>		/* INT_MAX etc. */
#include <assert.h>

#ifndef UTILS_DEFINED
#define UTILS_DEFINED
typedef char BOOL ;
#define TRUE 1
#define FALSE 0

typedef int8_t I8 ;
const static I8 I8MAX = 0x7f ;
typedef int16_t I16 ;
const static I16 I16MAX = 0x7fff ;
typedef int32_t I32 ;
const static I32 I32MAX = 0x7fffffff ;
typedef int64_t I64 ;
const static I64 I64MAX = 0x7fffffffffffffff ;

typedef uint8_t U8 ;
const static U8 U8MAX = 0xff ;
typedef uint16_t U16 ;
const static U16 U16MAX = 0xffff ;
typedef uint32_t U32 ;
const static U32 U32MAX = 0xffffffff ;
typedef uint64_t U64 ;
const static U64 U64MAX = 0xffffffffffffffff ;
#endif

#include "array.h"
#include "dict.h"
#include "hash.h"

void die (char *format, ...) ;
void *myalloc (size_t size) ;
void *mycalloc (size_t number, size_t size) ;
#define	new(n,type)	(type*)myalloc((n)*sizeof(type))
#define	new0(n,type)	(type*)mycalloc((n),sizeof(type))
#define resize(x,nOld,nNew,T) { T* z = new((nNew),T) ; if (nOld < nNew) memcpy(z,x,(nOld)*sizeof(T)) ; else memcpy(z,x,(nNew)*sizeof(T)) ; free(x) ; x = z ; }

char *fgetword (FILE *f) ;	/* not threadsafe */
FILE *fopenTag (char* root, char* tag, char* mode) ;
void timeUpdate (FILE *f) ;	/* print time usage since last call to file */
void timeTotal (FILE *f) ;	/* print full time usage since first call to timeUpdate */

/************************/
