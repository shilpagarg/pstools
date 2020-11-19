/*  File: dict.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2011
 *-------------------------------------------------------------------
 * Description: based on acedb code from Jean Thierry-Mieg and Richard Durbin 1999-2004
 * -------------------------------------------------------------------
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jan 24 07:06 2019 (rd109)
 * Created: July 2003 (rd)
 *-------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

/****************************************/

static void* remap (void *old, int oldSize, int newSize)
{
  void* new_buf = mycalloc (newSize, 1) ;
  memcpy(new_buf, old, oldSize) ;
  free (old) ;
  return new_buf ;
}

/****************************************/

static int hashString (char *cp, int n, int isDiff)
{
  int i ;
  unsigned int j, x = 0 ;
  int rotate = isDiff ? 21 : 13 ;
  int leftover = 8 * sizeof(int) - rotate ;

  while (*cp)
    x = (*cp++) ^ ((x >> leftover) | (x << rotate)) ;

  for (j = x, i = n ; i < sizeof(int) ; i += n)
    j ^= (x >> i) ;
  j &= (1 << n) - 1 ;

  if (isDiff)
    j |= 1 ;

  return j ;
}

/*****************************/

DICT *dictCreate (int size)
{
  DICT *dict = (DICT*) mycalloc (1, sizeof(DICT)) ;

  for (dict->dim = 10, dict->size = 1024 ; dict->size < size ; ++dict->dim, dict->size *= 2) ;
  dict->table = (int*) mycalloc (dict->size, sizeof(int)) ;
  dict->names = (char**) mycalloc (dict->size/2, sizeof(char*)) ;
  return dict ; 
}

/*****************************/

void dictDestroy (DICT *dict)
{
  int i ;
  for (i = 1 ; i <= dict->max ; ++i) free (dict->names[i]) ;
  free (dict->names) ;
  free (dict->table) ;
  free (dict) ;
}

/*****************************/

BOOL dictWrite (DICT *dict, FILE *f)
{
  if (fwrite (&dict->dim,sizeof(int),1,f) != 1) return FALSE ;
  if (fwrite (&dict->max,sizeof(int),1,f) != 1) return FALSE ;
  if (fwrite (dict->table,sizeof(int),dict->size,f) != dict->size) return FALSE ;
  if (fwrite (dict->names,sizeof(char*),dict->max+1,f) != dict->max+1) return FALSE ;
  int i ;
  for (i = 1 ; i <= dict->max ; ++i)
    { int len = strlen(dict->names[i]) ;
      if (fwrite (&len,sizeof(int),1,f) != 1) return FALSE ;
      if (fwrite (dict->names[i],1,len,f) != len) return FALSE ;
    }
  return TRUE ;
}
  
DICT *dictRead (FILE *f)
{
  int dim ; if (fread (&dim,sizeof(int),1,f) != 1) return 0 ;
  DICT *dict = dictCreate (1 << dim) ;
  if (fread (&dict->max,sizeof(int),1,f) != 1) return 0 ;
  if (fread (dict->table,sizeof(int),dict->size,f) != dict->size) return 0 ;
  if (fread (dict->names,sizeof(char*),dict->max+1,f) != dict->max+1) return 0 ;
  int i ;
  for (i = 1 ; i <= dict->max ; ++i)
    { int len ;
      if (fread (&len,sizeof(int),1,f) != 1) return 0 ;
      dict->names[i] = new0 (len+1, char) ;
      if (fread (dict->names[i],1,len,f) != len) return 0 ;
    }
  return dict ;
}

/*****************************/

static int newPos ;		/* communication between dictFind() and dictAdd() */

BOOL dictFind (DICT *dict, char *s, int *ip)
{
  int i, x, d ;

  if (!dict) die ("dictAdd received null dict\n") ;
  if (!s) die ("dictAdd received null string\n") ;

  x = hashString (s, dict->dim, 0) ;
  if (!(i = dict->table[x]))
    { newPos = x ; 
      return FALSE ; 
    }
  else if (!strcmp (s, dict->names[i]))
    { if (ip) *ip = i-1 ; 
      return TRUE ; 
    }
  else
    { d = hashString (s, dict->dim, 1) ;
      while (1)
	{ x = (x + d) & ((1 << dict->dim) - 1) ;
	  if (!(i = dict->table[x]))
	    { newPos = x ; 
	      return FALSE ; 
	    }
	  else if (!strcmp (s, dict->names[i]))
	    { if (ip) *ip = i-1 ; 
	      return TRUE ; 
	    }
	}
    }
}

/*****************************/

BOOL dictAdd (DICT *dict, char *s, int *ip)
{
  int i, x ;

  if (dictFind (dict, s, ip)) return FALSE ;

  i = ++dict->max ;
  dict->table[newPos] = i ;
  dict->names[i] = (char*) myalloc (strlen(s) + 1) ;
  strcpy (dict->names[i], s) ;
  if (ip) *ip = i-1 ;

  if (dict->max > 0.3 * dict->size) /* double table size and remap */
    { int *newTable ;
      ++dict->dim ; dict->size *= 2 ;
      dict->names = (char**) remap (dict->names, (dict->max+1)*sizeof(char*), (dict->size/2)*sizeof(char*)) ;
      newTable = (int*) mycalloc (dict->size, sizeof(int)) ;
      for (i = 1 ; i <= dict->max ; ++i)
	{ s = dict->names[i] ;
	  x = hashString (s, dict->dim, 0) ;
	  if (!newTable[x])
	    newTable[x] = i ;
	  else
	    { int d = hashString (s, dict->dim, 1) ;
	      while (1)
		{ x = (x + d) & ((1 << dict->dim) - 1) ;
		  if (!newTable[x])
		    { newTable[x] = i ; break ; }
		}
	    }
	}
      free (dict->table) ; dict->table = newTable ;
    }

  return TRUE ;
}

/*****************************/

char* dictName (DICT *dict, int i)
{ return dict->names[i+1] ; }

/*********** end of file ***********/
