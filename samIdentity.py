#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("samfile", help="input sam or bam file" )
#parser.add_argument("",help="",default=None)
parser.add_argument('-d', action="store_true", default=False)
parser.add_argument("--header", action="store_true", default=False)
parser.add_argument("--tag", help="ID tag to add to each row", default=None)
parser.add_argument("--mask", help="a list of query names to add a mask flag == True for", default=None)
args = parser.parse_args()
DEBUG=args.d

import re
import os
import sys
import numpy as np
import pysam 
import pandas as pd

print("Modules Loaded", file=sys.stderr)

def makeHeader():
	rtn = ["query_name", "query_start", "query_end", "query_length",
			"reference_name", "reference_start", "reference_end", "reference_length",
			"perID_by_matches", "perID_by_events", "perID_by_all", 
			"matches", "mismatches", 
			"insertions", "deletions", "insertion_events", "deletion_events"]

	return(rtn)

def perId(matches, mismatch, ins, dele, insEvent, delEvent):
	if( (matches + mismatch) == 0):
		return((0,0,0))

	bymatches = (100.0 * matches)/(matches + mismatch)
	byevents = (100.0 * matches)/(matches + mismatch + insEvent + delEvent)
	byall = (100.0 * matches)/(matches + mismatch + ins + dele)
	return((bymatches, byevents, byall))


# /(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)+/
# :[0-9]+ represents an identical block,
# -ata represents a deltion,i
# +gtc an insertion and
# *at indicates reference base a is substituted with a query base t.
# It is similar to the MD SAM tag but is standalone and easier to parse.

def formatRead(read, ref_length):
	
	if(read.flag == 4):
		return(None)
	
	tags = read.get_tags()
	if(len(tags) > 0):
		tags, values = zip(*tags)
	if( 'cs' in tags):
		cs = read.get_tag('cs')
		pattern = ":[0-9]+"
		match = re.findall(pattern, cs)
		if(match is not None):
			match = [ int(x[1:]) for x in match ]
			matches = sum(match)
		
		pattern = "-[a-z]+"
		match = re.findall(pattern, cs)
		if(match is not None):
			#print(match)
			match = [ len((x[1:])) for x in match ]
			delEvent = len(match)
			dele = sum(match)

		pattern = "\+[a-z]+"
		match = re.findall(pattern, cs)
		if(match is not None):
			#print(match)
			match = [ len((x[1:])) for x in match ]
			insEvent = len(match)
			ins = sum(match)

		mismatch = cs.count("*")
		#print(cs)
	else:
		counts, events = read.get_cigar_stats()
		matches = counts[7]
		mismatch = counts[8]
		ins = counts[1]
		dele = counts[2]
		insEvent = events[1]
		delEvent = events[2]
	
	bymatches, byevents, byall = perId(matches, mismatch, ins, dele, insEvent, delEvent)

	rtn = [read.query_name, read.query_alignment_start, read.query_alignment_end, read.infer_query_length(),
			read.reference_name, read.reference_start, read.reference_end, ref_length,
			bymatches, byevents, byall,
			matches, mismatch, ins, dele, insEvent, delEvent
			]
	return(rtn)





out = []
samfile = pysam.AlignmentFile(args.samfile)
for read in samfile.fetch(until_eof=True):
	ref_length = samfile.get_reference_length(read.reference_name)
	add = formatRead(read, ref_length)
	if(add is not None):
		out.append(add)

samfile.close()

df = pd.DataFrame.from_records(out)
df.columns = makeHeader()

if(args.tag is not None):
	df["id"] = args.tag

if(args.mask is not None):
	tomask = set([ name.strip() for name in open(args.mask).readlines() ])
	df["mask"] = False
	for idx, row in df.iterrows():
		if(row["query_name"] in tomask):
			df["mask"].iat[idx] = True

df.to_csv(sys.stdout, sep ="\t", header=args.header, index=False)