#!/usr/bin/env python3
"""Map degenerate primers to FASTQ/FASTQ.gz; one row per matching location."""
import argparse,csv
from primer_common import fastq,read_primers,matches
p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--primers',required=True);p.add_argument('--output',required=True);p.add_argument('--max-mismatches',type=int,default=0);a=p.parse_args()
if a.max_mismatches<0:p.error('mismatches must be nonnegative')
primers=list(read_primers(a.primers))
with open(a.output,'w',newline='') as out:
 w=csv.writer(out,delimiter='\t');w.writerow(['read_id','primer_id','orientation','start_0based','end_0based_exclusive','mismatches','read_length'])
 for ident,seq in fastq(a.input):
  for name,primer in primers:
   for orient,start,mm in matches(seq,primer,a.max_mismatches):w.writerow([ident,name,orient,start,start+len(primer),mm,len(seq)])
