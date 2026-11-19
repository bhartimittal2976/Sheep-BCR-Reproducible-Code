#!/usr/bin/env python3
"""Dependency-free IUPAC primer matching (both orientations, approximate mismatches)."""
import csv,gzip
from pathlib import Path
IUPAC={'A':set('A'),'C':set('C'),'G':set('G'),'T':set('T'),'R':set('AG'),'Y':set('CT'),'S':set('GC'),'W':set('AT'),'K':set('GT'),'M':set('AC'),'B':set('CGT'),'D':set('AGT'),'H':set('ACT'),'V':set('ACG'),'N':set('ACGT')}
COMP=str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')
def revcomp(s):return s.upper().translate(COMP)[::-1]
def read_primers(path):
 with open(path,newline='') as f:
  for r in csv.DictReader(f):
   s=r['sequence'].strip().upper().replace('U','T'); assert s and all(x in IUPAC for x in s),f"Invalid primer {r}"
   yield r['primer_id'],s
def open_text(path):return gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)
def fasta(path):
 with open_text(path) as f:
  ident=None; parts=[]
  for line in f:
   line=line.strip()
   if line.startswith('>'):
    if ident is not None:yield ident,''.join(parts).upper()
    ident=line[1:].split()[0];parts=[]
   elif line:parts.append(line)
  if ident is not None:yield ident,''.join(parts).upper()
def fastq(path):
 with open_text(path) as f:
  while True:
   head=f.readline()
   if not head:break
   seq=f.readline().strip();plus=f.readline();qual=f.readline().strip()
   if not head.startswith('@') or not plus.startswith('+') or len(seq)!=len(qual):raise ValueError('Malformed FASTQ')
   yield head[1:].split()[0],seq.upper()
def matches(seq,primer,max_mm):
 for orient,p in [('+',primer),('-',revcomp(primer))]:
  for start in range(len(seq)-len(p)+1):
   mismatches=sum(base not in IUPAC[code] for base,code in zip(seq[start:start+len(p)],p))
   if mismatches<=max_mm:yield orient,start,mismatches
