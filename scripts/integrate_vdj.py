#!/usr/bin/env python3
"""Integrate Cell Ranger filtered contig annotations by sample and barcode.
Standard library only; no biological reference or private input files bundled.
"""
import argparse,collections,csv,itertools,json
from pathlib import Path
REGIONS=['fwr1','cdr1','fwr2','cdr2','fwr3','cdr3','fwr4']
CHAIN_TYPES={'IGH','IGK','IGL'}
ISOTYPES={'IGHG':'IgG','IGHM':'IgM','IGHA':'IgA','IGHD':'IgD','IGHE':'IgE'}
ORDER=['IgG only','IgM only','IgA only','IgD only','IgE only','Other only','IgL only','IgK only']+[i+'/'+l for i in ['IgG','IgM','IgA','IgD','IgE','Other'] for l in ['IgL','IgK']]+['multiple chains','no productive chains']
CODON={'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W','CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R','ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
def iso(c):
 c=(c or '').upper();return next((v for k,v in ISOTYPES.items() if c.startswith(k)),'Other')
def reconstructed(r):
 nt=[(r.get(k+'_nt') or '').upper() for k in REGIONS]
 aa=[r.get(k) or '' for k in REGIONS]
 if not all(nt) or not all(aa):return '', ''.join(aa), 'INCOMPLETE_REGION'
 # Cell Ranger sometimes includes an extra terminal nucleotide in FWR4.
 if len(nt[-1])==3*len(aa[-1])+1:nt[-1]=nt[-1][:-1]
 sequence=''.join(nt);protein=''.join(aa)
 if len(sequence)%3:return sequence,protein,'NOT_CODON_ALIGNED'
 translated=''.join(CODON.get(sequence[i:i+3],'X') for i in range(0,len(sequence),3))
 return sequence,protein,'PASS' if translated==protein else 'TRANSLATION_MISMATCH'
def key(r):
 nt,aa,status=reconstructed(r)
 if nt and status=='PASS':return (r['chain'],nt,r.get('c_gene',''))
 return (r['chain'],'INCOMPLETE_OR_QC_FAIL',r.get('cdr3_nt',''),r.get('v_gene',''),r.get('j_gene',''),r.get('c_gene',''),r.get('contig_id',''))
def category(chains):
 H=[r for r in chains if r['chain']=='IGH'];L=[r for r in chains if r['chain'] in ('IGK','IGL')]
 if len(H)>1 or len(L)>1:return 'multiple chains'
 if H and L:return iso(H[0].get('c_gene'))+'/'+('IgK' if L[0]['chain']=='IGK' else 'IgL')
 if H:return iso(H[0].get('c_gene'))+' only'
 if L:return 'IgK only' if L[0]['chain']=='IGK' else 'IgL only'
 return 'no productive chains'
def write(path,rows,fields):
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',required=True);p.add_argument('--outdir',required=True);a=p.parse_args()
 manifest=Path(a.manifest).resolve();out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
 with manifest.open(newline='') as f: entries=list(csv.DictReader(f))
 assert entries and all(set(['run','sample','stage','contigs'])<=set(x) for x in entries)
 runs=set();data=collections.defaultdict(list);observations=collections.defaultdict(lambda:collections.Counter())
 for e in entries:
  run,sample,stage=e['run'],e['sample'],e['stage'].lower()
  if run in runs:raise ValueError('Duplicate run '+run)
  runs.add(run)
  if stage not in ('seq1','seq2'):raise ValueError('stage must be seq1 or seq2')
  path=Path(e['contigs']);path=path if path.is_absolute() else manifest.parent.parent/path
  with path.open(newline='') as f:
   for r in csv.DictReader(f):
    if r.get('is_cell','').lower()!='true' or r.get('productive','').lower()!='true' or r.get('chain') not in CHAIN_TYPES:continue
    bc=r['barcode'];r['_run']=run;r['_sample']=sample;r['_stage']=stage
    data[(sample,stage,bc)].append(r)
  # Library-level barcode observations computed from productive eligible records only.
  for (s,st,bc),rows in data.items():
   if s==sample and st==stage and any(x['_run']==run for x in rows):observations[(sample,stage,run)]['_barcodes']+=0
 def analyze(stage):
  cells=[];chains_out=[];pairs=[];summary=collections.defaultdict(collections.Counter);run_bcs=collections.defaultdict(set)
  grouped=collections.defaultdict(list)
  for (sample,st,bc),rows in data.items():
   if stage=='seq2' and st!='seq2':continue
   if stage=='seq1' and st!='seq1':continue
   grouped[(sample,bc)].extend(rows)
   for r in rows:run_bcs[(sample,r['_run'])].add(bc)
  for (sample,bc),rows in sorted(grouped.items()):
   bykey=collections.defaultdict(list)
   for r in rows:bykey[key(r)].append(r)
   unique=[]
   for k,rs in bykey.items():
    r=max(rs,key=lambda x:(reconstructed(x)[2]=='PASS',int(x.get('umis') or 0),int(x.get('reads') or 0)))
    nt,aa,qc=reconstructed(r)
    sources=';'.join(sorted(set(x['_run'] for x in rs)))
    chainrow={'sample':sample,'barcode':bc,'cell_id':sample+'|'+bc,'chain':r['chain'],'isotype':iso(r.get('c_gene')) if r['chain']=='IGH' else '', 'c_gene':r.get('c_gene',''),'v_gene':r.get('v_gene',''),'d_gene':r.get('d_gene',''),'j_gene':r.get('j_gene',''),'cdr3_aa':r.get('cdr3',''),'cdr3_nt':r.get('cdr3_nt',''),'source_runs':sources,'source_contigs':';'.join(sorted(x['_run']+':'+x.get('contig_id','') for x in rs)),'variable_nt':nt,'variable_aa':aa,'translation_qc':qc}
    for region in REGIONS:chainrow[region+'_aa']=r.get(region,'');chainrow[region+'_nt']=r.get(region+'_nt','')
    unique.append(chainrow);chains_out.append(chainrow)
   cat=category(unique);H=[r for r in unique if r['chain']=='IGH'];L=[r for r in unique if r['chain']!='IGH']
   summary[sample][cat]+=1;summary[sample]['total barcodes']+=1
   if len(H)==len(L)==1:summary[sample]['strict 1H1L']+=1
   if H and L:summary[sample]['any H+L']+=1
   cells.append({'sample':sample,'barcode':bc,'category':cat,'n_heavy':len(H),'n_light':len(L),'n_IGK':sum(x['chain']=='IGK' for x in L),'n_IGL':sum(x['chain']=='IGL' for x in L),'source_runs':';'.join(sorted(set(x['_run'] for x in rows)))})
   for h,l in itertools.product(H,L):
    pairs.append({'sample':sample,'barcode':bc,'status':'STRICT_1H1L' if len(H)==len(L)==1 else 'MULTICHAIN_CANDIDATE','heavy_isotype':h['isotype'],'light_chain':l['chain'],**{'heavy_'+k:h.get(k,'') for k in h if k not in ('sample','barcode','cell_id')},**{'light_'+k:l.get(k,'') for k in l if k not in ('sample','barcode','cell_id')}})
  prefix=stage
  write(out/f'{prefix}_cells.tsv',cells,['sample','barcode','category','n_heavy','n_light','n_IGK','n_IGL','source_runs'])
  cf=['sample','barcode','cell_id','chain','isotype','c_gene','v_gene','d_gene','j_gene','cdr3_aa','cdr3_nt','source_runs','source_contigs','variable_nt','variable_aa','translation_qc']+[x+'_'+y for x in REGIONS for y in ('aa','nt')]
  write(out/f'{prefix}_unique_chains.tsv',chains_out,cf)
  pf=['sample','barcode','status','heavy_isotype','light_chain']+['heavy_'+x for x in cf if x not in ('sample','barcode','cell_id')]+['light_'+x for x in cf if x not in ('sample','barcode','cell_id')]
  write(out/f'{prefix}_all_candidate_pairs.tsv',pairs,pf)
  write(out/f'{prefix}_strict_1H1L_pairs.tsv',(r for r in pairs if r['status']=='STRICT_1H1L'),pf)
  return summary,run_bcs,cells
 results={};bcsets={};cellstages={}
 for stage in ['seq1','seq2','merged']:
  results[stage],bcsets[stage],cellstages[stage]=analyze(stage)
 # Precompute individual-library counts once (library-level observations).
 per_run=collections.defaultdict(collections.Counter)
 for (sample,st,bc),rs in data.items():
  byrun=collections.defaultdict(list)
  for r in rs:byrun[r['_run']].append(r)
  for run,rr in byrun.items():
   distinct={key(r):r for r in rr};chains=list(distinct.values());c=category(chains)
   per_run[(sample,run)][c]+=1;per_run[(sample,run)]['total barcodes']+=1
   H=[r for r in chains if r['chain']=='IGH'];L=[r for r in chains if r['chain']!='IGH']
   if len(H)==len(L)==1:per_run[(sample,run)]['strict 1H1L']+=1
   if H and L:per_run[(sample,run)]['any H+L']+=1
 summary=[]
 samples=sorted(set(e['sample'] for e in entries))
 for stage in ['seq1','seq2','merged']:
  for category_name in ['total barcodes']+ORDER+['strict 1H1L','any H+L']:
   row={'stage':stage,'category':category_name}
   for sample in samples:
    row[sample+'_unique']=results[stage][sample][category_name]
    relevant=[e['run'] for e in entries if e['sample']==sample and (stage=='merged' or e['stage'].lower()==stage)]
    row[sample+'_observations']=sum(per_run[(sample,run)][category_name] for run in relevant)
   row['total_observations']=sum(row[s+'_observations'] for s in samples)
   row['total_unique']=sum(row[s+'_unique'] for s in samples)
   summary.append(row)
 write(out/'summary_observed_vs_unique.tsv',summary,['stage','category']+[k for s in samples for k in (s+'_observations',s+'_unique')]+['total_observations','total_unique'])
 overlaps=[]
 for sample in samples:
  seq2runs=[e['run'] for e in entries if e['sample']==sample and e['stage'].lower()=='seq2']
  if len(seq2runs)==2:
   a,b=seq2runs;A=bcsets['seq2'][(sample,a)];B=bcsets['seq2'][(sample,b)]
   overlaps.append({'sample':sample,'run1':a,'run2':b,'run1_barcodes':len(A),'run2_barcodes':len(B),'shared':len(A&B),'unique_to_run1':len(A-B),'unique_to_run2':len(B-A),'union':len(A|B)})
 write(out/'seq2_barcode_overlap.tsv',overlaps,['sample','run1','run2','run1_barcodes','run2_barcodes','shared','unique_to_run1','unique_to_run2','union'])
 change=[]
 for sample in samples:
  d={st:{r['barcode']:r for r in cellstages[st] if r['sample']==sample} for st in ['seq1','seq2','merged']}
  for bc in sorted(d['merged']):
   change.append({'sample':sample,'barcode':bc,'seq1_category':d['seq1'].get(bc,{}).get('category','ABSENT'),'seq2_category':d['seq2'].get(bc,{}).get('category','ABSENT'),'merged_category':d['merged'][bc]['category']})
 write(out/'barcode_category_changes.tsv',change,['sample','barcode','seq1_category','seq2_category','merged_category'])
 print(json.dumps({'output':str(out),'samples':samples,'merged_unique':sum(results['merged'][s]['total barcodes'] for s in samples),'merged_strict_pairs':sum(results['merged'][s]['strict 1H1L'] for s in samples)},indent=2))
if __name__=='__main__':main()
