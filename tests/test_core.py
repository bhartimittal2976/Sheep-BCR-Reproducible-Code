import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from primer_common import revcomp,matches
from integrate_vdj import iso,category,reconstructed
class CoreTests(unittest.TestCase):
 def test_primers(self):
  self.assertEqual(revcomp('ATGR'),'YCAT')
  self.assertTrue(any(x[0]=='+' for x in matches('AACGT','ACGN',0)))
 def test_isotypes(self):
  self.assertEqual(iso('IGHG1'),'IgG');self.assertEqual(iso('IGHD'),'IgD');self.assertEqual(iso(''),'Other')
 def test_categories(self):
  self.assertEqual(category([{'chain':'IGH','c_gene':'IGHG'},{'chain':'IGL'}]),'IgG/IgL')
  self.assertEqual(category([{'chain':'IGH','c_gene':'IGHG'},{'chain':'IGH','c_gene':'IGHM'},{'chain':'IGL'}]),'multiple chains')
 def test_translation(self):
  r={x:'A' for x in ['fwr1','cdr1','fwr2','cdr2','fwr3','cdr3','fwr4']}
  r.update({x+'_nt':'GCT' for x in r})
  self.assertEqual(reconstructed(r)[2],'PASS')
if __name__=='__main__':unittest.main()
