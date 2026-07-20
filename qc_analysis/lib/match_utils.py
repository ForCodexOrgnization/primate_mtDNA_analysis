"""Small standard-library helpers shared by downstream QC annotation modules."""
from __future__ import annotations
import csv, gzip
from pathlib import Path

def scalar(v):
 v=v.strip()
 if not v or v.lower() in {'null','none','~'}: return None
 if (v[:1],v[-1:]) in [("'","'"),('"','"')]: return v[1:-1]
 if v.lower() in {'true','false'}: return v.lower()=='true'
 try:return int(v)
 except ValueError:
  try:return float(v)
  except ValueError:return v

def yaml(path):
 root={}; stack=[(-1,root)]
 for raw in Path(path).open():
  line=raw.split('#',1)[0].rstrip()
  if not line.strip(): continue
  n=len(line)-len(line.lstrip()); key,val=line.strip().split(':',1); val=val.strip()
  while n<=stack[-1][0]: stack.pop()
  if val: stack[-1][1][key]=scalar(val)
  else:
   d={}; stack[-1][1][key]=d; stack.append((n,d))
 return root

def info_parse(s):
 return {} if s in {'','.','None'} else {x.split('=',1)[0]:x.split('=',1)[1] if '=' in x else True for x in s.split(';')}
def info_format(d): return ';'.join(k if v is True else f'{k}={v}' for k,v in d.items()) or '.'
def source(info):
 """Return original source chrom, position, ref, alt across supported liftover INFO conventions."""
 def g(a,b): return info.get(a,info.get(b,''))
 try: pos=int(g('SRC_POS','MTLIFT_ORIG_POS'))
 except (ValueError,TypeError): pos=None
 return g('SRC_CHROM','MTLIFT_ORIG_CHROM'),pos,g('SRC_REF','MTLIFT_ORIG_REF'),g('SRC_ALT','MTLIFT_ORIG_ALT')
def human_pos(fields,info):
 try:return int(info.get('MTLIFT_HUMAN_POS',fields[1]))
 except (ValueError,TypeError): return None
def open_text(p,mode='rt'): return gzip.open(p,mode) if str(p).endswith('.gz') else open(p,mode)
def rows(path):
 with open_text(path) as f: return list(csv.DictReader(f,delimiter='\t'))
def write_summary(path, row):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(row));w.writeheader();w.writerow(row)
def sample_names(cfg):
 p=cfg.get('coordinate_liftover',{}).get('paths',{}).get('sample_ref_file')
 if not p or not Path(p).exists(): return []
 with open(p) as f:return [r.get('sample','') for r in csv.DictReader(f,delimiter='\t') if r.get('sample')]
def inject_headers(header, fields, prefix):
 seen='\n'.join(header); out=[]
 for name,desc in fields:
  if f'ID={name},' not in seen: out.append(f'##INFO=<ID={name},Number=1,Type=String,Description="{desc}">\n')
 for i,x in enumerate(header):
  if x.startswith('#CHROM'): return header[:i]+out+header[i:]
 return header+out

def normalize_rna_base(base):
 """Return a single RNA base (DNA thymine is represented as uracil)."""
 b=str(base or '').strip().upper().replace('T','U')
 return b if b in {'A','C','G','U'} else None

def pair_type(base1,base2):
 a,b=normalize_rna_base(base1),normalize_rna_base(base2)
 if not a or not b:return 'NA'
 if (a,b) in {('A','U'),('U','A'),('G','C'),('C','G')}:return 'WC'
 if (a,b) in {('G','U'),('U','G')}:return 'GU_wobble'
 return 'non_WC'

def pair_state(kind):
 return 'NA' if kind in {'',None,'.','NA'} else ('paired' if kind in {'WC','GU_wobble','non_WC'} else str(kind))

def pair_effect(ref_pair_type,alt_pair_type):
 if ref_pair_type in {'',None,'.','NA'} or alt_pair_type in {'',None,'.','NA'}:return 'NA'
 return 'unchanged' if ref_pair_type==alt_pair_type else f'{ref_pair_type}_to_{alt_pair_type}'

def compare_values(a,b):
 if a in {'',None,'.','NA'} or b in {'',None,'.','NA'}:return '.'
 return 'yes' if str(a)==str(b) else 'no'

def load_coordinate_map(path):
 """Index an existing liftover map by original source position."""
 result={}
 if not path or not Path(path).exists():return result
 for row in rows(path):
  try: result[int(row.get('species_pos_original',row.get('source_pos','')))] = row.get('human_pos_canonical',row.get('human_pos',''))
  except (TypeError,ValueError):pass
 return result

def lift_source_pos_to_human(pos,coordinate_map):
 try:return str(coordinate_map.get(int(pos),'.') or '.')
 except (TypeError,ValueError):return '.'
