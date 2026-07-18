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
