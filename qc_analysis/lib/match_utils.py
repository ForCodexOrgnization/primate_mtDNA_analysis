"""Small standard-library helpers shared by downstream QC annotation modules."""
from __future__ import annotations
import csv, gzip, os, tempfile
from pathlib import Path

RESOLVED_DNA_BASES = frozenset("ACGT")
IUPAC_DNA_BASES = frozenset("ACGTRYSWKMBDHVN")
RESOLVED_RNA_BASES = frozenset("ACGU")
IUPAC_RNA_BASES = frozenset("ACGURYSWKMBDHVN")
DNA_IUPAC_COMPLEMENT = {
 'A':'T','T':'A','C':'G','G':'C','R':'Y','Y':'R','S':'S','W':'W',
 'K':'M','M':'K','B':'V','V':'B','D':'H','H':'D','N':'N',
}
RNA_IUPAC_STATES = {
 'A':frozenset('A'),'C':frozenset('C'),'G':frozenset('G'),'U':frozenset('U'),
 'R':frozenset('AG'),'Y':frozenset('CU'),'S':frozenset('CG'),'W':frozenset('AU'),
 'K':frozenset('GU'),'M':frozenset('AC'),'B':frozenset('CGU'),
 'D':frozenset('AGU'),'H':frozenset('ACU'),'V':frozenset('ACG'),
 'N':frozenset('ACGU'),
}

def _upper_symbol(base): return str(base or '').strip().upper()
def is_resolved_dna_base(base): return _upper_symbol(base) in RESOLVED_DNA_BASES
def is_iupac_dna_base(base): return _upper_symbol(base) in IUPAC_DNA_BASES
def is_resolved_rna_base(base): return _upper_symbol(base).replace('T','U') in RESOLVED_RNA_BASES
def is_iupac_rna_base(base): return _upper_symbol(base).replace('T','U') in IUPAC_RNA_BASES

def normalize_rna_symbol(base):
 """Return an RNA IUPAC symbol, or ``None`` for missing/invalid input."""
 b=_upper_symbol(base).replace('T','U')
 return b if b in IUPAC_RNA_BASES else None

def orient_dna_base_to_rna(base, strand='+'):
 """Convert one DNA IUPAC symbol to mature-RNA orientation."""
 b=_upper_symbol(base)
 if b not in IUPAC_DNA_BASES:return None
 if strand=='-':b=DNA_IUPAC_COMPLEMENT[b]
 return b.replace('T','U')

def rna_symbols_compatible(left,right):
 """Whether two valid RNA IUPAC symbols have at least one shared state."""
 a,b=normalize_rna_symbol(left),normalize_rna_symbol(right)
 return bool(a and b and RNA_IUPAC_STATES[a] & RNA_IUPAC_STATES[b])

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
 fd,tmp=tempfile.mkstemp(prefix=f'.{path.name}.',suffix='.tmp',dir=path.parent,text=True)
 with os.fdopen(fd,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(row),delimiter='\t');w.writeheader();w.writerow(row)
  f.flush(); os.fsync(f.fileno())
 os.replace(tmp,path)
def sample_names(cfg):
 p=cfg.get('coordinate_liftover',{}).get('paths',{}).get('sample_ref_file')
 if not p or not Path(p).exists(): return []
 with open(p) as f:return [r.get('sample','') for r in csv.DictReader(f,delimiter='\t') if r.get('sample')]
def inject_headers(header, fields, prefix):
 seen='\n'.join(header); out=[]
 for spec in fields:
  if len(spec)==2: name,desc=spec; number,type_='1','String'
  elif len(spec)==4: name,number,type_,desc=spec
  else: raise ValueError(f'Invalid INFO field specification: {spec!r}')
  if f'ID={name},' not in seen: out.append(f'##INFO=<ID={name},Number=1,Type=String,Description="{desc}">\n')
  if out and out[-1].startswith(f'##INFO=<ID={name},'):
   out[-1]=f'##INFO=<ID={name},Number={number},Type={type_},Description="{desc}">\n'
 for i,x in enumerate(header):
  if x.startswith('#CHROM'): return header[:i]+out+header[i:]
 return header+out

def normalize_rna_base(base):
 """Return a fully resolved RNA base; ambiguity is deliberately unresolved."""
 b=normalize_rna_symbol(base)
 return b if b in RESOLVED_RNA_BASES else None

def pair_type(base1,base2):
 a,b=normalize_rna_symbol(base1),normalize_rna_symbol(base2)
 if not a or not b:return 'NA'
 if a not in RESOLVED_RNA_BASES or b not in RESOLVED_RNA_BASES:return 'ambiguous'
 if (a,b) in {('A','U'),('U','A'),('G','C'),('C','G')}:return 'WC'
 if (a,b) in {('G','U'),('U','G')}:return 'GU_wobble'
 return 'non_WC'

def pair_state(kind):
 return 'NA' if kind in {'',None,'.','NA','ambiguous'} else ('paired' if kind in {'WC','GU_wobble','non_WC'} else str(kind))

def pair_effect(ref_pair_type,alt_pair_type):
 if ref_pair_type in {'',None,'.','NA','ambiguous'} or alt_pair_type in {'',None,'.','NA','ambiguous'}:return 'NA'
 return 'unchanged' if ref_pair_type==alt_pair_type else f'{ref_pair_type}_to_{alt_pair_type}'

def rrna_pair_type(base1,base2):
 """Return explicit RNA base-pair labels for rRNA structure annotation."""
 a,b=normalize_rna_symbol(base1),normalize_rna_symbol(base2)
 if not a or not b:return '.'
 if a not in RESOLVED_RNA_BASES or b not in RESOLVED_RNA_BASES:return 'other'
 label=f'{a}-{b}'
 return label if label in {'A-U','U-A','G-C','C-G','G-U','U-G'} else 'other'

def rrna_pair_state(kind,struct_class='stem'):
 """Classify an explicit rRNA pair as canonical, wobble, noncanonical, or unknown."""
 k=str(kind or '').strip()
 c=str(struct_class or '').strip().lower()
 if c=='loop':return 'unpaired'
 if c=='unknown' or k in {'','.','NA','unknown'}:return 'unknown'
 if k in {'A-U','U-A','G-C','C-G'}:return 'canonical'
 if k in {'G-U','U-G'}:return 'wobble'
 return 'noncanonical'

def rrna_pair_effect(ref_pair_type,alt_pair_type):
 """Describe how an alternate allele changes the human-reference pair category."""
 ref_state=rrna_pair_state(ref_pair_type)
 alt_state=rrna_pair_state(alt_pair_type)
 if 'unknown' in {ref_state,alt_state}:return 'NA'
 return 'unchanged' if ref_state==alt_state else f'{ref_state}_to_{alt_state}'

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
