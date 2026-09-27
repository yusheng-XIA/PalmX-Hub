#!/usr/bin/env python3
"""Normalize oil palm SV VCFs from SyRI, SVIM-asm and CuteSV for SURVIVOR.

- Converts contig names to Africa_hap2 canonical names chr01B..chr16B.
- Converts SyRI record IDs/ALT into SVTYPE and SVLEN.
- Keeps BND/TRA-like records even without length; filters length-bearing DEL/INS/DUP/INV shorter than min size.
- For CuteSV, also rewrites BND ALT mate contig names.
"""
import argparse, gzip, re, sys
from pathlib import Path

DROP_SYRI = {"SYN", "NOTAL", "SNP"}
SIZED = {"DEL", "INS", "DUP", "INV"}
TYPE_MAP = {"TRA":"TRA", "TRANS":"TRA", "BND":"BND", "HDR":"HDR", "CPG":"CPG", "CPL":"CPL", "TDM":"DUP", "DUP:TANDEM":"DUP", "DUP:INT":"DUP"}
END_RE = re.compile(r'(?:^|;)END=([0-9]+)(?:;|$)')
ID_PREFIX_RE = re.compile(r'([A-Za-z]+)')
INFO_ID_RE = re.compile(r'##INFO=<ID=([^,>]+)')
CONTIG_ID_RE = re.compile(r'##contig=<ID=([^,>]+)(.*)$')

def build_contig_map(path):
    m={}
    lengths={}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            a=line.rstrip('\n').split('\t')
            m[a[0]]=a[1]
            if len(a) >= 3 and a[2].isdigit():
                lengths[a[1]] = a[2]
    # Built-in aliases observed in this workspace.
    for i in range(1,17):
        m.setdefault(f'SL_af_Chr{i}', f'chr{i:02d}B')
        m.setdefault(f'chr{i:02d}', f'chr{i:02d}B')
        m.setdefault(f'chr{i}', f'chr{i:02d}B')
    return m, lengths

def op(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)

def canon(name, cmap):
    return cmap.get(name, name)

def rewrite_alt(alt, cmap):
    # Rewrite breakend mate names inside ALT, e.g. ]SL_af_Chr13:117[ -> ]chr13B:117[
    def repl(m):
        return m.group(1)+canon(m.group(2), cmap)+m.group(3)
    return re.sub(r'([\[\]])([^:\[\]]+)(:[0-9]+[\[\]])', repl, alt)

def get_info(fields, key):
    for f in fields:
        if f.startswith(key+'='):
            return f.split('=',1)[1]
    return None

def set_info(fields, key, val):
    out=[]; seen=False
    for f in fields:
        if f.startswith(key+'='):
            out.append(f'{key}={val}'); seen=True
        else:
            out.append(f)
    if not seen: out.append(f'{key}={val}')
    return out

def infer_type(rid, ref, alt):
    m=ID_PREFIX_RE.match(rid or '')
    if m:
        t=m.group(1).upper()
        return TYPE_MAP.get(t,t)
    a0=alt.split(',')[0]
    if a0.startswith('<') and a0.endswith('>'):
        t=a0[1:-1].upper()
        return TYPE_MAP.get(t,t)
    if '[' in a0 or ']' in a0: return 'BND'
    if len(ref)>len(a0): return 'DEL'
    if len(a0)>len(ref): return 'INS'
    return 'UNKNOWN'

def norm_type(t):
    if not t: return None
    t=t.upper()
    return TYPE_MAP.get(t,t.split(':',1)[0] if t.startswith('DUP:') else t)

def infer_svlen(pos, ref, alt, info):
    a0=alt.split(',')[0]
    if '[' in a0 or ']' in a0:
        return 0
    if a0.startswith('<'):
        m=END_RE.search(info)
        return abs(int(m.group(1))-int(pos)+1) if m else 0
    return abs(len(a0)-len(ref))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--caller', required=True, choices=['syri','svimasm','cutesv'])
    ap.add_argument('--in', dest='inp', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--sample', required=True)
    ap.add_argument('--contig-map', required=True)
    ap.add_argument('--min-svlen', type=int, default=50)
    ap.add_argument('--stats')
    args=ap.parse_args()
    cmap, contig_lengths=build_contig_map(args.contig_map)
    have_svtype=have_svlen=False
    seen_contigs=set(); n_in=n_out=n_drop_type=n_drop_size=0; by_type={}
    with op(args.inp) as fin, open(args.out,'w') as fout:
        for line in fin:
            if line.startswith('##'):
                if line.startswith('##INFO=<ID=SVTYPE,'): have_svtype=True
                if line.startswith('##INFO=<ID=SVLEN,'): have_svlen=True
                if line.startswith('##contig=<ID='):
                    mm=CONTIG_ID_RE.match(line.rstrip('\n'))
                    if mm:
                        cid=canon(mm.group(1), cmap)
                        if cid in seen_contigs: continue
                        seen_contigs.add(cid)
                        fout.write(f'##contig=<ID={cid}{mm.group(2)}\n')
                        continue
                fout.write(line)
                continue
            if line.startswith('#CHROM'):
                for cid in sorted(contig_lengths):
                    if cid not in seen_contigs:
                        fout.write(f'##contig=<ID={cid},length={contig_lengths[cid]}>\n')
                        seen_contigs.add(cid)
                if not have_svtype:
                    fout.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Normalized SV type">\n')
                if not have_svlen:
                    fout.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Normalized SV length">\n')
                cols = line.rstrip('\n').split('\t')
                if len(cols) > 9:
                    cols[9:] = [args.sample]
                    fout.write('\t'.join(cols) + '\n')
                else:
                    fout.write(line)
                continue
            if line.startswith('#'):
                fout.write(line); continue
            n_in += 1
            cols=line.rstrip('\n').split('\t')
            if len(cols)<8: continue
            chrom,pos,rid,ref,alt=cols[0],cols[1],cols[2],cols[3],cols[4]
            if args.caller=='cutesv':
                alt=rewrite_alt(alt, cmap); cols[4]=alt
            info=cols[7]
            fields=[] if info in ('','.') else info.split(';')
            if args.caller=='syri':
                svtype=infer_type(rid, ref, alt)
                if svtype in DROP_SYRI or svtype.endswith('AL'):
                    n_drop_type += 1; continue
            else:
                svtype=norm_type(get_info(fields,'SVTYPE')) or infer_type(rid, ref, alt)
            fields=set_info(fields,'SVTYPE',svtype)
            raw=get_info(fields,'SVLEN')
            if raw is None:
                svlen=infer_svlen(pos,ref,alt,info)
                fields=set_info(fields,'SVLEN',svlen)
            else:
                try: svlen=int(str(raw).split(',')[0])
                except Exception: svlen=0
            if svtype in SIZED and abs(svlen)<args.min_svlen:
                n_drop_size += 1; continue
            cols[0]=canon(chrom,cmap)
            cols[7]=';'.join(fields) if fields else '.'
            fout.write('\t'.join(cols)+'\n')
            n_out += 1
            by_type[svtype]=by_type.get(svtype,0)+1
    summary=f"{args.sample}\t{args.caller}\t{n_in}\t{n_out}\t{n_drop_type}\t{n_drop_size}\t" + ','.join(f'{k}:{v}' for k,v in sorted(by_type.items()))
    sys.stderr.write('[normalize]\t'+summary+'\n')
    if args.stats:
        with open(args.stats,'a') as s: s.write(summary+'\n')
if __name__=='__main__': main()
