#!/usr/bin/env python3
import argparse, csv, json
from bisect import bisect_left
from pathlib import Path

def read_tsv(path):
    with open(path) as fh: return list(csv.DictReader(fh, delimiter='\t'))

def write_tsv(path, rows, fields):
    path=Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with open(path,'w',newline='') as fh:
        w=csv.DictWriter(fh, fieldnames=fields, delimiter='\t'); w.writeheader(); w.writerows(rows)

def load_local(path, phase_col='local_phase_state', conf_col='confidence'):
    out=[]
    for r in read_tsv(path):
        try: out.append({'pos':int(r['pos']),'block_id':r.get('block_id',''),'phase_state':int(r.get(phase_col,0) or 0),'confidence':float(r.get(conf_col,0.0) or 0.0)})
        except: pass
    return sorted(out,key=lambda x:x['pos'])

def load_scaffold(path):
    out=[]
    for r in read_tsv(path):
        try: out.append({'pos':int(r['pos']),'phase_state':int(r.get('scaffold_phase_state',0) or 0),'confidence':float(r.get('scaffold_confidence',0.0) or 0.0)})
        except: pass
    return sorted(out,key=lambda x:x['pos'])

def nearest(scaffold,pos,max_gap):
    positions=[r['pos'] for r in scaffold if r['phase_state']!=0]
    idx=bisect_left(positions,pos)
    left=positions[idx-1] if idx>0 else None
    right=positions[idx] if idx<len(positions) else None
    if left is not None and pos-left>max_gap: left=None
    if right is not None and right-pos>max_gap: right=None
    return left,right

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--local-calls-tsv', required=True)
    ap.add_argument('--scaffold-calls-tsv', required=True)
    ap.add_argument('--min-local-confidence', type=float, default=0.70)
    ap.add_argument('--max-anchor-gap', type=int, default=250000)
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--out-summary-json', required=True)
    args=ap.parse_args()
    local=load_local(args.local_calls_tsv); scaffold=load_scaffold(args.scaffold_calls_tsv); sby={r['pos']:r for r in scaffold}
    rows=[]; q=0
    for lr in local:
        pos=lr['pos']; already=bool(sby.get(pos,{}).get('phase_state',0))
        unresolved=(lr['phase_state']==0) or (lr['confidence']<args.min_local_confidence) or (not already)
        if not unresolved: continue
        left,right=nearest(scaffold,pos,args.max_anchor_gap)
        left_state=sby.get(left,{}).get('phase_state',0) if left else 0
        right_state=sby.get(right,{}).get('phase_state',0) if right else 0
        bridge_hint=left_state if left_state!=0 and left_state==right_state else 0
        qid=f'q_{pos}'; q+=1
        common={'query_id':qid,'pos':pos,'local_block_id':lr['block_id'],'local_phase_state':lr['phase_state'],'local_confidence':f"{lr['confidence']:.6f}",'already_scaffolded':int(already),'left_anchor_pos':left or '','right_anchor_pos':right or '','left_anchor_state':left_state,'right_anchor_state':right_state,'bridge_hint_state':bridge_hint,'distance_left_anchor':(pos-left) if left else '','distance_right_anchor':(right-pos) if right else ''}
        for cid,state,ctype in [('A',1,'hap1'),('B',-1,'hap2'),('ABSTAIN',0,'abstain')]:
            r=dict(common); r.update({'candidate_id':cid,'candidate_type':ctype,'candidate_phase_state':state}); rows.append(r)
    fields=list(rows[0].keys()) if rows else []
    write_tsv(args.out_tsv, rows, fields)
    summary={'queries':q,'candidate_rows':len(rows)}
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True); json.dump(summary, open(args.out_summary_json,'w'), indent=2); print(json.dumps(summary, indent=2))

if __name__=='__main__': main()
