#!/usr/bin/env python3
"""Generate Supplementary Figure S1 from the authoritative 207-system CSV."""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, leaves_list, to_tree

ROOT = Path(__file__).resolve().parent.parent
CSV = ROOT / "source_data" / "Supplementary_Data_S1_system_inventory.csv"
OUT = ROOT / "figures"
OUT.mkdir(parents=True, exist_ok=True)

F = ["Gi", "Gs", "Gq", "G12-13"]
FL = {"Gi": "Gi/o", "Gs": "Gs", "Gq": "Gq/11", "G12-13": "G12/13"}
C = {"Gi": "#2F8F6B", "Gs": "#2C6FB3", "Gq": "#C0741A", "G12-13": "#8A4AA0"}
INK = "#20242B"
MUTED = "#68707A"
plt.rcParams.update({
 "font.family": "sans-serif",
 "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
 "font.size": 8.2,
 "pdf.fonttype": 42,
 "ps.fonttype": 42,
})

s = pd.read_csv(CSV)
assert len(s) == 207
assert s.gpcr_class.value_counts().to_dict() == {"A": 181, "B": 26}
assert s.g_protein_family.value_counts().to_dict() == {"Gi": 95, "Gs": 65, "Gq": 41, "G12-13": 6}
assert s.receptor_name.nunique()==174 and s.receptor_uniprot.nunique()==173
assert s.loc[s.receptor_uniprot.isna(),'system_id'].tolist()==['Gs_8HTI']
g=s.groupby('receptor_name');names=sorted(g.groups)
profiles={n:set(g.get_group(n).g_protein_family) for n in names}
X=np.array([[f in profiles[n] for f in F] for n in names],float)
Z=linkage(X,method='average',metric='jaccard');order=leaves_list(Z);root=to_tree(Z)
angles=np.zeros(len(names));angles[order]=np.pi/2-2*np.pi*np.arange(len(names))/len(names)
info={}
def walk(node):
 if node.is_leaf():info[id(node)]=(angles[node.id],angles[node.id],1.0);return info[id(node)]
 l,r=walk(node.left),walk(node.right);out=(l[0],r[1],.19+.77*(1-node.dist/(root.dist or 1)));info[id(node)]=out;return out
walk(root)
fig,ax=plt.subplots(figsize=(13.5,13.5));ax.set_aspect('equal');ax.axis('off');ax.set(xlim=(-1.33,1.33),ylim=(-1.33,1.33))
def arc(rad,a,b):
 t=np.linspace(a,b,50);ax.plot(rad*np.cos(t),rad*np.sin(t),c='#A9AFB7',lw=.55,zorder=1)
def draw(node):
 a,b,r=info[id(node)]
 if node.is_leaf():return
 for child in (node.left,node.right):
  ca,cb,cr=info[id(child)];mid=(ca+cb)/2
  ax.plot([r*np.cos(mid),cr*np.cos(mid)],[r*np.sin(mid),cr*np.sin(mid)],c='#A9AFB7',lw=.55,zorder=1);draw(child)
 arc(r,a,b)
draw(root)
unmapped_name='Consensus Olfactory Receptor OR52c'
for i,n in enumerate(names):
 a=angles[i];fs=profiles[n];x,y=np.cos(a),np.sin(a)
 col='#4A4F55' if len(fs)>1 else C[next(iter(fs))]
 ax.plot([.985*x,1.035*x],[.985*y,1.035*y],c=col,lw=2.0,solid_capstyle='butt',zorder=3)
 sub=g.get_group(n);fam=next(f for f in F if f in fs);pdb=str(sub[sub.g_protein_family.eq(fam)].sort_values('pdb_id').pdb_id.iloc[0])
 ha='left' if x>=0 else 'right';rot=np.degrees(a) if x>=0 else np.degrees(a)+180
 suffix=' *' if n==unmapped_name else ''
 ax.text(1.058*x,1.058*y,pdb+suffix,fontsize=8.2,rotation=rot,rotation_mode='anchor',ha=ha,va='center',color=col,fontweight='bold' if len(fs)>1 or n==unmapped_name else 'normal')
ax.add_patch(plt.Circle((0,0),.32,fc='white',ec='#D9DDE2',lw=.8,zorder=2))
ax.text(0,.055,'174',ha='center',va='center',fontsize=24,fontweight='bold',color=INK,zorder=4)
ax.text(0,-.025,'distinct receptor names',ha='center',va='center',fontsize=10.2,color=INK,zorder=4)
ax.text(0,-.09,'presence/absence profiles across\nfour G-protein families',ha='center',va='center',fontsize=8.4,color=MUTED,zorder=4)
handles=[Line2D([0],[0],color=C[f],lw=3,label=FL[f]) for f in F]+[Line2D([0],[0],color='#4A4F55',lw=3,label='>1 family')]
ax.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,-.01),ncol=5,frameon=False,fontsize=8.5,handlelength=1.7,columnspacing=1.5)
ax.text(-1.28,1.27,'S1',fontweight='bold',fontsize=15.5,ha='left',va='top')
ax.text(0,-1.22,'Leaf label = representative PDB identifier. * Gs_8HTI has no mapped UniProt accession.',ha='center',fontsize=8.2,color=INK)
ax.text(0,-1.27,'Branch proximity reflects release-coverage profiles only; it does not represent sequence, structural, or mechanistic similarity.',ha='center',fontsize=8.2,color=MUTED)
for ext in ('pdf','png'):fig.savefig(OUT/f'supplementary_figure_s1_receptor_profiles.{ext}',dpi=600 if ext=='png' else None,bbox_inches='tight',facecolor='white')
plt.close(fig)
print(f'S1 assertions passed: 207 systems; 174 receptor names; 173 mapped UniProt accessions; Gs_8HTI unmapped')
print(f'Wrote PDF/PNG pair to {OUT}')
