#!/usr/bin/env python3
"""Plot actual exported C/H atoms and explicit bonds from a diamond CVD run."""
import argparse
import json
import math
from pathlib import Path
import re

def read_xyz(path):
    lines=path.read_text().splitlines();n=int(lines[0]);cell=list(map(float,re.search(r'Lattice="([^"]+)"',lines[1]).group(1).split()))
    atoms={}
    for line in lines[2:2+n]:
        t=line.split();atoms[int(t[4])]={'element':t[0],'r':tuple(map(float,t[1:4])),'fixed':t[5]=='T','state':t[6],'parent':int(t[7]),'trapped':t[8]=='T'}
    if len(atoms)!=n:raise ValueError('incomplete or duplicate XYZ atoms')
    return atoms,[cell[0],cell[4],cell[8]]

def displacement(a,b,box):
    d=[b[k]-a[k] for k in range(3)]
    for k in [0,1]:d[k]-=round(d[k]/box[k])*box[k]
    return d

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True);p.add_argument('--output',type=Path);a=p.parse_args()
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
    except ImportError:p.error('Matplotlib is required')
    initial,box=read_xyz(a.results/'diamond_initial.xyz');final,_=read_xyz(a.results/'diamond_final.xyz')
    summary=json.loads((a.results/'diamond_summary.json').read_text())
    bond_sets=[[],[]];carbon=[i for i,v in initial.items() if v['element']=='C']
    for index,i in enumerate(carbon):
        for j in carbon[:index]:
            d=displacement(initial[i]['r'],initial[j]['r'],box)
            if math.sqrt(sum(v*v for v in d))<1.6:bond_sets[0].append((i,j,'carbon_sigma'))
    for i,v in initial.items():
        if v['element']=='H':bond_sets[0].append((v['parent'],i,'carbon_hydrogen'))
    for line in (a.results/'diamond_bonds.tsv').read_text().splitlines()[1:]:
        t=line.split();bond_sets[1].append((int(t[0]),int(t[1]),t[2]))
    fig=plt.figure(figsize=(12,7));fig.subplots_adjust(left=.02,right=.98,bottom=.14,top=.88,wspace=.08);zmax=max(v['r'][2] for v in final.values())
    for k,atoms in enumerate([initial,final]):
        ax=fig.add_subplot(1,2,k+1,projection='3d')
        for kind,color,width in [('carbon_sigma','#7f8e9a',1.2),('carbon_hydrogen','#d4a0a0',.8),('dimer','#9b429f',2.5)]:
            segments=[]
            for i,j,tag in bond_sets[k]:
                if tag!=kind:continue
                start=atoms[i]['r'];d=displacement(start,atoms[j]['r'],box);segments.append([start,tuple(start[c]+d[c] for c in range(3))])
            if segments:ax.add_collection3d(Line3DCollection(segments,colors=color,linewidths=width))
        categories=[('Initial diamond','#8997a2',lambda i,v:v['state']=='diamond' and i in initial),('Grown diamond','#247ab9',lambda i,v:v['state']=='diamond' and i not in initial),('Hydrocarbon C','#e99b32',lambda i,v:v['state']=='hydrocarbon'),('Hydrogen','#c74747',lambda i,v:v['element']=='H')]
        for label,color,select in categories:
            chosen=[v['r'] for i,v in atoms.items() if select(i,v)]
            if chosen:ax.scatter(*zip(*chosen),s=24 if label=='Hydrogen' else 45,color=color,label=label,depthshade=True)
        ax.set(xlim=(-.3,box[0]+.3),ylim=(-.3,box[1]+.3),zlim=(-.3,zmax+.7),xlabel='x (Å)',ylabel='y (Å)',zlabel='z (Å)',title='Initial H-passivated substrate' if k==0 else f"After {summary['time']:.5f} s: +{summary['deposited_diamond_carbon_atoms']} diamond C")
        ax.set_box_aspect((box[0],box[1],zmax+1));ax.view_init(elev=22,azim=-55)
        if k:ax.legend(loc='upper left',fontsize=8)
    fig.suptitle('Diamond (111) CVD: explicit C/H atoms and bonds')
    fig.text(.5,.035,f"{summary['events']:,} events · {summary['temperature_K']:.0f} K · stop: {summary['stop_reason']} · periodic x/y\nInitial mechanism run; published growth curves not validated",ha='center',fontsize=9)
    output=a.output or a.results/'diamond-growth.png';fig.savefig(output,dpi=180);plt.close(fig);print(output)
if __name__=='__main__':main()
