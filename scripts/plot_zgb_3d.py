#!/usr/bin/env python3
"""Render actual initial/final 3D XYZ particles; no height-map extrusion."""
import argparse
import json
from pathlib import Path
import re

def read_xyz(path):
    lines=path.read_text().splitlines()
    count=int(lines[0]);cell=[float(v) for v in re.search(r'Lattice="([^"]+)"',lines[1]).group(1).split()]
    particles=[]
    for row in lines[2:2+count]:
        t=row.split();particles.append({'position':tuple(float(v) for v in t[1:4]),'state':t[6],'role':t[7]})
    if len(particles)!=count:raise ValueError('truncated XYZ')
    return particles,cell

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,required=True)
    p.add_argument('--output',type=Path)
    p.add_argument('--substrate-layers',type=int,default=2)
    args=p.parse_args()
    if args.substrate_layers<1:p.error('--substrate-layers must be positive')
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
    except ImportError:p.error('Plotting requires matplotlib: python3 -m pip install matplotlib')
    prefix='surface' if (args.results/'surface_summary.json').exists() else 'zgb3d'
    summary=json.loads((args.results/(prefix+'_summary.json')).read_text())
    initial,cell=read_xyz(args.results/(prefix+'_initial.xyz'));final,_=read_xyz(args.results/(prefix+'_final.xyz'))
    a=cell[0]/summary['repetitions'][0]
    top=max(q['position'][2] for q in initial if q['state']=='substrate' or q['role']=='substrate')
    cutoff=top-(args.substrate_layers-1)*a-1e-8
    sets=[initial,final]
    zmax=max(q['position'][2] for q in final)
    fig=plt.figure(figsize=(12,6),layout='constrained')
    palette={'substrate':('#8b9299','o','Substrate'), 'AB2_solid':('#2575ae','o','AB2 solid'),
             'AB_adsorbate':('#df9b29','^','AB adsorbate'), 'B_adsorbate':('#b93840','^','B adsorbate')}
    if prefix=='surface':
        palette={}
        for q in final+initial:
            state=q['state'];key=(state,q['role'])
            if key in palette:continue
            color='#8b9299' if q['role']=='substrate' else plt.get_cmap('tab20')(len(palette)%20)
            palette[key]=(color,'^' if q['role']=='adsorbate' else 'o',state+' ('+q['role']+')')
    for index,particles in enumerate(sets):
        ax=fig.add_subplot(1,2,index+1,projection='3d')
        visible=[q for q in particles if q['position'][2]>=cutoff]
        points={tuple(round(v,7) for v in q['position']):q for q in visible}
        bonds=[]
        for q in visible:
            x,y,z=q['position'];below=(round(x,7),round(y,7),round(z-a,7))
            if below in points and points[below]['role'] in ('solid','film','substrate'):bonds.append([q['position'],points[below]['position']])
        ax.add_collection3d(Line3DCollection(bonds,colors='#7c8790',linewidths=0.6,alpha=0.6))
        for state,(color,marker,label) in palette.items():
            selected=[q for q in visible if ((q['state'],q['role'])==state if isinstance(state,tuple) else q['state']==state)]
            if selected:ax.scatter(*zip(*(q['position'] for q in selected)),s=13,marker=marker,color=color,label=label,depthshade=True)
        ax.set(xlim=(-a,cell[0]),ylim=(-a,cell[4]),zlim=(cutoff-a,zmax+a),xlabel='x (Å)',ylabel='y (Å)',zlabel='z (Å)',
               title='Initial substrate' if index==0 else f"After {summary['events']:,} events: {summary['solid_units']:,} solid units")
        ax.set_box_aspect((cell[0],cell[4],max(zmax-cutoff+2*a,5*a)))
        ax.view_init(elev=25,azim=-55)
        if index==1:ax.legend(loc='upper left',fontsize=8)
    fig.suptitle('Explicit 3D lattice occupancy and chemical states')
    fig.text(0.5,0.01,f'Lower substrate hidden; top {args.substrate_layers} substrate layers shown. Lines: vertical solid support. Markers: effective chemical units, not resolved molecules.',ha='center',fontsize=9)
    output=args.output or args.results/'surface-3d.png'
    fig.savefig(output,dpi=180);plt.close(fig);print(output)

if __name__=='__main__':main()
