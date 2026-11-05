#!/usr/bin/env python3
"""Plot the CSV reproduction outputs. Requires matplotlib (optional for the solver)."""
import argparse
import csv
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,required=True)
    p.add_argument('--output',type=Path)
    p.add_argument('--reference-surface',type=Path)
    args=p.parse_args()
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        p.error('Plotting requires matplotlib: python3 -m pip install matplotlib')
    rows=list(csv.DictReader((args.results/'comparison.csv').open()))
    output=args.output or args.results/'figures';output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    groups={variant:sorted([r for r in rows if r['variant']==variant],key=lambda r:float(r['y_ab'])) for variant in ['growth','catalytic']}
    growth=groups['growth']
    fig,ax=plt.subplots(figsize=(7,4.5),layout='constrained')
    for variant,items in groups.items():
        if not items:continue
        x=[float(r['y_ab']) for r in items]
        for key,color,label in [('coverage_ab','black','AB*'),('coverage_b','#b82d33','B*')]:
            ax.errorbar(x,[float(r[key+'_mean']) for r in items],yerr=[float(r[key+'_sem']) for r in items],
                        fmt='s' if variant=='growth' else '--',color=color,capsize=3,label=f'{label}, {variant}')
    for x in [0.389,0.525]:ax.axvline(x,color='0.65',lw=0.8,ls=':')
    ax.set(xlabel=r'$y_{AB}$',ylabel='Surface coverage',ylim=(-0.025,1.025),title='ZGB reference vs growth: trial means ± SEM')
    ax.legend(fontsize=9);fig.savefig(output/'coverage.png',dpi=180);plt.close(fig)
    if growth:
        x=[float(r['y_ab']) for r in growth]
        fig,axes=plt.subplots(1,2,figsize=(9,3.7),layout='constrained')
        for ax,key,label in zip(axes,['mean_growth_rate','time_average_height'],['Mean growth rate (ML/s)','Time-averaged mean height (ML)']):
            ax.errorbar(x,[float(r[key+'_mean']) for r in growth],yerr=[float(r[key+'_sem']) for r in growth],fmt='o-',capsize=3,color='black')
            ax.set(xlabel=r'$y_{AB}$',ylabel=label)
        fig.suptitle('Growth observables averaged over time until poisoning')
        fig.savefig(output/'growth.png',dpi=180);plt.close(fig)
        comparable=[r for r in growth if r['paper_table2_time_s']]
        if comparable:
            fig,ax=plt.subplots(figsize=(7,4.5),layout='constrained')
            x=[float(r['y_ab']) for r in comparable]
            ax.plot(x,[float(r['paper_table2_time_s']) for r in comparable],'s--',label='Paper Table 2',color='0.4')
            ax.errorbar(x,[float(r['time_mean']) for r in comparable],yerr=[float(r['time_sem']) for r in comparable],fmt='o-',capsize=3,label='New runs ± SEM',color='#24618b')
            ax.set(xlabel=r'$y_{AB}$',ylabel='Absorption time (s)',title='Poisoning-time comparison; paper uncertainty is unavailable')
            ax.legend();fig.savefig(output/'poisoning-time.png',dpi=180);plt.close(fig)
    # Figure-3-style first-trial snapshots, with exact coordinates from output.
    candidate=args.results/'growth-y0.2-trial1'/'zgb_final.csv'
    if candidate.exists():
        surfaces=[(candidate,'Growth: surface species','species'),(candidate,'Growth: column heights','height')]
        if args.reference_surface:surfaces.insert(0,(args.reference_surface,'Catalytic reference: species','species'))
        fig,axes=plt.subplots(1,len(surfaces),figsize=(4.3*len(surfaces),4),layout='constrained',squeeze=False)
        from matplotlib.colors import ListedColormap,BoundaryNorm
        for ax,(path,title,kind) in zip(axes[0],surfaces):
            cells=list(csv.DictReader(path.open()));nx=1+max(int(c['x']) for c in cells);ny=1+max(int(c['y']) for c in cells)
            matrix=[[0]*nx for _ in range(ny)]
            for c in cells:matrix[int(c['y'])][int(c['x'])]=int(c['height']) if kind=='height' else {'vacant':0,'B':1,'AB':2}[c['adsorbate']]
            image=ax.imshow(matrix,origin='lower',interpolation='nearest',**({'cmap':'viridis'} if kind=='height' else {'cmap':ListedColormap(['#eeeeee','#929292','#202020']),'norm':BoundaryNorm([-0.5,0.5,1.5,2.5],3)}))
            ax.set(xlabel='x site',ylabel='y site',title=title)
            bar=fig.colorbar(image,ax=ax,shrink=0.75)
            if kind=='species':bar.set_ticks([0,1,2],labels=['vacant','B*','AB*'])
            else:bar.set_label('Height (ML)')
        fig.savefig(output/'surface-snapshots.png',dpi=180);plt.close(fig)
    print('Plots saved to '+str(output))

if __name__=='__main__':main()
