#!/usr/bin/env python3
"""Create the Battaile/Srolovitz/Butler (1997) diamond (111) reference input.

All 12 Table II rows are supplied as 22 explicit directional processes. Geometry
is fixed carbon sites with an explicit hydrogen/bond-state simulation backend.
Gas temperature and lattice constant are stated implementation assumptions.
"""
import argparse
import json
import math
from pathlib import Path

# Exact fractional basis in an orthorhombic (111) cell. In conventional cubic
# coordinates its vectors are a/2[1,-1,0], a/2[1,1,-2], and a[1,1,1].
BASIS=[(0,0,0),(.5,.5,0),(0,0,.25),(.5,.5,.25),
       (0,1/3,1/3),(.5,5/6,1/3),(0,1/3,7/12),(.5,5/6,7/12),
       (.5,1/6,2/3),(0,2/3,2/3),(.5,1/6,11/12),(0,2/3,11/12)]
# Table II: A, n, E [kcal/mol], delta H [kcal/mol], delta S [cal/(mol K)].
TABLE={1:(1.3e14,0,7.3,-9.9,5.3),2:(1e13,0,0,-96.9,-32.8),
 3:(3e13,0,0,-24.6,7.9),4:(5e12,0,0,-70.9,-42.0),5:(8e10,0,7.7,-19.0,-10.7),
 6:(2.8e7,2,7.7,-11.3,6.6),7:(1e13,0,0,-83.0,-34.1),
 8:(9e6,2,5.0,-8.9,8.7),9:(2e13,0,0,-47.7,-36.2),
 10:(5e12,0,0,0,0),11:(1e13,0,0,4.9,.4),12:(2e13,0,8.8,0,0)}

def make_case(nx=2,ny=2,cells=6,temperature=1200,gas_temperature=1200,acetylene=.04,events=200000,duration=2):
    a=3.567;lx=a/math.sqrt(2);ly=a*math.sqrt(6)/2;lz=a*math.sqrt(3)
    points=[]
    for z in range(cells):
        for layer in range(6):
            for y in range(ny):
                for x in range(nx):
                    for f in BASIS[2*layer:2*layer+2]:
                        points.append({'position':[(x+f[0])*lx,(y+f[1])*ly,(z+f[2])*lz],
                                       'occupied':z==0 and layer<5,'fixed':z==0 and layer<2})
    j={'schema_version':1,'description':'Diamond (111), all Table II chemical channels; kinetic/bond assumptions documented in doc/diamond-cvd-1997.md',
       'model':{'type':'diamond_cvd_1997','parameters':{'initial_passivation':True,'site_volume':a**3/8,'ch_bond_length':1.09,
                  'reference':'Battaile, Srolovitz and Butler, J. Appl. Phys. 82, 6293 (1997)',
                  'surface_carbon_atoms':2*nx*ny,'lattice_constant_angstrom':a}},
       'lattice':{'representation':'atomistic','type':'Custom','parameters':{'cell':[nx*lx,ny*ly,cells*lz],
                   'periodic':[True,True,False],'neighbor_cutoff':a*math.sqrt(3)/4*1.01,'site_volume':a**3/8,'initial_species':'C','sites':points}},
       'species':{'C':{'parameters':{'element':'C'}},'H':{'parameters':{'element':'H'}}},
       'system':{'parameters':{'temperature':temperature,'gas_temperature':gas_temperature,
                              'partial_pressures_torr':{'H':.04,'H2':18,'CH3':.004,'C2H2':acetylene}}},
       'simulation':{'duration':duration,'random_seed':1234,'max_events':events},
       'processes':[],'output':{'log_interval':.1,'lattice_interval':.5,'write_xyz':True}}
    def channel(row,reverse,event,domain='',hydrogens=-1,gas=''):
        A,n,E,H,S=TABLE[row]
        correction=1.016e-5 if row in [2,4,5,7,9] else 1
        j['processes'].append({'id':'r'+str(row)+('_reverse' if reverse else '_forward'),'type':'reaction',
             'event':event,'parameters':{'domain':domain,'hydrogens':hydrogens,'table_row':row},
             'rate':{'model':'thermochemical','parameters':{'A':A,'temperature_power':n,
                'activation_energy_kcal_mol':E,'reaction_enthalpy_kcal_mol':H,'reaction_entropy_cal_mol_K':S,
                'reverse':reverse,'standard_concentration_mol_cm3':correction,'gas_reactant':gas}}})
    channel(1,False,'hydrogen_removal','diamond',gas='H');channel(1,True,'hydrogen_addition','diamond',gas='H2')
    channel(2,False,'hydrogen_addition','diamond',gas='H');channel(2,True,'hydrogen_removal','diamond')
    channel(3,False,'desorb_hydrocarbon','C1',2,'H');channel(3,True,'adsorb_carbon','C1',2,'CH3')
    channel(4,False,'adsorb_carbon','C1',3,'CH3');channel(4,True,'desorb_hydrocarbon','C1',3)
    channel(5,False,'adsorb_carbon','C2',2,'C2H2');channel(5,True,'desorb_hydrocarbon','C2',2)
    channel(6,False,'hydrogen_removal','C1',gas='H');channel(6,True,'hydrogen_addition','C1',gas='H2')
    channel(7,False,'hydrogen_addition','C1',gas='H');channel(7,True,'hydrogen_removal','C1')
    channel(8,False,'hydrogen_removal','C2',gas='H');channel(8,True,'hydrogen_addition','C2',gas='H2')
    channel(9,False,'hydrogen_addition','C2',gas='H');channel(9,True,'hydrogen_removal','C2')
    channel(10,False,'append_methyl','C1',gas='CH3')
    channel(11,False,'dimer_opening');channel(11,True,'dimer_closing');channel(12,False,'carbon_incorporation')
    return j

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--nx',type=int,default=2);p.add_argument('--ny',type=int,default=2);p.add_argument('--cells',type=int,default=6)
    p.add_argument('--temperature',type=float,default=1200);p.add_argument('--gas-temperature',type=float,default=1200)
    p.add_argument('--acetylene-pressure-torr',type=float,default=.04);p.add_argument('--duration',type=float,default=2)
    p.add_argument('--max-events',type=int,default=200000)
    a=p.parse_args()
    if a.nx<2 or a.ny<2 or a.cells<2 or not all(math.isfinite(v) and v>0 for v in [a.temperature,a.gas_temperature,a.duration]) or not math.isfinite(a.acetylene_pressure_torr) or a.acetylene_pressure_torr<0 or a.max_events<0:p.error('invalid geometry, temperatures, duration or event limit')
    if 12*a.nx*a.ny*a.cells>1000000:p.error('at most one million candidate sites')
    if a.output.exists():p.error('output already exists; choose a new input path')
    j=make_case(a.nx,a.ny,a.cells,a.temperature,a.gas_temperature,a.acetylene_pressure_torr,a.max_events,a.duration)
    a.output.write_text(json.dumps(j,indent=2)+'\n');print(a.output)

if __name__=='__main__':main()
