
from ase.io import read, write
from ase.db import connect
from echo.utils import get_nelect_neu_pp, get_GCDFT
import os

# functions for interfacing with vasp

def get_bluemoon():
    # read from a single-step REPORT file with IBLUEOUT=T
    # returns: cv, lambda, |z|^(-1/2), GkT, |z|^(-1/2)*(lambda+GkT), T_inst
    data_raw = open('REPORT').readlines()
    return [eval(l.split()[2]) for l in data_raw if 'cc>' in l][0],\
        [eval(l.split()[1]) for l in data_raw if 'b_m>' in l][0],\
        [eval(l.split()[2]) for l in data_raw if 'b_m>' in l][0],\
        [eval(l.split()[3]) for l in data_raw if 'b_m>' in l][0],\
        [eval(l.split()[4]) for l in data_raw if 'b_m>' in l][0],\
        [eval(l.split()[2]) for l in data_raw if 'tmprt>' in l][0]

def slow_growth(
    atoms,
    nsteps = 10,
    pot_target=0.0,
    pot_ref=4.6,
    pot_step=None,
):
    write('ini.vasp', atoms)
    nelect_neu = get_nelect_neu_pp(atoms)
    if pot_step is None:
        pot_step = 0.002 * nelect_neu
        print(f'Potentiostat step set to {pot_step:.6f} |e|/V (NELECT_neu={nelect_neu:8.3f} |e|)')
    nelect_net_now = 0

    for i in range(nsteps):
        os.system('$VASP_COMMAND > vasp.out')
        # GCDFT data extraction
        nelect_net_now, pot_now, energy_gcdft_now = get_GCDFT(pot_ref, dirName='.')
        nelect_net_new = nelect_net_now - (pot_target-pot_now)*pot_step
        # Blue moon ensemble data extraction
        os.system('cat REPORT >> bluemoon.txt')
        cv, bm1, bm2, bm3, bm4, temp = get_bluemoon() 
        # write to a trajectory db
        geom_fin = read('vasprun.xml')
        with connect(f'sg_traj.db', append=True) as db:
            db.write(
                geom_fin,
                charge_net=-nelect_net_now,
                nelect_net = nelect_net_now,
                nelect = nelect_net_now + nelect_neu,
                pot=pot_now,
                gcfe=energy_gcdft_now,
                potene=geom_fin.get_potential_energy(),
                temperature=temp,
                z=bm2**-2,
                cv=cv,
                fe_grad=bm4/bm2,
            )
        print(f'cv= {cv:.4f} T= {temp} pot= {pot_now:.3f} chg_net {-nelect_net_now:.3f} gcfe= {energy_gcdft_now:.3f} epot= {geom_fin.get_potential_energy():.3f} z= {bm2**-2:.3f} fe_grad= {bm4/bm2:.3f}')
        # update the NELECT in INCAR (1st line)
        incar = open('INCAR').readlines()
        incar[0] = f'NELECT={nelect_net_new + nelect_neu}\n'
        with open('INCAR', 'w') as f:
            for l in incar:
                f.write(l)
        # update the POSCAR input structure
        os.system('cp CONTCAR POSCAR')
        
    




