from echo.utils import get_GCDFT, get_nelect_neu_pp, get_nelect_incar 
from ase.db import connect
from ase.io import read, write
from ase.calculators.vasp import Vasp
import numpy as np
import os
import json


def sp_fix_pot(atoms, label='tmp', pot_target=0.0, pot_ref=4.6, pot_step=None, geom_step=1, calculator=None, pot_conv=0.001):
    if f'{label}_conv.db' in os.listdir() and f'{label}_gcdft.json' in os.listdir():
        print('Converged already')
        return read(f'{label}_conv.db'), json.load(open(f'{label}_gcdft.json')) 
    
    atoms_opt = atoms.copy()

    # # Re-initialize calculator
    # mycalc = Vasp()
    # mycalc.fromdict(calculator.asdict())
    mycalc = calculator
    mycalc.set(directory=label, txt='vasp.out')
    print(f'VASP settings:\n{mycalc.todict()}')

    atoms_opt.set_calculator(mycalc)
    nelect_neu = get_nelect_neu_pp(atoms_opt)
    try:
        nelect_now = get_nelect_incar(label)
        charge_now = nelect_neu-nelect_now
        print(f'Charge state read from last INCAR: {charge_now:.4f} |e|')
    except:
        charge_now = 0
        print('Starting from zero net charge...')

    print(f'\nFIX-POTENTIAL Single Point (pot_target= {pot_target}, pot_conv= {pot_conv})')
    if pot_step is None:
        pot_step = 0.002 * nelect_neu
        print(f'Potentiostat step set to {pot_step:.6f} |e|/V (NELECT_neu={nelect_neu:8.3f} |e|)')
    
    # Set calculator and start the run!
    atoms_opt.calc.set(charge=charge_now)
    traj_nelect = []
    traj_pot = []
    traj_gcfe = []
    nsteps_opt_now = geom_step
    converged_forces = False
    converged_pot = False
    while not converged_pot:
        sp_tmp = atoms_opt.get_potential_energy()
        nelect_net_now, pot_now, energy_gcdft_now = get_GCDFT(pot_ref, dirName=label)
        traj_nelect.append(nelect_net_now)
        traj_pot.append(pot_now)
        traj_gcfe.append(energy_gcdft_now)
        np.savetxt(f'{label}/log_nelect.txt', np.array(traj_nelect))
        np.savetxt(f'{label}/log_potential.txt', np.array(traj_pot))
        np.savetxt(f'{label}/log_gcfe.txt', np.array(traj_gcfe))
        print(f'Step {len(traj_nelect):4}: NELECT_net= {traj_nelect[-1]:8.4f} |e|;  U_she= {traj_pot[-1]:8.4f} V; GCFE_el= {traj_gcfe[-1]:12.4f} eV')

        # Convergence check and optimizer adjusting
        print(f'Conv {len(traj_nelect):4}: ', end='')
        u_diff = pot_now - pot_target
        if np.abs(u_diff) < pot_conv:
            converged_pot=True
            print(f'Potential [o] U_diff= {u_diff:8.4f} V | ', end='')
        else:
            converged_pot=False
            print(f'Potential [x] U_diff= {u_diff:8.4f} V | ', end='')
        
        if converged_pot:
            print('CONVERGED')
            break
        else:
            nelect_net_new = nelect_net_now - (pot_target-pot_now)*pot_step
            atoms_opt.calc.set(charge=-nelect_net_new)

    with connect(f'{label}_conv.db', append=False) as db:
        db.write(
            atoms_opt,
            charge_net=-traj_nelect[-1],
            nelect_net = traj_nelect[-1],
            pot=traj_pot[-1],
            gcfe=traj_gcfe[-1],
            potene=atoms_opt.get_potential_energy()
        )
    return atoms_opt
    


