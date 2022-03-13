from echo.utils import get_GCDFT, get_fmax, get_nelect_neu_ase, get_nelect_incar 
from ase.db import connect
from ase.io import read, write
from ase.calculators.vasp import Vasp
import numpy as np
import os
import json


def opt_fix_chg(atoms, label='tmp', net_charge=0, calculator=None, optimizer=None, fmax=0.01, steps=1000, must_run=False):
    if not must_run and f'{label}_conv.db' in os.listdir():
        print('Converged already')
        return read(f'{label}_conv.db')

    # # Re-initialize calculator
    # mycalc = Vasp()
    # mycalc.fromdict(calculator.asdict())
    mycalc=calculator
    mycalc.set(directory=label, charge=net_charge, txt='vasp.out')
    print(f'VASP settings:\n{mycalc.todict()}')

    # Read geometry if trajectory file exists
    restart_traj = False
    fn_traj = f'{label}.traj'
    try:
        atoms_opt = read(fn_traj)
        print(f'Reading geometry from {fn_traj}...')
        restart_traj=True
    except:
        atoms_opt = atoms.copy()
        atoms_opt.set_calculator()
        print('Starting from input geometry...')

    # End if converged already
    if restart_traj and atoms_opt.calc is not None:
        print('Checking forces from restart file...')
        if get_fmax(atoms_opt.get_forces()) < fmax:
            print('Converged already!')
            return atoms_opt

    
    # Set calculator and optimizer & RUN!
    print(f'\nFIX-CHARGE GEOMOPT (fmax={fmax})')
    atoms_opt.set_calculator(mycalc)
    geomopt = optimizer(atoms_opt, trajectory=fn_traj, logfile='-')
    geomopt.run(fmax=fmax, steps=steps)
    print('CONVERGED')
    write(f'{label}_conv.db', atoms_opt)
    return atoms_opt


def opt_fix_pot(atoms, label='tmp', pot_target=0.0, pot_ref=4.6, pot_step=None, geom_step=1, calculator=None, optimizer=None, fmax=0.01, pot_conv=0.001, geom_boost=3):
    if f'{label}_conv.db' in os.listdir() and f'{label}_gcdft.json' in os.listdir():
        print('Converged already')
        return read(f'{label}_conv.db'), json.load(open(f'{label}_gcdft.json')) 

    # # Re-initialize calculator
    # mycalc = Vasp()
    # mycalc.fromdict(calculator.asdict())
    mycalc = calculator
    mycalc.set(directory=label, txt='vasp.out')
    print(f'VASP settings:\n{mycalc.todict()}')

    # Read geometry if trajectory file exists
    restart_traj = False
    fn_traj = f'{label}.traj'
    try:
        atoms_opt = read(fn_traj)
        print(f'Reading geometry from {fn_traj}...')
        restart_traj = True
    except:
        atoms_opt = atoms.copy()
        atoms_opt.set_calculator(mycalc)
        print('Starting from input geometry...')

    # Get charge state and set calculator
    if restart_traj and 'charge' in atoms_opt.calc.todict().keys():
        charge_now = atoms_opt.calc.todict()['charge']
        atoms_opt.set_calculator(mycalc)
        nelect_neu = get_nelect_neu_ase(atoms_opt)
        print(f'Charge state read from trajectory: {charge_now:.4f} |e|')
    else:
        atoms_opt.set_calculator(mycalc)
        nelect_neu = get_nelect_neu_ase(atoms_opt)
        try:
            nelect_now = get_nelect_incar(label)
            charge_now = nelect_neu-nelect_now
            print(f'Charge state read from last INCAR: {charge_now:.4f} |e|')
        except:
            charge_now = 0
            print('Starting from zero net charge...')

    print(f'\nFIX-POTENTIAL GEOMOPT (pot_target= {pot_target}, fmax={fmax}, pot_conv= {pot_conv})')
    if pot_step is None:
        pot_step = 0.002 * nelect_neu
        print(f'Potentiostat step set to {pot_step:.6f} |e|/V (NELECT_neu={nelect_neu:8.3f} |e|)')
    atoms_opt.calc.set(charge=charge_now)
    traj_nelect = []
    traj_pot = []
    traj_gcfe = []
    nsteps_opt_now = geom_step
    converged_forces = False
    converged_pot = False
    while not (converged_forces and converged_pot):
        geomopt = optimizer(atoms_opt, trajectory=fn_traj, logfile=f'{label}.log')
        geomopt.run(fmax=fmax, steps=nsteps_opt_now)
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
            nsteps_opt_now = geom_boost * geom_step
            print(f'Potential [o] U_diff= {u_diff:8.4f} V | ', end='')
        else:
            converged_pot=False
            nsteps_opt_now = geom_step
            print(f'Potential [x] U_diff= {u_diff:8.4f} V | ', end='')
            
        fmax_now = get_fmax(atoms_opt.get_forces())
        if fmax_now < fmax:
            converged_forces=True
            print(f'Forces [o] F_max= {fmax_now:8.4f} eV/Å')
        else:
            nsteps_opt_now *= geom_boost
            converged_forces=False
            print(f'Forces [x] F_max= {fmax_now:8.4f} eV/Å')
        
        if converged_forces and converged_pot:
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
    