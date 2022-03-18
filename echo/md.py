from echo.utils import get_GCDFT, get_fmax, get_nelect_neu_pp, get_nelect_incar 
from ase.calculators.vasp import Vasp
from ase.io import read, write
from ase.db import connect
from ase.optimize import FIRE
import numpy as np
from ase import units

from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution


def langevin_fixpot(
    atoms=None,
    timestep_fs=1,
    temperature_K=300,
    nsteps = 1000,
    friction=0.1,
    label='',
    pot_target=0.0,
    pot_ref=4.6,
    pot_step=None,
    potentiostat_nsteps=10,
    calculator=None,
):
    mycalc = calculator
    mycalc.set(directory=label, txt='vasp.out')
    print(f'VASP settings:\n{mycalc.todict()}')

    # Read geometry if trajectory file exists
    restart_traj = False
    fn_traj = f'{label}.traj'
    try:
        geom_md = read(fn_traj)
        print(f'Reading geometry from {fn_traj}...')
        restart_traj = True
    except:
        geom_md = atoms.copy()
        geom_md.set_calculator(mycalc)
        print('Starting from input geometry...')

    # Get charge state and set calculator
    if restart_traj and 'charge' in geom_md.calc.todict().keys():
        charge_now = geom_md.calc.todict()['charge']
        geom_md.set_calculator(mycalc)
        nelect_neu = get_nelect_neu_pp(geom_md)
        print(f'Charge state read from trajectory: {charge_now:.4f} |e|')
    else:
        geom_md.set_calculator(mycalc)
        nelect_neu = get_nelect_neu_pp(geom_md)
        try:
            nelect_now = get_nelect_incar(label)
            charge_now = nelect_neu-nelect_now
            print(f'Charge state read from last INCAR: {charge_now:.4f} |e|')
        except:
            charge_now = 0
            print('Starting from zero net charge...')

    dyn = Langevin(
        atoms=geom_md,
        timestep=timestep_fs * units.fs,
        temperature_K=temperature_K,
        friction=friction,
        logfile=f'{label}.log',
        trajectory=f'{label}.traj',
    )

    print(f'\nFIX-POTENTIAL BOMD (pot_target= {pot_target}, timestep={timestep_fs} fs, potentialstat per {potentiostat_nsteps} fs)')
    if pot_step is None:
        pot_step = 0.002 * nelect_neu
        print(f'Potentiostat step set to {pot_step:.6f} |e|/V (NELECT_neu={nelect_neu:8.3f} |e|)')

    geom_md.calc.set(charge=charge_now)
    traj_nelect = []
    traj_pot = []
    traj_gcfe = []
    nsteps_done = 0
    while nsteps_done < nsteps:
        dyn.run(potentiostat_nsteps)
        nsteps_done += potentiostat_nsteps
        nelect_net_now, pot_now, energy_gcdft_now = get_GCDFT(pot_ref, dirName=label)
        traj_nelect.append(nelect_net_now)
        traj_pot.append(pot_now)
        traj_gcfe.append(energy_gcdft_now)
        np.savetxt(f'{label}/log_nelect.txt', np.array(traj_nelect))
        np.savetxt(f'{label}/log_potential.txt', np.array(traj_pot))
        np.savetxt(f'{label}/log_gcfe.txt', np.array(traj_gcfe))
        print(f'Step {len(traj_nelect):4}: NELECT_net= {traj_nelect[-1]:8.4f} |e|;  U_she= {traj_pot[-1]:8.4f} V; GCFE_el= {traj_gcfe[-1]:12.4f} eV')
        with connect(f'{label}_traj.db', append=True) as db:
            db.write(
                geom_md,
                charge_net=-traj_nelect[-1],
                nelect_net = traj_nelect[-1],
                pot=traj_pot[-1],
                gcfe=traj_gcfe[-1],
                potene=geom_md.get_potential_energy()
            )
        if nsteps_done >= nsteps:
            break
        else:
            nelect_net_new = nelect_net_now - (pot_target-pot_now)*pot_step
            geom_md.calc.set(charge=-nelect_net_new)
    return geom_md
