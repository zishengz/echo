

from echo.utils import get_rssd, get_nelect_neu_ase, get_GCDFT_neb
from ase.io import read, write
from ase.db import connect
from ase.calculators.vasp import Vasp
from ase.neb import NEB
import numpy as np
import os


def neb_fix_chg(img_ini, img_fin, n_images=None, band_inp=None, net_charge=0, climb=False, calculator=None, optimizer=None, fmax=0.01, steps=1000):
    print('\nTransition state search')
    if n_images == None:
        rssd_rxn = get_rssd(img_ini, img_fin)
        n_images = max(int(get_rssd(img_ini, img_fin) / 0.8), 3)
        print(f'n_images= {n_images} (from RXN RSSD={rssd_rxn:.4f} A)')
    else:
        print(f'n_images= {n_images}')
    workdirs = [f'img{i+1:03}' for i in range(n_images)]
    print(f'CINEB directories: {workdirs}')

    restart_neb = False
    if band_inp is not None:
        band = band_inp.copy()
        restart_neb = True
        print('Starting from provided band')
    else:
        try:
            band = read('neb.traj', index=f'{-n_images-2}:')
            restart_neb = True
            print('Restarting band from neb.traj')
        except:
            band = [img_ini]
            band += [img_ini.copy() for i in range(n_images)]
            band += [img_fin]
            print('Band will be initialized by IDPP interpolation')
    
    # Set separate calculators for each image
    for i in range(len(band)-2):
        band[i+1].calc = Vasp()
        band[i+1].calc.fromdict(calculator.asdict())
        band[i+1].calc.set(charge=net_charge, directory=workdirs[i], txt='vasp.out', lreal='Auto')

    cineb = NEB(band, climb=climb)
    if not restart_neb:
        cineb.interpolate(method='idpp')
    optimizer = optimizer(cineb, trajectory='neb.traj', logfile='-')
    optimizer.run(fmax=fmax, steps=steps)
    print('CONVERGED')
    write('neb-conv.db', band)
    return band


def neb_fix_pot(img_ini, img_fin, n_images=None, band_inp=None, chg_inp=None, pot_ini=0, pot_fin=0, gcfe_ini=0, gcfe_fin=0, nelect_net_ini=0, nelect_net_fin=0, pot_target=0.0, pot_ref=4.6, pot_step=None, geom_step=1, climb=False, calculator=None, optimizer=None, fmax=0.01, pot_conv=0.001, geom_boost=3):
    print(
    f'Ini state: NELECT_net= {nelect_net_ini:8.4f} |e|;  U_she= {pot_ini:8.4f} V; GCFE_el= {gcfe_ini:12.4f} eV')
    print(
    f'Fin state: NELECT_net= {nelect_net_fin:8.4f} |e|;  U_she= {pot_fin:8.4f} V; GCFE_el= {gcfe_fin:12.4f} eV')
    print(f'Target potential: {pot_target} V (SHE refernce: {pot_ref} V)')
    
    # Determine number of images along the band
    if n_images == None:
        rssd_rxn = get_rssd(img_ini, img_fin)
        n_images = max(int(get_rssd(img_ini, img_fin) / 0.8), 3)
        print(f'n_images= {n_images} (from RXN RSSD={rssd_rxn:.4f} A)')
    else:
        print(f'n_images= {n_images}')
    workdirs = [f'img{i+1:03}' for i in range(n_images)]
    print(f'CINEB directories: {workdirs}')

    # Geometry along the band
    restart_band = False
    if band_inp is None:
        try:
            # Read the restart.db
            band = [img_ini]+read('restart.db', index=f'1:{n_images+1}')+[img_fin]
            print('Geometry: Restarting from the last band in restart.db')
            restart_band = True
        except:
            # Start from scratch
            band = [img_ini]
            band += [img_ini.copy() for i in range(n_images)]
            band += [img_fin]
            print('Geometry: Starting from scratch (IDPP interpolation)')
    else:
        print('Geometry: Starting from provided band')
        # Use the provided band geometry (excluding IS/FS)
        if len(band_inp) == n_images:
            band = [img_ini] + band_inp + [img_fin]
        elif len(band_inp) == n_images+2:
            band = [img_ini] + band_inp[1:n_images+1] + [img_fin]
        else:
            print('Input band has wrong size! Check it pls')
            exit()

    
    # Get neutrla state charge
    geom_tmp = img_ini.copy()
    geom_tmp.calc = Vasp()
    geom_tmp.calc.fromdict(calculator.todict())
    nelect_neu = get_nelect_neu_ase(geom_tmp)

    # Charge state along the band
    if restart_band:
        nelect_net_guess = []
        with connect('restart.db') as db:
            for r in db.select():
                try:
                    nelect_net_guess.append(r.nelect_net)
                except:
                    nelect_net_guess.append(-r.charge_net)
        print(f'Charge states: Restarting from restart.db')
    else:
        try:
            nelect_net_guess = np.loadtxt('restart.nelect_net')
            print(f'Charge states: Restarting from restart.nelect_net')
        except:
            nelect_net_guess = np.linspace(
                nelect_net_ini, nelect_net_fin, n_images+2)[1:-1]
            print(f'Charge states: Restarting from scratch (linear CT interpolation)')
        print(f'Net charges along the band: {nelect_net_guess}')

    # Set CINEB band!
    cineb = NEB(band, climb=climb)
    if not restart_band and band_inp is None:
        cineb.interpolate(method='idpp')
    write('band_init.db', band, append=False)
    rssd_img = [get_rssd(band[0], band[m])
                for m in range(n_images+2)]
    print(
        f'Reaction coordinate (RSSD): {[round(n, 4) for n in rssd_img]} Å')    

    # Set separate calculators for each image
    print('\nTransition state search')
    for i in range(n_images):
        band[i+1].calc = Vasp()
        band[i+1].calc.fromdict(calculator.asdict())
        band[i+1].calc.set(charge=-nelect_net_guess[i], directory=workdirs[i], txt='vasp.out', lreal='Auto')

    # Starting the run!
    print(f'\nFIX-POTENTIAL NEB (pot_target= {pot_target}, fmax={fmax}, pot_conv= {pot_conv})')
    if pot_step is None:
        # choose potentionstat step
        pot_step = 0.002 * nelect_neu
        print(f'Potentiostat step set to {pot_step:.6f} |e|/V (NELECT_neu={nelect_neu:8.3f} |e|)')

    traj_nelect = []
    traj_pot = []
    traj_gcfe = []
    nsteps_opt_now = geom_step
    converged_forces = False
    converged_pot = False
    while not (converged_forces and converged_pot):
        optimizer = optimizer(cineb,
                        trajectory='neb.traj', logfile='neb.log')
        optimizer.run(fmax=fmax, steps=nsteps_opt_now)
        nelect_net_now, pot_now, energy_gcdft_now = get_GCDFT_neb(workdirs, pot_ref)

        traj_nelect.append(nelect_net_now)
        traj_pot.append(pot_now)
        traj_gcfe.append(energy_gcdft_now)

        # Write restart files
        with connect('restart.db', append=False) as db:
            for ii in range(n_images):
                if ii == 0 or ii == n_images-1:
                    db.write(band[ii], 
                        charge_net=-nelect_net_ini, nelect_net=nelect_net_ini,
                        pot=pot_ini,
                        gcfe=gcfe_ini)
                elif ii == n_images-1:
                    db.write(band[ii],
                        charge_net=-nelect_net_fin,
                        nelect_net=nelect_net_fin,
                        pot=pot_fin,
                        gcfe=gcfe_fin)
                else:
                    db.write(band[ii],
                        charge_net=-nelect_net_now[ii-1],
                        nelect_net=nelect_net_now[ii-1],
                        pot=pot_now[ii-1], gcfe=energy_gcdft_now[ii-1],
                        potene=band[ii].get_potential_energy())
        np.savetxt('restart.nelect_net', np.array(nelect_net_now))

        print(
            f'\nStep {len(traj_nelect):4}: NELECT_net: {[round(n,4) for n in nelect_net_now]}')
        print(f'Step {len(traj_nelect):4}: U_she: {[round(n,4) for n in pot_now]}')
        print(
            f'Step {len(traj_nelect):4}: GCFE_el: {[round(n,4) for n in energy_gcdft_now]}')

        # Convergence check and optimizer adjusting
        print(f'Conv {len(traj_nelect):4}: ', end='')
        u_diff = (np.abs(np.array(pot_now) - pot_target)).max()
        if u_diff < pot_conv:
            converged_pot = True
            nsteps_opt_now = geom_boost * geom_step
            print(f'Potential [o] U_diff= {u_diff:8.4f} V | ', end='')
        else:
            converged_pot = False
            nsteps_opt_now = geom_step
            print(f'Potential [x] U_diff= {u_diff:8.4f} V | ', end='')

        #fmax = cineb.get_residual()
        fmax_now = (np.sqrt((cineb.get_forces()**2).sum(axis=1))).max()
        if fmax_now < fmax:
            converged_forces = True
            print(f'Forces [o] F_max= {fmax_now:8.4f} eV/Å')
        else:
            nsteps_opt_now *= geom_boost
            converged_forces = False
            print(f'Forces [x] F_max= {fmax_now:8.4f} eV/Å')

        if converged_forces and converged_pot:
            print('CONVERGED')
            break
        else:
            for i in range(len(band)-2):
                nelect_net_new = nelect_net_now[i] - (pot_target-pot_now[i])*pot_step
                cineb.images[i+1].calc.set(charge=-nelect_net_new)
    
    os.system('cp restart.db neb_conv.db')
    return band

