from ase.utils.forcecurve import fit_raw
from ase.db import connect


def fit_images(neb_db):
    R = []
    E = []
    F = []
    A = []
    pbc = []
    for r in neb_db.select():
        R.append(r.toatoms().positions)
        E.append(r.gcfe)
        F.append(r.toatoms().get_forces())
        A = r.toatoms().cell
        pbc = r.toatoms().pbc
    return fit_raw(E, F, R, A, pbc)


class NEBTools:
    """Class to make many of the common tools for NEB analysis available to
    the user. Useful for scripting the output of many jobs. Initialize with
    list of images which make up one or more band of the NEB relaxation."""

    def __init__(self, images):
        self.images = images

    def get_fit(self):
        return fit_images(self.images)

    def get_barrier(self, fit=True, raw=False):
        """Returns the barrier estimate from the NEB, along with the
        Delta E of the elementary reaction. If fit=True, the barrier is
        estimated based on the interpolated fit to the images; if
        fit=False, the barrier is taken as the maximum-energy image
        without interpolation. Set raw=True to get the raw energy of the
        transition state instead of the forward barrier."""
        forcefit = fit_images(self.images)
        energies = forcefit.energies
        fit_energies = forcefit.fit_energies
        dE = energies[-1] - energies[0]
        if fit:
            barrier = max(fit_energies)
        else:
            barrier = max(energies)
        if raw:
            barrier += self.images[0].get_potential_energy()
        return barrier, dE


    def get_fmax(self, **kwargs):
        """Returns fmax, as used by optimizers with NEB."""
        neb = NEB(self.images, **kwargs)
        forces = neb.get_forces()
        return np.sqrt((forces ** 2).sum(axis=1).max())

    def plot_band(self, ax=None):
        """Plots the NEB band on matplotlib axes object 'ax'. If ax=None
        returns a new figure object."""
        forcefit = fit_images(self.images)
        ax = forcefit.plot(ax=ax)
        return ax.figure


    def plot_bands(self, constant_x=False, constant_y=False,
                   nimages=None, label='nebplots'):
        """Given a trajectory containing many steps of a NEB, makes
        plots of each band in the series in a single PDF.

        constant_x: bool
            Use the same x limits on all plots.
        constant_y: bool
            Use the same y limits on all plots.
        nimages: int
            Number of images per band. Guessed if not supplied.
        label: str
            Name for the output file. .pdf will be appended.
        """
        from matplotlib import pyplot
        from matplotlib.backends.backend_pdf import PdfPages
        if nimages is None:
            nimages = self._guess_nimages()
        nebsteps = len(self.images) // nimages
        if constant_x or constant_y:
            sys.stdout.write('Scaling axes.\n')
            sys.stdout.flush()
            # Plot all to one plot, then pull its x and y range.
            fig, ax = pyplot.subplots()
            for index in range(nebsteps):
                images = self.images[index * nimages:(index + 1) * nimages]
                NEBTools(images).plot_band(ax=ax)
                xlim = ax.get_xlim()
                ylim = ax.get_ylim()
            pyplot.close(fig)  # Reference counting "bug" in pyplot.
        with PdfPages(label + '.pdf') as pdf:
            for index in range(nebsteps):
                sys.stdout.write('\rProcessing band {:10d} / {:10d}'
                                 .format(index, nebsteps))
                sys.stdout.flush()
                fig, ax = pyplot.subplots()
                images = self.images[index * nimages:(index + 1) * nimages]
                NEBTools(images).plot_band(ax=ax)
                if constant_x:
                    ax.set_xlim(xlim)
                if constant_y:
                    ax.set_ylim(ylim)
                pdf.savefig(fig)
                pyplot.close(fig)  # Reference counting "bug" in pyplot.
        sys.stdout.write('\n')


    def _guess_nimages(self):
        """Attempts to guess the number of images per band from
        a trajectory, based solely on the repetition of the
        potential energy of images. This should also work for symmetric
        cases."""
        e_first = self.images[0].get_potential_energy()
        nimages = None
        for index, image in enumerate(self.images[1:], start=1):
            e = image.get_potential_energy()
            if e == e_first:
                # Need to check for symmetric case when e_first = e_last.
                try:
                    e_next = self.images[index + 1].get_potential_energy()
                except IndexError:
                    pass
                else:
                    if e_next == e_first:
                        nimages = index + 1  # Symmetric
                        break
                nimages = index  # Normal
                break
        if nimages is None:
            sys.stdout.write('Appears to be only one band in the images.\n')
            return len(self.images)
        # Sanity check that the energies of the last images line up too.
        e_last = self.images[nimages - 1].get_potential_energy()
        e_nextlast = self.images[2 * nimages - 1].get_potential_energy()
        if not (e_last == e_nextlast):
            raise RuntimeError('Could not guess number of images per band.')
        sys.stdout.write('Number of images per band guessed to be {:d}.\n'
                         .format(nimages))
        return nimages



    