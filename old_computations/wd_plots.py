import numpy as np
import astropy.units as u
import astropy.constants as ac
import matplotlib.pyplot as plt
from matplotlib import ticker
from astropy import visualization
from pathlib import Path

from lgwa_skyloc.plotting import FIGS, CMAP
from lgwa_skyloc.population_time_to_merger import time_to_merger
from lgwa_skyloc.white_dwarfs import merger_gw_frequency, wd_radius

PLOT_NUM = 150

def merger_freq_plot(
    mass_range = (.2, 1.2),
    freq_range = (1e-2, 1.),
    fig_name = 'bwd_merger_frequencies_by_mass',
    log_mass = False):

    if log_mass:
        masses = np.geomspace(*mass_range, num=PLOT_NUM) * u.Msun
    else:
        masses = np.linspace(*mass_range, num=PLOT_NUM) * u.Msun

    plt.plot(merger_gw_frequency(masses, masses), masses, c='white', lw=3)
    plt.xscale('log')
    if log_mass:
        plt.yscale('log')

    freqs = np.geomspace(*freq_range, num=PLOT_NUM)
    F, M = np.meshgrid(freqs, masses.value)
    T = np.maximum(time_to_merger(F, M, M) / 3.154e+7, 1./365./24.)
    plt.contourf(F, M, T,
                locator=ticker.LogLocator(subs=(1, np.sqrt(10))),
                cmap=CMAP)
    # plt.quiver(F, M, (1/T)**(1/3), np.zeros_like(T))
    # plt.legend()
    def fmt(x, pos):
        a, b = f'{x:.2e}'.split('e')
        b = int(b)
        return r'$10^{{{}}}$'.format(b)


    plt.colorbar(label='Time to merger [yr]', format=ticker.FuncFormatter(fmt))
    plt.xlabel('Frequency [Hz, log]')
    plt.ylabel('Star mass [$M_{\odot}$]')
    plt.savefig(FIGS / (fig_name + '.pdf'))
    plt.close()

if __name__ == '__main__':
    merger_freq_plot()
    merger_freq_plot(
        mass_range=(.2, 100.),
        freq_range=(1e-2, 1e4), 
        fig_name='merger_freq_plot_wd_to_bh', 
        log_mass=True
    )