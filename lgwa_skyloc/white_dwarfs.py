import numpy as np
import astropy.units as u
import astropy.constants as ac
import matplotlib.pyplot as plt
from astropy import visualization
from pathlib import Path

FIGS = Path(__file__).parent.parent / 'figs'

def wd_radius(mass: u.Quantity):
    return 0.012 * u.Rsun * (
        + (mass / 1.44 / u.Msun)**(-2/3)
        - (mass / 1.44 / u.Msun)**(2/3)
    )**(1/2)

def merger_gw_frequency(mass_1: u.Quantity, mass_2: u.Quantity):
    return 2 * np.sqrt(
        ac.G * (mass_1 + mass_2) / 4 / np.pi**2 
        / ((wd_radius(mass_1) + wd_radius(mass_2)) / 2)**3
    ).si

def simulate_bwd(rng, size):
    mass_1 = rng.normal(loc=.6, scale=.1, size=size) * u.Msun
    mass_2 = rng.normal(loc=.6, scale=.1, size=size) * u.Msun

    return mass_1, mass_2

if __name__ == '__main__':
    
    rng = np.random.default_rng()
    
    freqs = merger_gw_frequency(*simulate_bwd(rng, size=100_000))
    
    with visualization.quantity_support():
        plt.hist(freqs, bins=100, density=True)
    plt.ylabel('Probability density')
    plt.xlabel('Merger frequency [Hz]')
    plt.savefig(FIGS / 'bwd_merger_frequencies.pdf')
    plt.close()
    
    masses = np.geomspace(0.1, 1.2, num=400) * u.Msun
    with visualization.quantity_support():
        plt.plot(masses, merger_gw_frequency(masses, masses))
    plt.ylabel('Merger frequency [Hz]')
    plt.xlabel('White dwarf mass [$M_{\odot}$]')
    plt.savefig(FIGS / 'bwd_merger_frequencies_by_mass.pdf')
    plt.close()

    with visualization.quantity_support():
        plt.plot(masses, wd_radius(masses))
    plt.savefig(FIGS / 'wd_radii_by_mass.pdf')
    plt.close()

    with visualization.quantity_support():
        plt.loglog(masses, (wd_radius(masses) / masses / ac.G * ac.c**2).si )
    plt.savefig(FIGS / 'wd_compactness_by_mass.pdf')
    plt.close()
