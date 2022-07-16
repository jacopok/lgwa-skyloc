import numpy as np
import astropy.units as u
import astropy.constants as ac
import matplotlib.pyplot as plt
from astropy import visualization
from pathlib import Path

from .plotting import FIGS

NS_RADIUS = (12 * u.km).to(u.Rsun)

def bh_radius(mass: u.Quantity):
    return (mass * 2 * ac.G / ac.c**2).to(u.Rsun)

def wd_radius(mass: u.Quantity):
    return np.maximum(
        0.012 * u.Rsun * np.nan_to_num((
        + (mass / 1.44 / u.Msun)**(-2/3)
        - (mass / 1.44 / u.Msun)**(2/3)
    )**(1/2), nan=0.), np.maximum(NS_RADIUS, bh_radius(mass)))

def merger_gw_frequency(mass_1: u.Quantity, mass_2: u.Quantity):
    return 2 * np.sqrt(
        ac.G * (mass_1 + mass_2) / 4 / np.pi**2 
        / ((wd_radius(mass_1) + wd_radius(mass_2)) / 2)**3
    ).to(u.Hz)

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

    # with visualization.quantity_support():
    #     plt.plot(masses, wd_radius(masses))
    # plt.savefig(FIGS / 'wd_radii_by_mass.pdf')
    # plt.close()

    # with visualization.quantity_support():
    #     plt.loglog(masses, (wd_radius(masses) / masses / ac.G * ac.c**2).si )
    # plt.savefig(FIGS / 'wd_compactness_by_mass.pdf')
    # plt.close()
