import numpy as np
import astropy.units as u
import astropy.constants as ac
import matplotlib.pyplot as plt
from astropy import visualization

def wd_radius(mass: u.Quantity):
    return 0.012 * u.Rsun * (
        + (mass / 1.44 / u.Msun)**(-2/3)
        - (mass / 1.44 / u.Msun)**(2/3)
    )**(1/2)

def merger_gw_frequency_hz(mass_1: u.Quantity, mass_2: u.Quantity):
    return 2 * np.sqrt(
        ac.G * (mass_1 + mass_2) / 4 / np.pi**2 
        / ((wd_radius(mass_1) + wd_radius(mass_2)) / 2)**3
    ).si.value

def simulate_bwd(rng, size):
    mass_1 = rng.normal(loc=.6, scale=.1, size=size) * u.Msun
    mass_2 = rng.normal(loc=.6, scale=.1, size=size) * u.Msun

    return mass_1, mass_2

if __name__ == '__main__':
    
    rng = np.random.default_rng()
    
    freqs = merger_gw_frequency_hz(*simulate_bwd(rng, size=100_000))
    
    with visualization.quantity_support():
        plt.hist(freqs, bins=100, density=True)
    plt.ylabel('Probability density')
    plt.xlabel('Merger frequency [Hz]')
    plt.savefig('bwd_merger_frequencies.pdf')