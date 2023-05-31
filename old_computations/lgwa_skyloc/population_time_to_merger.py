import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib import ticker
from pathlib import Path

from .white_dwarfs import merger_gw_frequency
from .plotting import FIGS, CMAP

# this value refers to the eta=1/4, m_tot=2.8 case
# it is the seglen, in seconds, for a signal starting at 20Hz
# computed at the PN level
SEGLEN_20_HZ = 157.86933774

@np.vectorize
def base_time_to_merger(
    f_0: float,
    m_1: float,
    m_2: float,
) -> float:
    r"""
    The seglen has a closed-form expression in the Newtonian limit,
    see e.g. Maggiore (2007), eq. 4.21:

    :math:`t = 5/256 (\pi f)^{-8/3} (G M \eta ^{3/5} / c^3)^{-5/3}`
    """

    mass_ratio = m_1 / m_2
    if mass_ratio < 1:
        mass_ratio = 1 / mass_ratio
    m_tot = m_1 + m_2
    eta = mass_ratio / (1 + mass_ratio) ** 2

    return (
        SEGLEN_20_HZ * (f_0 / 20) ** (-8 / 3) * (m_tot / 2.8) ** (-5 / 3) / (4 * eta)
    )

@np.vectorize
def time_to_merger(
    f_0: float,
    m_1: float,
    m_2: float,

) -> float:
    
    merger_freq = merger_gw_frequency(m_1 * u.Msun, m_2 * u.Msun).value
    
    t0 = base_time_to_merger(merger_freq, m_1, m_2)
    
    return max(base_time_to_merger(f_0, m_1, m_2) - t0, 0)

def plot_time_to_merger(p_min: float, file_name: Path):
    m1_range = np.linspace(.2, 1.2, num=30)
    m2_range = np.linspace(.2, 1.2, num=30)
    
    t = np.zeros((m1_range.shape + m2_range.shape))
    
    for (i, j), time  in np.ndenumerate(t):
        t[i, j] = time_to_merger(2/ p_min / 60, m1_range[i], m2_range[j]) / 3600 / 24 / 365

    plt.contourf(m1_range, m2_range, t, levels=100, cmap=CMAP)
    plt.xlabel(r"Primary mass [$M_{\odot}]$")
    plt.ylabel(r"Secondary mass [$M_{\odot}]$")
    # plt.yscale('log')
    # plt.xscale('log')
    plt.title(f'Time to merge for a WD binary with $P = {p_min}' + r'\ \text{min}$')
    plt.colorbar(label='Time to merge [yr]')
    plt.savefig(file_name)
    plt.close()

if __name__ == '__main__':
    
    plot_time_to_merger(6, FIGS / 'bwd_time_to_merger_six_minutes.pdf')
    
    plot_time_to_merger(.5, FIGS / 'bwd_time_to_merger_half_minute.pdf')
