import matplotlib.pyplot as plt
from lgwa_skyloc.psds_against_bns_signal import plot_snr
from lgwa_skyloc.white_dwarfs import FIGS
import numpy as np

if __name__ == '__main__':
    parameters = {
        'phase': 0.,
        # 'geocent_time': 1577491218., # 1st of January 2030
        'redshift': .01,
        'mass_1': 1.4,
        'mass_2': 1.4,
        'luminosity_distance': 40,
        'theta_jn': 5/6*np.pi,
        'ra': 3.45,
        'dec': -0.41,
        'psi': 1.6 ,
        'phase': 0 ,
        'geocent_time': 1187008882, # GW170817 time
        'lambda_1': 400,
        'lambda_2': 400,
        'a_1': 0,
        'a_2': 0,

    }

    plot_snr(parameters)
    plt.savefig(FIGS / 'ET_LGWA_BNS.pdf')