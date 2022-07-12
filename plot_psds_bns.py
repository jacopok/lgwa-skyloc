import matplotlib.pyplot as plt
from lgwa_skyloc.psds_against_bns_signal import plot_snr
from lgwa_skyloc.white_dwarfs import FIGS

if __name__ == '__main__':
    parameters = {
        'phase': 0.,
        'geocent_time': 1577491218., # 1st of January 2030
        'redshift': .01,
        'luminosity_distance': 40,
        'theta_jn': 1.,
        'mass_1': 1.5,
        'mass_2': 1.5,
        'ra': 0.,
        'dec': 0.,
        'psi': 0.,
        'lambda_1': 500., 
        'lambda_2': 500., 
        'a_1': 0., 
        'a_2': 0., 
    }

    plot_snr(parameters)
    plt.savefig(FIGS / 'ET_LGWA_BNS.pdf')