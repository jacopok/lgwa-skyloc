import matplotlib.pyplot as plt
from lgwa_skyloc.psds_against_bns_signal import plot_snr
from lgwa_skyloc.white_dwarfs import merger_gw_frequency, u
from lgwa_skyloc.plotting import FIGS
import numpy as np

if __name__ == '__main__':
    parameters = {
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

    plot_snr(parameters, signal_name = 'GW170817', waveform_model='mlgw_bns')
    plt.savefig(FIGS / 'ET_LGWA_BNS_170817.pdf')

    parameters = {
        'redshift': .01,
        'mass_1': 36.2,
        'mass_2': 29.1,
        'luminosity_distance': 410,
        'theta_jn': 5/6*np.pi,
        'ra': 7/24 * 2 * np.pi,
        'dec': -75 / 180 * np.pi / 2,
        'psi': 1.6,
        'phase': 0 ,
        'geocent_time': 1126259462, # GW150914 time
        'a_1': -.1,
        'a_2': -.1,
    }

    plot_snr(parameters, signal_name = 'GW150914', waveform_model='lalsim_IMRPhenomD')
    plt.savefig(FIGS / 'ET_LGWA_BBH_150914.pdf')

    parameters = {
        'redshift': 1,
        'mass_1': 1e4,
        'mass_2': 1e4,
        'luminosity_distance': 6791,
        'theta_jn': 5/6*np.pi,
        'ra': 7/24 * 2 * np.pi,
        'dec': -75 / 180 * np.pi / 2,
        'psi': 1.6,
        'phase': 0 ,
        'geocent_time': 1126259462, # GW150914 time
        'a_1': -.1,
        'a_2': -.1,
    }

    plot_snr(parameters, signal_name = 'IMBH', waveform_model='lalsim_IMRPhenomD')
    plt.savefig(FIGS / 'ET_LGWA_BBH_IMBH.pdf')

    parameters = {
        'redshift': .005,
        'mass_1': .8,
        'mass_2': .8,
        'luminosity_distance': 22,
        'theta_jn': 5/6*np.pi,
        'ra': 7/24 * 2 * np.pi,
        'dec': -75 / 180 * np.pi / 2,
        'psi': 1.6,
        'phase': 0 ,
        'geocent_time': 1187008882, # GW170817 time
        'a_1': -.1,
        'a_2': -.1,
    }

    plot_snr(parameters, signal_name = 'BWD', waveform_model='lalsim_IMRPhenomD')
    # plt.axvline(merger_gw_frequency(.8*u.Msun, .8*u.Msun).value)    
    plt.savefig(FIGS / 'ET_LGWA_BBH_BWD.pdf')
