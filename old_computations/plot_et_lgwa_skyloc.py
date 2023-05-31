#! /usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from lgwa_skyloc.skyloc_against_time import compute_fisher

from lgwa_skyloc.plotting import CMAP, FIGS

def time(f):
    return (f / 20)**(-8/3) * 180

if __name__ == '__main__':

    params = {
        'mass_1': 1.4, 
        'mass_2': 1.4, 
        'redshift': 0.01,
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
    
    f_low = np.geomspace(.5, 1.6, num=200)
    f_high = np.geomspace(1.6, 10, num=20)
    
    all_f = np.concatenate([f_high])
        
    sls = []
    snrs = []
    
    sls_nolgwa = []
    snrs_nolgwa = []

    for f in f_high:
        sl, pe, snr = compute_fisher(f, params, detectors_ids=['LGWA', 'ET'])
        sls.append(sl* (180/np.pi)**2)
        snrs.append(snr)

    for f in f_high:
        sl, pe, snr = compute_fisher(f, params, detectors_ids=['ET'])
        sls_nolgwa.append(sl* (180/np.pi)**2)
        snrs_nolgwa.append(snr)
    
    plt.loglog(time(all_f) / 3600, sls, color=CMAP(.2), ls=':', label='Skyloc [square degrees], ET+LGWA')
    plt.loglog(time(f_high) / 3600, sls_nolgwa, color=CMAP(.8), ls=':', label='Skyloc [square degrees], ET')

    plt.loglog(time(all_f) / 3600, snrs, color=CMAP(.2), label='SNR, ET+LGWA')
    plt.loglog(time(f_high) / 3600, snrs_nolgwa, color=CMAP(.8), label='SNR, ET')

    FOVs = {
        # 'THESEUS or SVOM': 0.16,
        'Fermi-GBM': .75,
        'Swift-BAT': .11,
        # 'HERMES': 1.
    }
    for name, fov in FOVs.items():
        plt.axhline(fov * (180/np.pi)**2 * 4 * np.pi, lw=.5, color='black')
    plt.axhline(0.7, lw=.5, color='black')
    
    plt.legend()
    plt.xlabel('Time to merger [hours]')
    plt.savefig(FIGS / 'skyloc.pdf')
