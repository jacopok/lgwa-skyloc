import pandas as pd
from GWFish.modules.detection import Network 
from GWFish.modules.fishermatrix import compute_network_errors
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

BASE_PATH = Path(__file__).parent

params = pd.read_hdf(BASE_PATH / 'injections/GWTC3_cosmo_median.hdf5')

z = params.pop('redshift')
params['mass_1'] *= (1+z)
params['mass_2'] *= (1+z)
print(params)

# params = params.append([params]*999)
# plt.scatter(params['ra'], np.cos(np.pi/2 - params['dec']))
# plt.xlim(0, 2*np.pi)
# plt.ylim(-1, 1)
# plt.show()

# rng = np.random.default_rng(seed=1)
# params['ra'] = rng.uniform(0, 2 * np.pi, size=len(params))
# params['dec'] = np.pi/2 - np.arccos(rng.uniform(0, 1, size=len(params)))

# do a 1 parameter Fisher matrix for speed
fisher_params = ['mass_1']

network = Network(['LGWA'], detection_SNR=(0., 0.))

network_snr, parameter_errors, sky_localization = compute_network_errors(
    network,
    params,
    fisher_parameters=fisher_params, 
    waveform_model='IMRPhenomXPHM'
)

bins = np.arange(0, int(max(network_snr))+1)

detected_snrs = network_snr[network_snr >= 9]
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
})
plt.hist(network_snr, bins=bins, alpha=.5, color='tab:blue')
plt.hist(detected_snrs, bins=bins, alpha=.5, color='tab:blue')
plt.xlabel('LGWA SNR')
plt.ylabel('Number of sources')
# plt.title(f'GWTC-3 detected sources = {len(detected_snrs)}/{len(network_snr)} $\\approx$ {len(detected_snrs) /len(network_snr):.0%}%')
plt.savefig('lgwa_gwtc3.pdf', dpi=100)
# plt.show()

