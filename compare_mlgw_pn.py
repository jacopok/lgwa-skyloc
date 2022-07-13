import numpy as np
import matplotlib.pyplot as plt
from lgwa_skyloc.skyloc_against_time import compute_fisher, DEFAULT_FISHER_PARAMETERS

RNG = np.random.default_rng()

def random_bns_params(rng=RNG):
    
    d = rng.uniform(20., 100.)
    return {
        'mass_1': rng.uniform(1., 1.8),
        'mass_2': rng.uniform(1., 1.8),
        'redshift': (70/300_000) * d,
        'luminosity_distance': d,
        'theta_jn': np.arccos(rng.uniform(-1., 1.)),
        'ra': rng.uniform(0, 2. * np.pi),
        'dec': np.arccos(rng.uniform(-1., 1.)) - np.pi / 2.,
        'psi': rng.uniform(0, 2. * np.pi),
        'phase': rng.uniform(0, 2. * np.pi),
        'geocent_time': 1187008882 + rng.uniform(0, 3600*24*365),
        'lambda_1': rng.uniform(10., 1000.),
        'lambda_2': rng.uniform(10., 1000.),
        'a_1': rng.uniform(-.1, .1),
        'a_2': rng.uniform(-.1, .1),
    }


if __name__ == '__main__':

    N = 200
    
    ratios = np.zeros((N, len(DEFAULT_FISHER_PARAMETERS)+2))
    for i in range(N):
        params = random_bns_params()
        sltf2, petf2, snrtf2 = compute_fisher(1024., params, detectors_ids=['ET'], waveform_model='lalsim_IMRPhenomD_NRTidal')
        x_tf2 = np.concatenate(([snrtf2], [sltf2], petf2))
        slmb, pemb, snrmb = compute_fisher(1024., params, detectors_ids=['ET'], waveform_model='mlgw_bns')
        x_mb = np.concatenate(([snrmb], [slmb], pemb))
        
        ratios[i] = x_tf2 / x_mb
        print(x_tf2 / x_mb)
    
    names = ['SNR', 'skyloc'] + DEFAULT_FISHER_PARAMETERS
    
    for i, name in enumerate(names):
        plt.hist(np.log(ratios[:, i]), bins=30)
        plt.title(name)
        plt.xlabel('log (sigma_1 / sigma_2)')
        plt.savefig(f'figs/param_comparison/{name}.png')
        plt.close()
