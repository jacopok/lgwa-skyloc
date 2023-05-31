import numpy as np
import matplotlib.pyplot as plt
from lgwa_skyloc.skyloc_against_time import compute_fisher, DEFAULT_FISHER_PARAMETERS

RNG = np.random.default_rng(seed=1)

fisher_parameters = DEFAULT_FISHER_PARAMETERS + [
    'lambda_1',
    'lambda_2',
    'a_1',
    'a_2',
]

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

    N = 20
    
    w_1 = 'lalsim_IMRPhenomD_NRTidal'
    w_2 = 'mlgw_bns'
    
    import warnings; warnings.filterwarnings("ignore", category=UserWarning)
    
    ratios = np.zeros((N, len(fisher_parameters)+2))
    for i in range(N):
        params = random_bns_params()
        sl1, pe1, snr1 = compute_fisher(
            1024., 
            params, 
            detectors_ids=['ET'], 
            waveform_model=w_1, 
            fisher_parameters=fisher_parameters
        )
        x_1 = np.concatenate(([snr1], [sl1], pe1))

        sl2, pe2, snr2 = compute_fisher(
            1024., 
            params, 
            detectors_ids=['ET'], 
            waveform_model=w_2, 
            fisher_parameters=fisher_parameters
        )
        x_2 = np.concatenate(([snr2], [sl2], pe2))
        
        ratios[i] = x_1 / x_2
        print(x_1 / x_2)
    
    names = ['SNR', 'skyloc'] + fisher_parameters
    
    for i, name in enumerate(names):
        plt.hist(np.log(ratios[:, i]), bins=30)
        plt.title(name)
        plt.xlabel(f'log ({w_1} / {w_2})')
        plt.savefig(f'figs/param_comparison/{name}.png')
        plt.close()
