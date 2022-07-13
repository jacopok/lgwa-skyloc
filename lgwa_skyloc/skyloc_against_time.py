import numpy as np
import pandas as pd

import GWFish.modules as gw

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

fisher_parameters = [
    'geocent_time',
    'ra', 
    'dec', 
    'psi', 
    'theta_jn', 
    'luminosity_distance',
    'mass_1',
    'mass_2',
    'phase',
    # 'lambda_1',
    # 'lambda_2',
    # 'a_1',
    # 'a_2'
]
threshold_SNR = np.array([0., 0.1])
duty_cycle = False  # whether to consider the duty cycle of detectors

def compute_fisher(fmax, params, detectors_ids = ['ET', 'LGWA']):
    
    networks_ids = [list(range(len(detectors_ids)))]

    parameters = pd.DataFrame.from_dict(
        {k: v*np.ones((1,)) for k, v in params.items()}
    )
    
    network = gw.detection.Network(detectors_ids, detection_SNR=threshold_SNR, parameters=parameters,
                                   fisher_parameters=fisher_parameters, fmax=fmax)

    waveform_model = 'gwfish_TaylorF2'
    # waveform_model = 'mlgw_bns'

    parameter_values = parameters.iloc[0]

    networkSNR_sq = 0
    for det in network.detectors:
        wave, t_of_f = gw.waveforms.hphc_amplitudes(waveform_model, parameter_values,
                                                    det.frequencyvector)
                                                    #plot=det.plotrange)
        signal = gw.detection.projection(parameter_values, det, wave, t_of_f)

        SNRs = gw.detection.SNR(det, signal, duty_cycle=duty_cycle)
        networkSNR_sq += np.sum(SNRs ** 2)
        det.SNR[0] = np.sqrt(np.sum(SNRs ** 2))

        det.fisher_matrix[0, :, :] = \
            gw.fishermatrix.FisherMatrix(waveform_model, parameter_values, fisher_parameters, det)

    network.SNR[0] = np.sqrt(networkSNR_sq)

    gw.detection.analyzeDetections(network, parameters, 'single_binary', networks_ids)
    sl, pe = gw.fishermatrix.analyzeFisherErrors(network, parameters, fisher_parameters, 'single_binary', networks_ids)

    return sl[0], pe[0], network.SNR[0]

