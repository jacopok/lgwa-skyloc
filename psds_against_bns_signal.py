import numpy as np
import matplotlib.pyplot as plt
from GWFish.modules.waveforms import hphc_amplitudes
from GWFish.modules import detection

from scipy.interpolate import interp1d

from white_dwarfs import FIGS

class SignalInDetector:
    
    def __init__(self, parameters: dict, detector: detection.Detector):
        self.parameters = parameters
        self.detector = detector

        self.detector.frequencyvector = self.frequencyvector

    
    @property    
    def psd(self):
        return self.detector.components[0].psd_data.T[1]

    @property
    def frequencyvector(self):
        return self.detector.components[0].psd_data.T[0]

    @property
    def projection(self):
        polarizations, timevector = hphc_amplitudes(
            'gwfish_TaylorF2', 
            self.parameters, 
            self.frequencyvector[:, np.newaxis], plot=None)

        projections_set = detection.projection(
            self.parameters, 
            self.detector, 
            polarizations, 
            timevector
        )
        
        return np.sqrt(np.sum(abs(projections_set)**2, axis=1))

    @property
    def characteristic_strain_signal(self):
        return characteristic_strain_from_fourier_amplitude(
                self.frequencyvector, self.projection)

    @property
    def characteristic_strain_noise(self):
        return characteristic_strain_from_psd(
                self.frequencyvector, self.psd
            )

    @property
    def snr_per_efolding(self):
        return self.characteristic_strain_signal / self.characteristic_strain_noise

    @property
    def cumulative_snr(self):
        d_log_f = np.gradient(self.frequencyvector) / self.frequencyvector

        return np.sqrt(np.cumsum(self.snr_per_efolding**2 * d_log_f))


class SignalAcrossDetectors:
    
    def __init__(self, parameters: dict, detectors: list[detection.Detector]):
        self.parameters = parameters
        self.detectors = detectors
        
        self.signals = [SignalInDetector(parameters, detector) for detector in self.detectors]
        
        self.polarizations, self.timevector = hphc_amplitudes(
            'mlgw_bns', 
            parameters, 
            self.all_freqs, plot=None)

        tc = self.timevector[-1, 0]
        self.time_left = tc - self.timevector[:,0]
        # print(time_left)
        # index = next((i for i, t in enumerate(time_left) if t < 0), 0)
        # self.time_left = time_left[:index]

    @property        
    def all_freqs(self):
        return np.sort(
            np.concatenate([sig.frequencyvector for sig in self.signals])
            )[:, np.newaxis]
    
    @property
    def cumulative_snr(self):
        interpolants = [
            interp1d(
                sig.frequencyvector, 
                sig.cumulative_snr, 
                fill_value=(0, sig.cumulative_snr[-1]),
                bounds_error=False)
            for sig in self.signals
        ]
        
        interpolated_snrs = [
            interp(self.all_freqs) for interp in interpolants
        ]
        
        return np.sqrt(np.sum(np.array(interpolated_snrs)**2, axis=0))


def characteristic_strain_from_psd(frequencies, psd):
    return np.sqrt(frequencies * psd)

def characteristic_strain_from_fourier_amplitude(frequencies, amplitude):
    return 2 * frequencies * amplitude

def main():
    cmap = plt.get_cmap('inferno')
    lgwa_color = cmap(.8)
    et_color = cmap(.5)
    lisa_color = cmap(.2)
    signal_color = cmap(0.)
    
    # config = 'GWFish/detectors.yaml'
    
    et = detection.Detector('ET', parameters= [None], fisher_parameters= [None])
    lgwa = detection.Detector('LGWA', parameters= [None], fisher_parameters= [None])
    # lisa = detection.Detector('LISA', parameters= [None], fisher_parameters= [None], config=config)
    
    parameters = {
        'phase': 0.,
        'geocent_time': 1577491218., # 1st of January 2030
        'redshift': .01,
        'luminosity_distance': 40,
        'theta_jn': 1.,
        'mass_1': 1.4,
        'mass_2': 1.4,
        'ra': 0.,
        'dec': 0.,
        'psi': 0.,
        'lambda_1': 500., 
        'lambda_2': 500., 
        'a_1': 500., 
        'a_2': 500., 
    }
    
    signal_across_detectors = SignalAcrossDetectors(parameters, detectors=[
        # lisa, 
        lgwa, 
        et
    ])
    colors = [
        # lisa_color, 
        lgwa_color, 
        et_color,
    ]
    
    fix, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 6))

    for signal, color in zip(signal_across_detectors.signals, colors):
        
        axs[0].loglog(
            signal.frequencyvector,
            signal.characteristic_strain_signal,
            color=color, 
            ls='-', 
            lw =.3,
            alpha=.5)

        axs[0].loglog(
            signal.frequencyvector, 
            signal.characteristic_strain_noise, 
            label=signal.detector.name, 
            color=color)
    
        axs[1].loglog(
            signal.frequencyvector, 
            signal.snr_per_efolding * np.log(10),
            label=f'{signal.detector.name} SNR per decade',
            color=color, 
            ls=':'
        )

        axs[1].loglog(
            signal.frequencyvector, 
            signal.cumulative_snr,
            label=f'{signal.detector.name} cumulative SNR',
            color=color, 
        )

    times = {
        '1 second': 1,
        '1 minute': 60,
        '1 hour': 60*60,
        '1 day': 60*60*24,
        '1 month': 60*60*24*30,
        '1 year': 60*60*24*365,
    }

    axs[0].set_ylabel('Characteristic strain')
    axs[0].legend()
    
    axs[1].loglog(
        signal_across_detectors.all_freqs, 
        signal_across_detectors.cumulative_snr, 
        color=signal_color,
        label='Overall cumulative SNR'
    )

    axs[1].loglog(
        signal_across_detectors.all_freqs,
        10 * np.ones_like(signal_across_detectors.all_freqs),
        color='black',
        lw=.5,
        label='SNR threshold: 10'
    )
    
    signal_polarization = characteristic_strain_from_fourier_amplitude(
            signal_across_detectors.all_freqs[:, 0],
            abs(signal_across_detectors.polarizations[:, 0])
        )
    axs[0].loglog(
        signal_across_detectors.all_freqs[:, 0], 
        signal_polarization,
        label='GW170817',
        color=signal_color, 
        ls='--')


    for label, time in times.items(): 
        i = np.searchsorted(signal_across_detectors.time_left[::-1], time)
        i = len(signal_across_detectors.time_left) - i
        freq = signal_across_detectors.all_freqs[i, 0]
        strain = signal_polarization[i]
        
        axs[0].scatter([freq], [strain], color=signal_color, marker='+')
        axs[0].text(freq, strain*2, label, rotation=45.)

        # snr = np.sqrt(
        #     interp_et_cum_snr(freq)**2+
        #     interp_lgwa_cum_snr(freq)**2
        # )

        # axs[1].scatter([freq], [snr], color=signal_color, marker='x')

    axs[1].set_xlabel('Frequency [Hz]')
    axs[1].set_ylabel('SNR per decade / cumulative SNR')
    axs[0].set_xlim(1e-1, 2e3)
    axs[0].set_ylim(1e-24, 1e-18)
    axs[1].set_ylim(1, 1e4)
    axs[1].legend()

    axs[0].set_title('A best-case-scenario BNS detection')
    plt.tight_layout()
    plt.savefig(FIGS / 'ET_LGWA_BNS.pdf')

if __name__ == '__main__':
    main()