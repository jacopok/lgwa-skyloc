import json
import os
import warnings
from pathlib import Path

import astropy.constants as ac
import astropy.units as u
import GWFish.modules as gw
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from gw_landscape.einstein_telescope import (EinsteinTelescopeCryo,
                                             EinsteinTelescopeHF)
from gw_landscape.gwfish import GWFishDetector
from gw_landscape.lisa import LISA
from gw_landscape.plot import make_frequency_axis, make_time_axis_fancy
from gw_landscape.plot_signals import (chirp_mass, get_projection,
                                       make_time_to_merger_axis,
                                       plot_characteristic_noise_strain,
                                       plot_characteristic_signal_strain, set_color_cycle)
from scipy.interpolate import interp1d

from lgwa_skyloc import FIG_DIR
from colors import get_color_dict
from parameters import BNS_PARAMS


def plot_bns_lgwa_et(fig_path):
    warnings.filterwarnings("ignore")

    detector_list = [
        GWFishDetector("LGWA"),
        EinsteinTelescopeCryo(),
        EinsteinTelescopeHF(),
        LISA(),
    ]
    color_dict = get_color_dict()
    color_list = [
        color_dict['LGWA'],
        color_dict['ET-LF'],
        color_dict['ET-HF'],
        color_dict['LISA'],
    ]

    detector_list[3].annotation_place = (2e-3, 5e-21)
    detector_list[0].annotation_place = (2e-3, 4e-17)
    detector_list[1].annotation_place = (2e3, 4e-19)
    detector_list[2].annotation_place = (2e3, 1.5e-21)

    plot_characteristic_noise_strain(detector_list, colors=color_list)

    approx = "IMRPhenomD_NRTidalv2"

    plot_characteristic_signal_strain(
        BNS_PARAMS,
        GWFishDetector("LGWA").gdet,
        waveform_model=approx,
        color=color_dict["LGWA"],
        ls="-",
        lw=1,
        alpha=0.8,
        label='LGWA projection'
    )
    plot_characteristic_signal_strain(
        BNS_PARAMS,
        GWFishDetector("ET").gdet,
        waveform_model=approx,
        color=color_dict["ET-LF"],
        ls="-",
        lw=1,
        alpha=0.8,
        label='ET projection'
    )

    plt.legend()
    ax_freq = make_frequency_axis()

    make_time_to_merger_axis(chirp_mass(BNS_PARAMS["mass_1"], BNS_PARAMS["mass_2"]))

    plt.ylabel('Characteristic signal/noise strain')
    plt.ylim(1e-24, 1e-16)
    plt.xlim(1e-3, 1e4)

    plt.savefig(fig_path, dpi=200)
    plt.close()


def plot_bns_lgwa_et_snr_integrand(fig_path):
    plt.rcParams.update({"text.usetex": True, "font.family": "Serif"})
    set_color_cycle()
    warnings.filterwarnings("ignore")
    lgwa = GWFishDetector("LGWA")
    et_lf = GWFishDetector("ET_LF", psd_path=Path(__file__).parent / "data")
    et_hf = GWFishDetector("ET_HF", psd_path=Path(__file__).parent / "data")
    detector_list = [
        lgwa,
        et_lf,
        et_hf,
    ]
    color_dict = get_color_dict()
    color_list = [
        color_dict['LGWA'],
        color_dict['ET-LF'],
        color_dict['ET-HF'],
    ]

    # colors = plot_characteristic_noise_strain(detector_list)

    params = BNS_PARAMS
    approx = "IMRPhenomD_NRTidalv2"

    frequencyvectors = []
    integrands = []
    for detector, color in zip(detector_list, color_list):
        proj, _ = get_projection(params, detector.gdet, waveform_model=approx)
        integrand = detector.snr_integrand(
            proj, detector.gdet.frequencyvector[:, 0], log_base=2
        )
        integrands.append(integrand / np.log(2))
        frequencyvectors.append(detector.gdet.frequencyvector[:, 0])
        plt.loglog(detector.gdet.frequencyvector[:, 0], np.sqrt(integrand), color=color)

    # et_proj, _ = get_projection(params, GWFishDetector('ET').gdet, waveform_model=approx, frequencies=lgwa)

    plt.ylabel("Colored: SNR per octave")
    make_frequency_axis()

    ax1_snr = plt.gca()
    ax2_snr = ax1_snr.twinx()

    f, i = compute_total_snr(frequencyvectors, integrands)

    ax2_snr.loglog(f, np.sqrt(i), c="black")
    ax2_snr.set_ylabel("Black: total accumulated SNR")
    subdivisions = {
        1.: 'second',
        60.: 'minute',
        3600.: 'hour',
        3600*24.: 'day',
        3600*24.*30.: 'month',
        3600*24*365.24: 'year',
    }


    make_time_to_merger_axis(
        chirp_mass(BNS_PARAMS["mass_1"], BNS_PARAMS["mass_2"]), subdivisions=subdivisions
    )
    # plt.legend()
    for ax in [ax1_snr, ax2_snr]:
        ax.set_xlim(5e-2, 1e3)
        ax.set_ylim(1, 1e3)

    plt.savefig(fig_path, dpi=200)
    plt.close()


def compute_total_snr(frequencyvectors, integrands):
    all_freqs = np.sort(np.concatenate(frequencyvectors))

    sum_squares = np.zeros_like(all_freqs)
    for f, i in zip(frequencyvectors, integrands):
        sum_squares += interp1d(f, i, fill_value=0.0, bounds_error=False)(all_freqs)

    to_integrate = sum_squares * np.gradient(np.log(all_freqs))
    return all_freqs, np.cumsum(to_integrate)


if __name__ == "__main__":
    plot_bns_lgwa_et(FIG_DIR / f"sensitivity_curves.pdf")
    plot_bns_lgwa_et_snr_integrand(FIG_DIR / f"snr_integrand.pdf")
