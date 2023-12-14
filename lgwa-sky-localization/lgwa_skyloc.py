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
from gw_landscape.plot import (make_frequency_axis,
                               make_time_axis_fancy, set_color_cycle)
from gw_landscape.plot_signals import (chirp_mass,
                                       make_time_to_merger_axis)

from colors import get_color_dict
from parameters import BNS_PARAMS

NINETY_PERCENT_INTERVAL = (180 / np.pi) ** 2 * 4.605170185988092

FISHER_PARAMETERS = [
    "mass_1",
    "mass_2",
    "luminosity_distance",
    "theta_jn",
    "dec",
    "ra",
    "psi",
    "phase",
    "geocent_time",
    "lambda_1",
    "lambda_2",
]

FIG_DIR = Path(__file__).parent.parent / "plots"
PSD_PATH = Path().parent / 'data'

def compute_sky_localization(cutoffs):
    population = "BNS"

    # detectors = ["ET", "LGWA"]
    detectors = ['LGWA']
    # detectors = ['LBI-SUS']

    # networks = "[[0], [1], [0, 1]]"
    networks = '[[0]]'

    detectors_ids = detectors
    networks_ids = json.loads(networks)
    duty_cycle = False
    

    waveform_model = "IMRPhenomD_NRTidalv2"

    params = BNS_PARAMS

    parameters = pd.DataFrame({k: v * np.ones_like(cutoffs) for k, v in params.items()})
    parameters["max_frequency_cutoff"] = cutoffs

    threshold_SNR = np.array([0.0, 8.0])
    fisher_parameters = FISHER_PARAMETERS
    network = gw.detection.Network(
        detectors_ids,
        detection_SNR=threshold_SNR,
        parameters=parameters,
        fisher_parameters=fisher_parameters,
    )
    waveform_class = gw.waveforms.LALFD_Waveform

    networkSNR_sq = 0
    for k in np.arange(len(parameters)):
        parameter_values = parameters.iloc[k]

        networkSNR_sq = 0
        for d in np.arange(len(network.detectors)):
            data_params = {
                "frequencyvector": network.detectors[d].frequencyvector,
                "f_ref": 50.0,
            }
            waveform_obj = waveform_class(waveform_model, parameter_values, data_params)
            wave = waveform_obj()
            t_of_f = waveform_obj.t_of_f

            signal = gw.detection.projection(
                parameter_values, network.detectors[d], wave, t_of_f
            )
            if any(np.isnan(signal[:, 0])):
                breakpoint()

            SNRs = gw.detection.SNR(network.detectors[d], signal, duty_cycle=duty_cycle)
            print(f"{network.detectors[d].name}: SNR={np.sqrt(np.sum(SNRs ** 2)):.1f}")
            networkSNR_sq += np.sum(SNRs**2)
            network.detectors[d].SNR[k] = np.sqrt(np.sum(SNRs**2))

            network.detectors[d].fisher_matrix[k, :, :] = gw.fishermatrix.FisherMatrix(
                waveform_model,
                parameter_values,
                fisher_parameters,
                network.detectors[d],
                waveform_class=waveform_class,
            ).fm

        network.SNR[k] = np.sqrt(networkSNR_sq)

    gw.detection.analyzeDetections(network, parameters, population, networks_ids)

    return gw.fishermatrix.analyzeFisherErrors(
        network, parameters, fisher_parameters, population, networks_ids
    )

def vary_mass(ns):
    
    rng = np.random.default_rng()
    
    mass = np.geomspace(1., 100., num=ns)
    distance = 30 * mass
    
    parameters = pd.DataFrame.from_dict(
        BNS_PARAMS | {
            "mass_1": mass,
            "mass_2": mass,
            "luminosity_distance": distance,
            # "redshift": 0.00980 * np.ones((ns,)),
            # "theta_jn": np.arccos(rng.uniform(-1.0, 1.0, size=(ns,))),
            # "dec": np.arccos(rng.uniform(-1.0, 1.0, size=(ns,))) - np.pi / 2.0,
            # "ra": rng.uniform(0, 2.0 * np.pi, size=(ns,)),
            # "psi": rng.uniform(0, 2.0 * np.pi, size=(ns,)),
            # "phase": rng.uniform(0, 2.0 * np.pi, size=(ns,)),
            # "geocent_time": rng.uniform(1735257618, 1766793618, size=(ns,)),
            # "a_1": 0.005136138323169717 * np.ones((ns,)),
            # "a_2": 0.003235146993487445 * np.ones((ns,)),
            # "lambda_1": 368.17802383555687 * np.ones((ns,)),
            # "lambda_2": 586.5487031450857 * np.ones((ns,)),
        }
    )
    
    population = 'mass_varied'
    return compute_fisher_multiple_params(parameters, population)

def vary_distance():
        
    distance = np.geomspace(40, 120, num=10)
    
    parameters = pd.DataFrame.from_dict(
        BNS_PARAMS | {
            "luminosity_distance": distance,
        }
    )
    
    population = 'varied_distance'
    return compute_fisher_multiple_params(parameters, population)


def randomize_sky_location(ns):
    rng = np.random.default_rng()

    parameters = pd.DataFrame.from_dict(
        {
            "mass_1": 1.4957673 * np.ones((ns,)),
            "mass_2": 1.24276395 * np.ones((ns,)),
            "redshift": 0.00980 * np.ones((ns,)),
            "luminosity_distance": 43.74755446 * np.ones((ns,)),
            "theta_jn": np.arccos(rng.uniform(-1.0, 1.0, size=(ns,))),
            "dec": np.arccos(rng.uniform(-1.0, 1.0, size=(ns,))) - np.pi / 2.0,
            "ra": rng.uniform(0, 2.0 * np.pi, size=(ns,)),
            "psi": rng.uniform(0, 2.0 * np.pi, size=(ns,)),
            "phase": rng.uniform(0, 2.0 * np.pi, size=(ns,)),
            "geocent_time": rng.uniform(1735257618, 1766793618, size=(ns,)),
            "a_1": 0.005136138323169717 * np.ones((ns,)),
            "a_2": 0.003235146993487445 * np.ones((ns,)),
            "lambda_1": 368.17802383555687 * np.ones((ns,)),
            "lambda_2": 586.5487031450857 * np.ones((ns,)),
        }
    )

    population = "BNS_randomized"
    return compute_fisher_multiple_params(parameters, population)

def compute_fisher_multiple_params(parameters, population):

    detectors = ["ET", "LGWA"]
    # detectors = ['LGWA']

    networks = "[[0], [1], [0, 1]]"
    # networks = '[[0, 1]]'

    detectors_ids = detectors
    networks_ids = json.loads(networks)
    duty_cycle = False

    waveform_model = "IMRPhenomD_NRTidalv2"

    threshold_SNR = np.array([0.0, 8.0])
    fisher_parameters = FISHER_PARAMETERS
    network = gw.detection.Network(
        detectors_ids,
        detection_SNR=threshold_SNR,
        parameters=parameters,
        fisher_parameters=fisher_parameters,
    )
    waveform_class = gw.waveforms.LALFD_Waveform

    networkSNR_sq = 0
    for k in np.arange(len(parameters)):
        parameter_values = parameters.iloc[k]

        networkSNR_sq = 0
        for d in np.arange(len(network.detectors)):
            data_params = {
                "frequencyvector": network.detectors[d].frequencyvector,
                "f_ref": 50.0,
            }
            waveform_obj = waveform_class(waveform_model, parameter_values, data_params)
            wave = waveform_obj()
            t_of_f = waveform_obj.t_of_f

            signal = gw.detection.projection(
                parameter_values, network.detectors[d], wave, t_of_f
            )

            SNRs = gw.detection.SNR(network.detectors[d], signal, duty_cycle=duty_cycle)
            print(f"{network.detectors[d].name}: SNR={np.sqrt(np.sum(SNRs ** 2)):.1f}")
            networkSNR_sq += np.sum(SNRs**2)
            network.detectors[d].SNR[k] = np.sqrt(np.sum(SNRs**2))

            network.detectors[d].fisher_matrix[k, :, :] = gw.fishermatrix.FisherMatrix(
                waveform_model,
                parameter_values,
                fisher_parameters,
                network.detectors[d],
                waveform_class=waveform_class,
            ).fm

        network.SNR[k] = np.sqrt(networkSNR_sq)

    gw.detection.analyzeDetections(network, parameters, population, networks_ids)

    return gw.fishermatrix.analyzeFisherErrors(
        network, parameters, fisher_parameters, population, networks_ids
    )


def skyloc_in_time_plot(fig_path, recompute=True):
    # cutoffs = [1.]
    cutoffs = np.geomspace(0.1, 30.0, num=40)

    if recompute:
        compute_sky_localization(cutoffs)
    # err_combined, = compute_sky_localization(cutoffs)
    err_et = pd.read_csv("Errors_ET_BNS_SNR8.0.txt", sep=" ")["err_sky_location"]
    err_lgwa = pd.read_csv("Errors_LGWA_BNS_SNR8.0.txt", sep=" ")["err_sky_location"]
    err_combined = pd.read_csv("Errors_ET_LGWA_BNS_SNR8.0.txt", sep=" ")[
        "err_sky_location"
    ]

    plt.rcParams.update({"text.usetex": True, "font.family": "Serif"})
    warnings.filterwarnings("ignore")

    colors = get_color_dict()
    
    plt.loglog(cutoffs[-len(err_et) :], err_et * NINETY_PERCENT_INTERVAL, c=colors["ET-LF"], label='ET')
    plt.loglog(cutoffs[-len(err_lgwa) :], err_lgwa * NINETY_PERCENT_INTERVAL, c=colors["LGWA"], label='LGWA')
    plt.loglog(cutoffs[-len(err_combined) :], err_combined * NINETY_PERCENT_INTERVAL, c=colors["all"], label='ET + LGWA')
    plt.legend()
    make_frequency_axis()
    # make_time_axis_fancy()
    subdivisions = {
        1.: 'second',
        60.: 'minute',
        3600.: 'hour',
        3600*24.: 'day',
        3600*24.*7.: 'week',
        3600*24.*30.: 'month',
        3600*24*365.24: 'year',
    }
    
    make_time_to_merger_axis(chirp_mass(BNS_PARAMS["mass_1"], BNS_PARAMS["mass_2"]), subdivisions=subdivisions)
    plt.ylabel('90\\% sky area [square degrees]')

    plt.savefig(fig_path)


def randomized_skyloc_plot(fig_path):
    # err_et, err_lgwa, err_combined = randomize_sky_location(1000)
    err_et = pd.read_csv("Errors_ET_BNS_randomized_SNR8.0.txt", sep=" ")["err_sky_location"]
    err_lgwa = pd.read_csv("Errors_LGWA_BNS_randomized_SNR8.0.txt", sep=" ")["err_sky_location"]
    err_combined = pd.read_csv("Errors_ET_LGWA_BNS_randomized_SNR8.0.txt", sep=" ")[
        "err_sky_location"
    ]
    plt.hist(np.log(err_lgwa * NINETY_PERCENT_INTERVAL), alpha=0.5, bins=30, label='LGWA')
    plt.hist(np.log(err_et * NINETY_PERCENT_INTERVAL), alpha=0.5, bins=30, label='ET')
    plt.hist(np.log(err_combined * NINETY_PERCENT_INTERVAL), alpha=0.5, bins=30, label='ET+LGWA')
    plt.gca().xaxis.set_major_formatter(lambda x, pos: f"${np.exp(x)*3600:.1f}$")

    plt.xlabel('90% sky area [square arcminutes]')
    plt.ylabel('Number of events')
    plt.legend()
    plt.savefig(fig_path)

def varied_mass_plot(fig_path):
    ns = 20
    # vary_mass(ns)
    err_et = pd.to_numeric(pd.read_csv("Errors_ET_mass_varied_SNR8.0.txt", sep=" ").replace('NAN', np.nan)["err_sky_location"])
    mass_et = pd.read_csv("Errors_ET_mass_varied_SNR8.0.txt", sep=" ").replace('NAN', np.nan)["mass_1"]
    err_lgwa = pd.read_csv("Errors_LGWA_mass_varied_SNR8.0.txt", sep=" ")["err_sky_location"]
    mass_lgwa = pd.read_csv("Errors_LGWA_mass_varied_SNR8.0.txt", sep=" ")["mass_1"]
    # err_combined = np.array(pd.read_csv("Errors_ET_LGWA_mass_varied_SNR8.0.txt", sep=" ")[
    #     "err_sky_location"
    # ])
    # mass_combined = np.array(pd.read_csv("Errors_ET_LGWA_mass_varied_SNR8.0.txt", sep=" ")["mass_1"])
    
    # breakpoint()
    
    plt.loglog(mass_et, err_et*NINETY_PERCENT_INTERVAL, label='ET')
    plt.loglog(mass_lgwa, err_lgwa*NINETY_PERCENT_INTERVAL, label='LGWA')
    # plt.loglog(mass_combined, err_combined*NINETY_PERCENT_INTERVAL, label='ET+LGWA')

    plt.xlabel('Mass [$M_\odot$]')
    plt.ylabel('90% sky area [square degrees]')
    plt.legend()
    plt.savefig(fig_path)

def varied_distance_plot(fig_path):
    vary_distance()
    err_et = pd.to_numeric(pd.read_csv("Errors_ET_varied_distance_SNR8.0.txt", sep=" ").replace('NAN', np.nan)["err_sky_location"])
    dist = pd.read_csv("Errors_ET_varied_distance_SNR8.0.txt", sep=" ").replace('NAN', np.nan)["luminosity_distance"]
    err_lgwa = pd.read_csv("Errors_LGWA_varied_distance_SNR8.0.txt", sep=" ")["err_sky_location"]
    err_combined = np.array(pd.read_csv("Errors_ET_LGWA_varied_distance_SNR8.0.txt", sep=" ")[
        "err_sky_location"
    ])
    
    # breakpoint()
    
    plt.loglog(dist, err_et*NINETY_PERCENT_INTERVAL, label='ET')
    plt.loglog(dist, err_lgwa*NINETY_PERCENT_INTERVAL, label='LGWA')
    plt.loglog(dist, err_combined*NINETY_PERCENT_INTERVAL, label='ET+LGWA')

    plt.xlabel('Distance [Mpc]')
    plt.ylabel('90% sky area [square degrees]')
    plt.legend()
    plt.savefig(fig_path)
    
if __name__ == "__main__":
    # skyloc_in_time_plot(FIG_DIR / "localizations.pdf", recompute=True)
    # randomized_skyloc_plot(FIG_DIR / "localizations_histogram.png")
    # varied_distance_plot(FIG_DIR / "varied_distance.png")
    # varied_mass_plot(FIG_DIR / "varied_mass.png")
    print(compute_sky_localization([10.]))