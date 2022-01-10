"""
Bernardo Pacini
Date Created: Feb. 3rd, 2021 Wed
Date Modified: Apr. 11th, 2021 Sun
Description: This script reads a json file with wind tunnel test data for the
    MR8x45 rotor and verifies HELIX against the recorded data. This script
    outputs plots for relevant recorded quantities.
"""
import json
import numpy as np
import matplotlib.pyplot as plt

import HELIX_component_ADAPTED as helix_comp

import niceplots

niceplots.setRCParams()

def main():
    # =========================================================================
    # Parameters
    # =========================================================================
    rpm_max = 10000.0

    tol = 0.1

    V_inf = np.linsapce([3, 16, 20])
    V_inf_tag = ["V_05"]

    AOA = np.array([0.0])
    AOA_tag = ["AOA_00"]

    N_RPM = 9

    N_VAL = 10
    # =========================================================================
    # Setup Data Arrays
    # =========================================================================
    rpm = np.zeros((np.size(V_inf), np.size(AOA), N_VAL))
    thrust = np.zeros((np.size(V_inf), np.size(AOA), N_VAL))

    # =========================================================================
    # Import data into dictionary
    # =========================================================================
    with open("/Users/joaquinexalto/Documents/TU_Delft/code/verification/HELIX 2/Data-MR8x45.json") as f:
        wt_data = json.load(f)

    # =========================================================================
    # Setup HELIX Case
    # =========================================================================
    # Pass V_inf, AOA, and RPM to function, return thrust
    for iCase in range(0, np.size(V_inf)):
        for iAOA in range(0, np.size(AOA)):
            # Pick Minimum Thrust
            iRPM = 0

            rpm_min = wt_data[V_inf_tag[iCase]][AOA_tag[iAOA]]["rpm"][iRPM]
            rpm[iCase, iAOA, :] = np.linspace(rpm_min, rpm_max, N_VAL)

            for iRPM in range(0, N_VAL):
                v = np.array(
                    [
                        -V_inf[iCase] * np.sin(AOA[iAOA] * np.pi / 180.0),
                        V_inf[iCase] * np.cos(AOA[iAOA] * np.pi / 180.0),
                        0.0,
                    ]
                )
                thrust[iCase, iAOA, iRPM] = helix_comp.run_helix(rpm[iCase, iAOA, iRPM], v)
    # =========================================================================
    # Plot Results
    # =========================================================================
    colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]
    labels = [
        r"$90^\degree$",
        r"$80^\degree$",
        r"$70^\degree$",
        r"$60^\degree$",
        r"$50^\degree$",
        r"$40^\degree$",
        r"$30^\degree$",
        r"$20^\degree$",
        r"$10^\degree$",
        r"$0^\degree$",
    ]

    # V = 5 m/s
    _, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(25,7))
    # _, ax1 = plt.subplots(figsize=(7,7))
    for iAOA in range(0, np.size(AOA)):
        # Plot Wind Tunnel Data
        ax1.scatter(
            wt_data[V_inf_tag[0]][AOA_tag[iAOA]]["rpm"][:],
            wt_data[V_inf_tag[0]][AOA_tag[iAOA]]["thrust"][:],
        )

        # Plot HELIX Result
        ax1.plot(
            rpm[0, iAOA, :],
            thrust[0, iAOA, :],
            colors[iAOA],
            # label=labels[iAOA],
        )
    # ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    ax1.set_xlim([-100,10000])
    ax1.set_ylim([-1,8])
    ax1.set_xlabel("RPM")
    ax1.set_ylabel(r"Thrust $[N]$")
    niceplots.adjust_spines(ax1, outward=True)
    # plt.tight_layout()

    # V = 10 m/s
    # _, ax2 = plt.subplots(figsize=(7,7))
    for iAOA in range(0, np.size(AOA)):
        # Plot Wind Tunnel Data
        ax2.scatter(
            wt_data[V_inf_tag[1]][AOA_tag[iAOA]]["rpm"][:],
            wt_data[V_inf_tag[1]][AOA_tag[iAOA]]["thrust"][:],
        )
        # Plot HELIX Result
        ax2.plot(
            rpm[1, iAOA, :],
            thrust[1, iAOA, :],
            colors[iAOA],
            # label=labels[iAOA],
        )
    # ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    ax2.set_xlim([-100,10000])
    ax2.set_ylim([-1,8])
    ax2.set_xlabel("RPM")
    ax2.set_ylabel(r"Thrust $[N]$")
    niceplots.adjust_spines(ax2, outward=True)
    # plt.tight_layout()

    # V = 15 m/s
    # _, ax3 = plt.subplots(figsize=(7,7))
    for iAOA in range(0, np.size(AOA)):
        # Plot Wind Tunnel Data
        ax3.scatter(
            wt_data[V_inf_tag[2]][AOA_tag[iAOA]]["rpm"][:],
            wt_data[V_inf_tag[2]][AOA_tag[iAOA]]["thrust"][:],
        )
        # Plot HELIX Result
        ax3.plot(
            rpm[2, iAOA, :],
            thrust[2, iAOA, :],
            colors[iAOA],
            label=labels[iAOA],
        )
    ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    ax3.set_xlim([-100,10000])
    ax3.set_ylim([-1,8])
    ax3.set_xlabel("RPM")
    ax3.set_ylabel(r"Thrust $[N]$")
    niceplots.adjust_spines(ax3, outward=True)
    # plt.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()
