import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "text.latex.preamble": r"\usepackage{mathtools}",
})
import periodictable as pt
import argparse


def parse_ensdf_levels_gammas(path):
    # Dictionary of the nuclear levels
    levels = {}
    parent = {}

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        # Read all lines and store them
        lines = f.readlines()
    # print(lines[0])

    Ng = 0
    Ndec = 0

    # Iterate over lines in the file
    for i, line in enumerate(lines):
        line = line.rstrip("\n")
        # print(line)

        # This is the record identifier (ex. L=level, G=gamma)
        rec = line[5:8].strip()
        # print(rec, line)

        # First line for parent and daughter
        if i == 0:
            # Daughter properties
            daughter_A = line[:3]
            daughter_sy = line[3].upper() + line[4].lower()
            daughter_Z = getattr(pt, daughter_sy).number
            # Parent properties
            parent_A = line[9:12]
            parent_sy = line[12].upper() + line[13].lower()
            parent_Z = getattr(pt, parent_sy).number

            parent.update({"d_A":daughter_A, "d_sy":daughter_sy, "d_Z":daughter_Z,
                           "p_A":parent_A, "p_sy":parent_sy, "p_Z":parent_Z})

        # Look for the parent first
        if rec == "P":
            # Energy in keV
            E_p = float(line[9:19].strip())
            # Spin-parity
            JP_p = line[21:39].strip()
            # Half life with unit
            T12_p = line[39:49].strip().lower()
            # Q-value in keV
            Q_p = float(line[64:74].strip())

            parent.update({"E":E_p, "JP":JP_p, "T12":T12_p, "Q":Q_p})
        
        # Look for the normalization
        if rec == "N":
            # The normalization converts relative gamma intensity to gamma intensity
            # So the result is gammas per 100 decays (if in percent)
            Ng = float(line[9:19].strip())
            # Also for beta decay
            Ndec = float(line[41:49].strip())

        # Look for a level record
        if rec == "L":
            # Energy in keV
            E = float(line[9:19].strip())
            # Spin-parity
            JP = line[21:39].strip()
            # Half life with unit
            T12 = line[39:49].strip().lower()
            # print(E, JP, T12)

            levels[E] = {"E":E, "JP":JP, "T12":T12, "G":[]}
        
        # Look for the decay record to the level record
        if rec == "B" or rec == "EC":
            # Intensity of the decay to reach this state
            Idec = float(line[21:29].strip()) * Ndec
            if Idec == "":
                Idec = 0
            # Type of decay
            dec = rec

            levels[E].update({"Idec":Idec, "dec":dec})
        
        # What about alpha decay, rec == "A"?

        # Look for gamma rays from this level
        if rec == "G":
            # print(line)
            # Energy in keV
            Eg = float(line[9:19].strip())
            # Relative photon intensity
            Ig = line[21:29].strip()
            if Ig != "":
                Ig = float(Ig) * Ng
            # Multipolarity
            M = line[31:41].strip()
            # Mixing ratio
            MR = line[41:49].strip()

            levels[E]["G"].append({"Eg":Eg, "Ig":Ig, "M":M, "MR":MR})

    # print(levels)
    # print(all_levels)
    return levels, parent


def limit_intensity(levels, I_min=0, I_max=100):
    # To store filtered levels and gammas
    new_levels = {}

    # Go though all data in the previous dictionary
    for key, value in levels.items():
        # Gamma rays that have intensity >= than the specified I_min
        gammas = [dic for dic in value["G"] if dic["Ig"] != "" and dic["Ig"] >= I_min and dic["Ig"] <= I_max]
        # print(gammas)
        # print(key, value)

        # Only use the intense gamma rays
        value["G"] = gammas

        # Only store the intense gamma ray levels and data
        if gammas != [] or key == 0.0:
            new_levels[key] = value

    # Make a list of all levels in order
    new_all_levels = sorted(new_levels.keys(), reverse=True)

    # print(new_levels)
    return new_levels, new_all_levels


def plot_decay_scheme(data_path, I_min=0, I_max=100, save_path="figure"):
    # Open, parse, and save data file contents
    levels, parent = parse_ensdf_levels_gammas(data_path)
    levels, all_levels = limit_intensity(levels, I_min, I_max)

    plt.close('all')
    inch_to_mm = 25.4
    color = plt.cm.tab10

    fig, ax = plt.subplots(figsize=(88/inch_to_mm,80/inch_to_mm))

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_position([0, 0, 1, 1], which="both")
    ax.set_frame_on(False)

    energy_span = max(levels.keys())
    spacing_arrow_background = energy_span / 20

    y_pad_under = 0.15
    y_pad_over = 0.6
    ax.set_ylim([-y_pad_under*energy_span, (1+y_pad_over)*energy_span])

    number_of_gammas = sum([len(level["G"]) for level in levels.values()])
    # print(number_of_gammas)

    x_min = 0.0
    x_max = 1.0
    x_gamma_min = 0.15
    x_gamma_max = 0.55
    x_step = (x_gamma_max - x_gamma_min) / (number_of_gammas - 1)
    x_pos = x_gamma_min

    x_pad_left = 0.3
    x_pad_right = 0.25
    ax.set_xlim([-x_pad_left*x_max, (1+x_pad_right)*x_max])

    # for key, value in levels.items():
    for E in all_levels:
        # E = key
        level = levels[E]
        JP = r"${}^{}$".format(level["JP"][:-1], level["JP"][-1])
        T12 = level["T12"]
        try:
            Idec = level["Idec"]
        except KeyError:
            Idec = "?"

        if E == 0.0:
            ax.hlines(E, x_min, x_max, color="black", linewidth=1.5, zorder=4)
            ax.text(x_max-0.01, E, f"{E} keV", ha="right", va="bottom", fontsize=8, zorder=4)
        else:
            ax.hlines(E, x_min, x_max, color="black", linewidth=1.0, zorder=1)
            ax.text(x_max-0.01, E, f"{E}", ha="right", va="bottom", fontsize=8, zorder=4)
        
        ax.text(x_max+0.02, E, f"{T12}", ha="left", va="center", fontsize=8, zorder=4)

        for gamma in level["G"]:
            Eg = gamma["Eg"]
            Ig = gamma["Ig"]

            ax.plot([x_pos, x_pos], [E-spacing_arrow_background, E-Eg+spacing_arrow_background], 
                    color="white", linewidth=4, zorder=2)
            ax.annotate(None, xy=(x_pos, E-Eg), xytext=(x_pos, E), 
                        arrowprops=dict(arrowstyle="->", lw=1.0, shrinkA=0, shrinkB=0),
                        zorder=3)
            ax.text(x_pos+0.01, E+0.005*energy_span, f"{Eg:.2f} keV {Ig:.2f}\%",
                    rotation=60, rotation_mode="anchor", ha="left", va="bottom", fontsize=8,
                    bbox=dict(boxstyle="square,pad=0.01", fc="white", ec="none"),
                    zorder=4)
            
            x_pos += x_step

        x_pos_decay_arrow_left = -x_pad_left * x_max + 0.04
        x_pos_decay_arrow_right = -0.01
        ax.annotate(None, xy=(x_pos_decay_arrow_right, E), xytext=(x_pos_decay_arrow_left, E), 
                    arrowprops=dict(arrowstyle="->", lw=1.0, shrinkA=0, shrinkB=0),
                    zorder=4)
        ax.text(x_pos_decay_arrow_right-0.04, E, f"{Idec}\%", ha="right", va="bottom", fontsize=8, zorder=4)
        ax.text(x_min+0.01, E, f"${JP}$", ha="left", va="bottom", fontsize=8, zorder=4)
    
    y_pad_over_space = 0.15
    y_pos_arrow = (1 + y_pad_over - y_pad_over_space) * energy_span
    ax.annotate(None, xy=(x_pos_decay_arrow_left, 0.0), xytext=(x_pos_decay_arrow_left, y_pos_arrow), 
                arrowprops=dict(arrowstyle="-", lw=1.0, shrinkA=0, shrinkB=0),
                zorder=4)
    ax.hlines(y_pos_arrow, x_pos_decay_arrow_left-0.003, x_pos_decay_arrow_right, color="black", linewidth=1.5, zorder=4)

    E_p = parent["E"]
    JP_p = r"${}^{}$".format(parent["JP"][:-1], parent["JP"][-1])
    T12_p = parent["T12"]
    Q_p = parent["Q"]
    parent_str = rf"$\mathrm{{\prescript{{{parent['p_A']}}}{{{parent['p_Z']}}}{{{parent['p_sy']}}}}}$"
    daughter_str = rf"$\mathrm{{\prescript{{{parent['d_A']}}}{{{parent['d_Z']}}}{{{parent['d_sy']}}}}}$"
    ax.text(x_pos_decay_arrow_left+0.01, y_pos_arrow, f"{JP_p}", ha="left", va="bottom", fontsize=8, zorder=4)
    ax.text(x_pos_decay_arrow_right-0.01, y_pos_arrow, f"{E_p}", ha="right", va="bottom", fontsize=8, zorder=4)
    ax.text((x_pos_decay_arrow_left+x_pos_decay_arrow_right)/2, y_pos_arrow+0.06*energy_span, f"{T12_p}", ha="center", va="bottom", fontsize=8, zorder=4)
    ax.text(x_pos_decay_arrow_left, y_pos_arrow-(y_pad_over-y_pad_over_space)/2*energy_span, f"//",
            rotation=90, rotation_mode="anchor", ha="center", va="center", fontsize=6,
            bbox=dict(boxstyle="square,pad=0.01", fc="white", ec="none"),
            zorder=4)
    # q_text=r"$Q_{\beta\text{-}}$"
    # ax.text(x_pos_decay_arrow_left+0.03, y_pos_arrow-(y_pad_over-y_pad_over_space)/2*energy_span, f"{q_text}={Q_p} keV", ha="left", va="center", fontsize=8, zorder=4)
    # ax.text((x_pos_decay_arrow_left+x_pos_decay_arrow_right)/2, y_pos_arrow-0.2*energy_span, f"{q_text}{Q_p} keV", ha="center", va="bottom", fontsize=8, zorder=4)

    ax.text(0.5, 0.0-0.03*energy_span, daughter_str, ha="center", va="top", fontsize=10, zorder=4)
    ax.text((x_pos_decay_arrow_left+x_pos_decay_arrow_right)/2, y_pos_arrow-0.03*energy_span, parent_str, ha="center", va="top", fontsize=10, zorder=4)

    ax.text((1+x_pad_right)*x_max-0.02, y_pos_arrow+0.06*energy_span, 
            rf"(shown: ${{{I_min:.0f}}}\% \leq I_{{\gamma}} \leq {{{I_max:.0f}}}\%$)", 
            ha="right", va="bottom", fontsize=8, zorder=4)

    plt.savefig(f'{save_path}.jpg', dpi=600)
    plt.savefig(f'{save_path}.pdf')

    plt.show()


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str, required=True, help="ENSDF file path")
parser.add_argument("-s", "--save", type=str, default="figure", help="save path (default: 'figure')")
parser.add_argument("-imin", "--imin", type=float, default=0.0, help="Minimum gamma intensity (default: 0.0)")
parser.add_argument("-imax", "--imax", type=float, default=100.0, help="Maximum gamma intensity (default: 100.0)")
args = parser.parse_args()

# print(levels)
# print(parent)
# data_path = "ensdf_files/140la.txt"
# plot_decay_scheme(data_path, I_min=10, I_max=100, save_path="figures/figure")

plot_decay_scheme(data_path=args.file, I_min=args.imin, I_max=args.imax, save_path=args.save)
