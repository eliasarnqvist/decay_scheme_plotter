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


def parse_ensdf_levels_gammas(path, nucid):
    # Get the nuclear id
    nucid = nucid[:5].upper()
    # Dictionary of the nuclear levels
    levels = {}
    parent = {}

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        # Read all lines and store them
        lines = f.readlines()
    # print(lines[0])

    # Iterate over lines in the file
    for i, line in enumerate(lines):
        line = line.rstrip("\n")
        # print(line)

        # This is the record identifier (ex. L=level, G=gamma)
        rec = line[5:8].strip()
        # print(rec, line)

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

            parent = {"E":E_p, "JP":JP_p, "T12":T12_p, "Q":Q_p}

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
        if rec == "B" or rec == "EC" or rec == "A":
            # Intensity of the decay to reach this state
            Idec = float(line[21:29].strip())
            # Type of decay
            dec = rec

            levels[E].update({"Idec":Idec, "dec":dec})
        
        # Look for gamma rays from this level
        if rec == "G":
            # print(line)
            # Energy in keV
            Eg = float(line[9:19].strip())
            # Relative photon intensity
            Ig = line[21:29].strip()
            if Ig != "":
                Ig = float(Ig)
            # Multipolarity
            M = line[31:41].strip()
            # Mixing ratio
            MR = line[41:49].strip()

            levels[E]["G"].append({"Eg":Eg, "Ig":Ig, "M":M, "MR":MR})

    # print(levels)
    # print(all_levels)
    return levels, parent


def limit_intensity(levels, I_min=0):
    # To store filtered levels and gammas
    new_levels = {}

    # Go though all data in the previous dictionary
    for key, value in levels.items():
        # Gamma rays that have intensity >= than the specified I_min
        gammas = [dic for dic in value["G"] if dic["Ig"] != "" and dic["Ig"] >= I_min]
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


def plot_decay_scheme(levels, all_levels, parent, parent_name="??", daughter_name="??", q_text=""):
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
    y_pad_over = 0.5
    ax.set_ylim([-y_pad_under*energy_span, (1+y_pad_over)*energy_span])

    number_of_gammas = sum([len(level["G"]) for level in levels.values()])
    # print(number_of_gammas)

    x_min = 0.0
    x_max = 1.0
    x_gamma_min = 0.2
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
        JP = level["JP"]
        T12 = level["T12"]
        Idec = level["Idec"]

        if E == 0.0:
            ax.hlines(E, x_min, x_max, color="black", linewidth=2.0, zorder=4)
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
            ax.text(x_pos, E+0.005*energy_span, f"{Eg:.2f} keV {Ig}\%",
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
    ax.hlines(y_pos_arrow, x_pos_decay_arrow_left-0.003, x_pos_decay_arrow_right, color="black", linewidth=2.0, zorder=4)

    ax.text(0.5, 0.0-0.03*energy_span, daughter_name, ha="center", va="top", fontsize=10, zorder=4)
    ax.text((x_pos_decay_arrow_left+x_pos_decay_arrow_right)/2, y_pos_arrow-0.03*energy_span, parent_name, ha="center", va="top", fontsize=10, zorder=4)

    E_p = parent["E"]
    JP_p = parent["JP"]
    T12_p = parent["T12"]
    Q_p = parent["Q"]
    ax.text(x_pos_decay_arrow_left+0.01, y_pos_arrow, f"{JP_p}", ha="left", va="bottom", fontsize=8, zorder=4)
    ax.text(x_pos_decay_arrow_right-0.01, y_pos_arrow, f"{E_p}", ha="right", va="bottom", fontsize=8, zorder=4)
    ax.text((x_pos_decay_arrow_left+x_pos_decay_arrow_right)/2, y_pos_arrow+0.06*energy_span, f"{T12_p}", ha="center", va="bottom", fontsize=8, zorder=4)
    ax.text(x_pos_decay_arrow_left, y_pos_arrow-(y_pad_over-y_pad_over_space)/2*energy_span, f"//",
            rotation=90, rotation_mode="anchor", ha="center", va="center", fontsize=6,
            bbox=dict(boxstyle="square,pad=0.01", fc="white", ec="none"),
            zorder=4)
    ax.text(x_pos_decay_arrow_left+0.03, y_pos_arrow-(y_pad_over-y_pad_over_space)/2*energy_span, f"{q_text}={Q_p} keV", ha="left", va="center", fontsize=8, zorder=4)
    # ax.text((x_pos_decay_arrow_left+x_pos_decay_arrow_right)/2, y_pos_arrow-0.2*energy_span, f"{q_text}{Q_p} keV", ha="center", va="bottom", fontsize=8, zorder=4)

    save_name = "figure"
    # plt.savefig(f'figures/{save_name}.jpg', dpi=600)
    plt.savefig(f'figures/{save_name}.pdf')

    plt.show()


levels, parent = parse_ensdf_levels_gammas("ensdf_files/140la.txt", nucid="140CE")
levels, all_levels = limit_intensity(levels, I_min=10)

print(levels)
print(parent)

plot_decay_scheme(levels, all_levels, parent,
                  parent_name=r"$\mathrm{\prescript{140}{57}{La}}$", 
                  daughter_name=r"$\mathrm{\prescript{140}{56}{Ce}}$",
                  q_text=r"$Q_{\beta\text{-}}$")
