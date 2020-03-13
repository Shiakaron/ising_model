import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import itertools

folder = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\data\\energy_data\\vs_temp"
folder2 = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\pngs\\"
texfolder = "C:\\Users\\savva\\OneDrive - University of Cambridge\\Part2\\Computational Project\\report\\texfigures\\"

def plot_1():
    """
    energy vs temperature with different L's
    """
    p_files = []
    dim = 2
    for file in sorted(os.listdir(folder)):
        if file.endswith(".txt") and not file.endswith(").txt"):
            p_files.append(os.path.join(folder,file))

    L_list = []
    fig, ax = plt.subplots(figsize=(12,8))
    ax.axvline(2.2692, label="$T_c$", linestyle="--",color="k")
    marker = itertools.cycle(('*', '+', '.', ',', 'o'))
    for p_file in p_files:
        L = (os.path.splitext(os.path.basename(p_file))[0]).split('_',3)[3]
        avgE = []
        errE = []
        T = []
        if (L not in L_list):
            L_list.append(L)
            with open(p_file) as csvfile:
                lines = csv.reader(csvfile, delimiter=' ')
                for row in lines:
                    T.append(float(row[0]))
                    avgE.append(float(row[1]))
                    errE.append(float(row[2]))
            ax.errorbar(T, avgE, errE, ls='',marker = next(marker), label="L = "+str(L))

    ax.set_title("<Energy> vs Temperature")
    ax.set_ylabel("e / J")
    ax.set_xlabel(r"T / $J/k_B$")
    ax.legend()

    fig.savefig(texfolder+"energy_vs_temp.pdf")
    fig.savefig(folder2+"energy_vs_temp.png")


def main():
    plot_1()
    plt.show()

if (__name__ == '__main__'):
    main()
