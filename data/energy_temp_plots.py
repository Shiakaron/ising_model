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
    p_files_dict = {}
    dim = 2
    L_list = []
    for file in sorted(os.listdir(folder)):
        if file.endswith(".txt") and not file.endswith(").txt"):
            L = int((os.path.splitext(os.path.basename(file))[0]).split('_',3)[3])
            if L not in L_list:
                p_files_dict[L] = os.path.join(folder,file)
                L_list.append(L)

    L_list.sort()
    fig, ax = plt.subplots()
    ax.axvline(2.2692, label="$T_c$", linestyle="--",color="k", alpha=0.6)
    marker = itertools.cycle(('*', '+', '.', ',', 'o'))
    for L in L_list:
        p_file = p_files_dict[L]
        avgE = []
        errE = []
        T = []
        with open(p_file) as csvfile:
            lines = csv.reader(csvfile, delimiter=' ')
            for row in lines:
                T.append(float(row[0]))
                avgE.append(float(row[1]))
                errE.append(float(row[2]))
        ax.errorbar(T, avgE, errE, ls='',marker = next(marker), label="L = "+str(L), alpha=0.8)

    # ax.set_title("<Energy> vs Temperature")
    ax.set_ylabel("Energy per link (J)")
    ax.set_xlabel(r"T ($J/k_B$)")
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.subplots_adjust(right=1,top=1)
    ax.legend()

    fig.savefig(texfolder+"energy_vs_temp.pdf")
    fig.savefig(folder2+"energy_vs_temp.png")


def main():
    plot_1()
    plt.show()

if (__name__ == '__main__'):
    main()
