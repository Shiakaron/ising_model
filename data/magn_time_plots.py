import numpy as np
import matplotlib.pyplot as plt
import csv
import os

# please fix filepaths if you are running this on different devices
folder = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\data\\magn_vs_time"
folder2 = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\pngs\\"
texfolder = "C:\\Users\\savva\\OneDrive - University of Cambridge\\Part2\\Computational Project\\report\\texfigures\\"
def thermalisation_1():
    """
    Plot at different temperatures but same L
    """
    p_files = []
    dim = 2
    L = 48
    filename = "magn_vs_time_"+str(dim)+"D_"+str(L)
    for file in sorted(os.listdir(folder)):
        if file.startswith(filename):
            p_files.append(os.path.join(folder,file))
    T_list = []
    fig, ax = plt.subplots()
    for p_file in p_files:
        T = None
        M = []
        t = []
        with open(p_file) as f:
                T = float(f.readline())
        if (T not in T_list):
            T_list.append(T)
            with open(p_file) as csvfile:
                lines = csv.reader(csvfile, delimiter=' ')
                sweep = -1
                for row in lines:
                    M.append(float(row[0]))
                    t.append(sweep)
                    sweep += 1
            # first data in txt file is the temperature
            ax.plot(t[1:2000], M[1:2000])
    ax.set_title("|Magnetisation| vs time, L = "+str(L)+", T = "+str(T))
    ax.set_ylabel("|M|")
    ax.set_xlabel("t / sweeps")

    fig.savefig(folder2+"thermalisation_1.png")
    fig.savefig(texfolder+"thermalisation_1.pdf")

def thermalisation_2():
    """
    Same temperature T and lattice sixe L to show randomness of the statistical model
    """
    p_files = []
    dim = 2
    L = 40
    T = "1.500000"
    filename = "magn_vs_time_"+str(dim)+"D_"+str(L)+"_"+T
    for file in sorted(os.listdir(folder)):
        if file.startswith(filename):
            p_files.append(os.path.join(folder,file))

    fig, ax = plt.subplots()
    for p_file in p_files:
        T = None
        M = []
        t = []
        with open(p_file) as csvfile:
            lines = csv.reader(csvfile, delimiter=' ')
            sweep = -1
            for row in lines:
                M.append(float(row[0]))
                t.append(sweep)
                sweep += 1
            # first data in txt file is the temperature
            T = M[0]
        ax.plot(t[1:600], M[1:600])

    ax.set_title("|Magnetisation| vs time, L ="+str(L)+", T = "+str(T))
    ax.set_ylabel("|M|")
    ax.set_xlabel("t / sweeps")
    # ax.legend()
    fig.savefig(folder2+"thermalisation_2.png")
    fig.savefig(texfolder+"thermalisation_2.pdf")


def thermalisation_3():
    """
    Plot at specific temperature but different L
    """
    p_files = []
    T = 1.5
    for file in sorted(os.listdir(folder)):
        if file.endswith(str(T)+"00000.txt"):
            p_files.append(os.path.join(folder,file))

    fig, ax = plt.subplots()
    for p_file in p_files:
        L = (os.path.splitext(os.path.basename(p_file))[0]).split('_',5)[4]
        N = int(L)*int(L)
        M = []
        t = []
        with open(p_file) as csvfile:
            lines = csv.reader(csvfile, delimiter=' ')
            sweep = -1
            for row in lines:
                M.append(float(row[0])/N)
                t.append(sweep)
                sweep += 1
        # first data in txt file is the temperature
        ax.plot(t[1:500], M[1:500], label="L = "+str(L))

    ax.set_title("<|Magnetisation|> vs time, T = "+str(T))
    ax.set_ylabel("m")
    ax.set_xlabel("t / sweeps")
    ax.legend()
    fig.savefig(folder2+"thermalisation_3.png")
    fig.savefig(texfolder+"thermalisation_3.pdf")


def main():
    thermalisation_1()
    thermalisation_2()
    thermalisation_3()
    plt.show()

if (__name__ == '__main__'):
    main()
