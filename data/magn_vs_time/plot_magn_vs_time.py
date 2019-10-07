import numpy as np
import matplotlib.pyplot as plt
import csv
import os

def plot_1():
    """
    Plot for same temperature T and lattice size L to show the randomness of the
    statistical model.
    """
    p_files = []
    dim = 2
    L = 64
    T = "1.500000"
    filename = "magn_vs_time_"+str(dim)+"D_"+str(L)+"_"+T
    folder = "."
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
        ax.plot(t[1:500], M[1:500], label="T = "+str(T))

    ax.set_title("|Magnetisation| vs time(sweeps)")
    ax.set_ylabel("|M|")
    ax.set_xlabel("t")
    ax.legend()
    fig.savefig(folder+"\\magn_vs_time_rand.png")
    plt.show()

def plot_2():
    """
    Similar plot but varying the temperature this time
    """
    p_files = []
    dim = 2
    L = 64
    filename = "magn_vs_time_"+str(dim)+"D_"+str(L)
    folder = "."
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
            ax.plot(t[1:500], M[1:500], label="T = "+str(T))

    ax.set_title("|Magnetisation| vs time(sweeps)")
    ax.set_ylabel("|M|")
    ax.set_xlabel("t")
    ax.legend()
    fig.savefig(folder+"\\magn_vs_time_temp.png")
    plt.show()

def main():
    plot_1()
    plot_2()

if (__name__ == '__main__'):
    main()
