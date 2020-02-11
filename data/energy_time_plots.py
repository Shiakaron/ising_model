import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import itertools

folder = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\data\\energy_data\\vs_time"
folder2 = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\pngs\\"
texfolder = "C:\\Users\\savva\\OneDrive - University of Cambridge\\Part2\\Computational Projects\\ising_model\\Report\\texfigures\\"

def plot_1():
    """
    plot <E> vs time for L=80 at different temperatures
    """
    p_files = []
    filename = "energy_data_2D_80"
    for file in sorted(os.listdir(folder)):
        if file.startswith(filename):
            p_files.append(os.path.join(folder,file))
    T_list = []
    fig, ax = plt.subplots()
    for p_file in p_files[3::3]:
        T = (os.path.splitext(os.path.basename(p_file))[0]).split('_',4)[4]
        #print(T)
        E = []
        t = []
        if (T not in T_list):
            T_list.append(T)
            with open(p_file) as csvfile:
                lines = csv.reader(csvfile, delimiter=' ')
                sweep = 0
                for row in lines:
                    E.append(float(row[0]))
                    t.append(sweep)
                    sweep += 1
            ax.plot(t[0:200], E[0:200],label="T = "+format(T[0:3]))
    ax.set_title("Energy per bond vs Time")
    ax.set_ylabel("e / J")
    ax.set_xlabel("t / sweeps")
    ax.legend()

    fig.savefig(folder2+"energy_vs_time.png")
    fig.savefig(texfolder+"energy_vs_time.pdf")

def main():
    plot_1()
    plt.show()

if (__name__ == '__main__'):
    main()
