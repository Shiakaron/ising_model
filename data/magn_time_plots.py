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
    L = 80
    filename = "magn_vs_time_"+str(dim)+"D_"+str(L)
    for file in sorted(os.listdir(folder)):
        if file.startswith(filename) and not file.endswith(").txt"):
            p_files.append(os.path.join(folder,file))
    T_list = []
    fig, ax = plt.subplots()
    for p_file in p_files:
        T = None
        M = []
        t = []
        with open(p_file) as f:
            T = float((os.path.splitext(os.path.basename(p_file))[0]).split('_',5)[5])
        if (T not in T_list):
            T_list.append(T)
            with open(p_file) as csvfile:
                lines = csv.reader(csvfile, delimiter=' ')
                sweep = 0
                for row in lines:
                    M.append(float(row[0]))
                    t.append(sweep)
                    sweep += 1
            # first data in txt file is the temperature
            #ax.plot(t[1:2000], M[1:2000])
            ax.plot(t[:3000],M[:3000],label="T = "+str(T),alpha=0.8)
    # ax.set_title("|Magnetisation| vs time, L = "+str(L)+", T = "+str(T))
    ax.set_ylabel("m")
    ax.set_xlabel("t (MCS)")
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.legend()
    plt.subplots_adjust(right=1,top=1)
    fig.savefig(folder2+"thermalisation_1.png")
    fig.savefig(texfolder+"thermalisation_1.pdf")

def thermalisation_2():
    """
    Same temperature T and lattice sixe L to show randomness of the statistical model
    """
    p_files = []
    dim = 2
    L = 80
    T = "1.500000"
    filename = "magn_vs_time_"+str(dim)+"D_"+str(L)+"_"+T
    for file in sorted(os.listdir(folder)):
        if file.startswith(filename):
            p_files.append(os.path.join(folder,file))

    fig, ax = plt.subplots()
    for p_file in p_files:
        M = []
        t = []
        with open(p_file) as csvfile:
            lines = csv.reader(csvfile, delimiter=' ')
            sweep = 0
            for row in lines:
                M.append(float(row[0]))
                t.append(sweep)
                sweep += 1
        if (M[0]==1):
            ax.plot(t[:2000], M[:2000],color="b")
        else:
            ax.plot(t[:2000], M[:2000],color="r",alpha=0.7)

    # ax.set_title("|Magnetisation| vs time, L ="+str(L)+", T = "+str(T))
    ax.set_ylabel("m")
    ax.set_xlabel("t (MCS)")
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.subplots_adjust(right=1,top=1)
    fig.savefig(folder2+"thermalisation_2.png")
    fig.savefig(texfolder+"thermalisation_2.pdf")

def thermalisation_3():
    """
    Plot at specific temperature but different L
    """
    p_files_dict = {}
    L_list = []
    T = 1.5
    for file in sorted(os.listdir(folder)):
        if file.endswith(str(T)+"00000.txt"):
            L = int((os.path.splitext(os.path.basename(file))[0]).split('_',5)[4])
            L_list.append(L)
            p_files_dict[L] = os.path.join(folder,file)
    L_list.sort()
    fig, ax = plt.subplots()
    for L in L_list:
        p_file = p_files_dict[L]
        M = []
        t = []
        with open(p_file) as csvfile:
            lines = csv.reader(csvfile, delimiter=' ')
            sweep = 0
            for row in lines:
                M.append(float(row[0]))
                t.append(sweep)
                sweep += 1
        # first data in txt file is the temperature
        ax.plot(t, M, label="L = "+str(L),alpha=0.9)

    # ax.set_title("<|Magnetisation|> vs time, T = "+str(T))
    ax.set_ylabel("m")
    ax.set_xlabel("t (MCS)")
    ax.legend()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.subplots_adjust(right=1,top=1)
    fig.savefig(folder2+"thermalisation_3.png")
    fig.savefig(texfolder+"thermalisation_3.pdf")

def main():
    thermalisation_1()
    thermalisation_2()
    thermalisation_3()
    plt.show()

if (__name__ == '__main__'):
    main()
