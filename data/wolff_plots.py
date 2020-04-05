import numpy as np
import matplotlib.pyplot as plt
import csv
import os

folder0 = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\data\\wolff_data\\"
folder2 = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\pngs\\"
texfolder = "C:\\Users\\savva\\OneDrive - University of Cambridge\\Part2\\Computational Project\\report\\texfigures\\"

def plot_1():
    """
    cluster size vs temperature
    """
    folder = folder0+"cluster_size\\"
    p_files_dict = {}
    L_list = []
    dim = 2
    # please change file path when running on different devices
    for file in sorted(os.listdir(folder)):
        if file.endswith(".txt") and not file.endswith(").txt"):
            D = int(((os.path.splitext(os.path.basename(file))[0]).split('_',3)[2])[0])
            if (D != dim):
                continue
            L = int((os.path.splitext(os.path.basename(file))[0]).split('_',3)[3])
            if (L not in L_list):
                L_list.append(L)
                #p_files.append(os.path.join(folder,file))
                p_files_dict[L] = os.path.join(folder,file)
    L_list.sort()
    fig, ax = plt.subplots()
    ax.axvline(2.2692, label="$T_c$", linestyle="--",color="k",alpha=0.7)
    for L in L_list:
        p_file = p_files_dict[L]
        avg_n = []
        err_n = []
        T = []
        with open(p_file) as csvfile:
            lines = csv.reader(csvfile, delimiter=' ')
            for row in lines:
                T.append(float(row[0]))
                avg_n.append(float(row[1]))
                err_n.append(float(row[2]))
        ax.errorbar(T, avg_n, err_n, ls='',marker='+', label="L = "+str(L))

    #ax.set_title("<|Magnetisation|> vs Temperature, 2D")
    ax.set_ylabel(r"$\langle n \rangle$ / N")
    ax.set_xlabel(r"T ($J/k_B$)")
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.subplots_adjust(right=1,top=1)
    ax.legend()
    fig.savefig(folder2+"clustersize_vs_temp.png")
    fig.savefig(texfolder+"clustersize_vs_temp.pdf")

def plot_2():
    """
    Plot tau_e vs temperature
    """
    folder = folder0+"autocorrelation\\"
    p_files_dict = {}
    L_list = []
    dim = 2
    # please change file path when running on different devices
    for file in sorted(os.listdir(folder)):
        if file.endswith(".txt") and not file.endswith(").txt"):
            # D = int(((os.path.splitext(os.path.basename(file))[0]).split('_',3)[2])[0])
            # if (D != dim):
            #     continue
            L = int((os.path.splitext(os.path.basename(file))[0]).split('_',3)[2])
            if (L not in L_list):
                L_list.append(L)
                #p_files.append(os.path.join(folder,file))
                p_files_dict[L] = os.path.join(folder,file)
    L_list.sort()
    fig, ax = plt.subplots()
    plt.subplots_adjust(right=1,top=1)
    ax.axvline(2.2692, label="$T_c$", linestyle="--",color="k",alpha=0.5)
    for L in L_list:
        p_file = p_files_dict[L]
        avg_tau = []
        err_tau = []
        T = []
        with open(p_file) as csvfile:
            lines = csv.reader(csvfile, delimiter=' ')
            for row in lines:
                T.append(float(row[0]))
                avg_tau.append(float(row[1]))
                err_tau.append(float(row[2]))
        ax.errorbar(T, avg_tau, err_tau, ls='',marker='+', label="L = "+str(L))
    #ax.set_title("Time lag vs Temperature")
    ax.set_ylabel(r"$\tau_e$ (1 WS)")
    ax.set_xlabel(r"T ($J/k_B$)")
    #ax.set_yscale("log")
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.legend()
    fig.savefig(folder2+"wolff_autorr_vs_temp.png")
    fig.savefig(texfolder+"wolff_autocorr_vs_temp.pdf")

def main():
    #plot_1()
    plot_2()
    plt.show()

if (__name__ == '__main__'):
    main()
