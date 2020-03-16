import numpy as np
import matplotlib.pyplot as plt
import csv
import os

folder = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\data\\magn_data"
folder2 = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\pngs\\"
texfolder = "C:\\Users\\savva\\OneDrive - University of Cambridge\\Part2\\Computational Project\\report\\texfigures\\"

def plot_1():
    """
    magnetisation vs temperature with different lattice sizes
    """
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
        avgM = []
        errM = []
        T = []
        with open(p_file) as csvfile:
            lines = csv.reader(csvfile, delimiter=' ')
            for row in lines:
                T.append(float(row[0]))
                avgM.append(float(row[1]))
                errM.append(float(row[2]))
        ax.errorbar(T, avgM, errM, ls='',marker='+', label="L = "+str(L))

    #ax.set_title("<|Magnetisation|> vs Temperature, 2D")
    ax.set_ylabel("m")
    ax.set_xlabel(r"T ($J/k_B$)")
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.subplots_adjust(right=1,top=1)
    ax.legend()
    fig.savefig(folder2+"magn_vs_temp.png")
    fig.savefig(texfolder+"magn_vs_temp.pdf")


def plot_2():
    """
    magnetisation vs temperature with different lattice sizes
    """
    p_files_dict = {}
    L_list = []
    dim = 3
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
    #ax.axvline(2.2692, label="$T_c$", linestyle="--",color="k")
    for L in L_list:
        p_file = p_files_dict[L]
        avgM = []
        errM = []
        T = []
        with open(p_file) as csvfile:
            lines = csv.reader(csvfile, delimiter=' ')
            for row in lines:
                T.append(float(row[0]))
                avgM.append(float(row[1]))
                errM.append(float(row[2]))
        ax.errorbar(T, avgM, errM, ls='',marker='+', label="L = "+str(L))


    #ax.set_title("<|Magnetisation|> vs Temperature, 3D")
    ax.set_ylabel("m")
    ax.set_xlabel(r"T ($J/k_B$)")
    ax.legend()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.subplots_adjust(right=1,top=1)
    fig.savefig(folder2+"magn_vs_temp_3D.png")
    fig.savefig(texfolder+"magn_vs_temp_3D.pdf")

def plot_3():
    """
    magnetisation vs temperature with different lattice sizes
    """
    p_files_dict = {}
    L_list = []
    dim = 4
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
    #ax.axvline(2.2692, label="$T_c$", linestyle="--",color="k")
    for L in L_list:
        p_file = p_files_dict[L]
        avgM = []
        errM = []
        T = []
        with open(p_file) as csvfile:
            lines = csv.reader(csvfile, delimiter=' ')
            for row in lines:
                T.append(float(row[0]))
                avgM.append(float(row[1]))
                errM.append(float(row[2]))
        ax.errorbar(T, avgM, errM, ls='',marker='+', label="L = "+str(L))


    #ax.set_title("<|Magnetisation|> vs Temperature, 4D")
    ax.set_ylabel("m")
    ax.set_xlabel(r"T ($J/k_B$)")
    ax.legend()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.subplots_adjust(right=1,top=1)
    fig.savefig(folder2+"magn_vs_temp_4D.png")
    fig.savefig(texfolder+"magn_vs_temp_4D.pdf")



def main():
    plot_1()
    plot_2()
    plot_3()
    plt.show()

if (__name__ == '__main__'):
    main()
