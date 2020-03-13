import numpy as np
import matplotlib.pyplot as plt
import csv
import os

folder = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\data\\n2n_data"
folder2 = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\pngs\\"
texfolder = "C:\\Users\\savva\\OneDrive - University of Cambridge\\Part2\\Computational Project\\report\\texfigures\\"

def plot_1():
    """
    magnetisation vs temperature
    """
    p_files_dict = {}
    R_list = []
    dim = 2
    L_ = 40
    # please change file path when running on different devices
    for file in sorted(os.listdir(folder)):
        if file.endswith(".txt") and not file.endswith(").txt"):
            D = int(((os.path.splitext(os.path.basename(file))[0]).split('_',4)[2])[0])
            L = int((os.path.splitext(os.path.basename(file))[0]).split('_',4)[3])
            R = float((os.path.splitext(os.path.basename(file))[0]).split('_',4)[4])
            if (D == dim) and (L == L_) and (R not in R_list):
                R_list.append(R)
                p_files_dict[R] = os.path.join(folder,file)
    R_list.sort()
    #print(R_list)
    fig, ax = plt.subplots(figsize=(12,8))
    fig1, ax1 = plt.subplots(figsize=(12,8))
    for key in R_list:
        p_file = p_files_dict[key]
        M = []
        errM = []
        E = []
        errE = []
        T = []
        with open(p_file) as csvfile:
            lines = csv.reader(csvfile, delimiter=' ')
            for row in lines:
                T.append(float(row[0]))
                M.append(float(row[1]))
                errM.append(float(row[2]))
                E.append(float(row[3]))
                errE.append(float(row[4]))
        ax.errorbar(T,M,errM,label="R = "+str(key),ls=" ",marker="+")
        ax1.errorbar(T,E,errE,label="R = "+str(key),ls=" ",marker="+")
    ax.legend()
    ax1.legend()


def main():
    plot_1()
    plt.show()

if (__name__ == '__main__'):
    main()
