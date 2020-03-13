import numpy as np
import matplotlib.pyplot as plt
import csv
import os

folder = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\data\\H_data"
folder2 = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\pngs\\"
texfolder = "C:\\Users\\savva\\OneDrive - University of Cambridge\\Part2\\Computational Project\\report\\texfigures\\"

def plot_1():
    """
    magnetisation vs external field at different temperatures
    """
    p_files_dict = {}
    L_ = 36
    dim = 2
    # please change file path when running on different devices
    for file in sorted(os.listdir(folder)):
        if file.endswith(".txt") and not file.endswith(").txt"):
            D = int(((os.path.splitext(os.path.basename(file))[0]).split('_',4)[2])[0])
            L = int((os.path.splitext(os.path.basename(file))[0]).split('_',4)[3])
            T = (os.path.splitext(os.path.basename(file))[0]).split('_',4)[4]
            #print(D,L,T)
            if (D != dim) or (L != L_):
                continue
            if T not in p_files_dict:
                p_files_dict[T] = os.path.join(folder,file)
    #print(p_files_dict)
    fig, axs = plt.subplots(1,2,figsize=(16,8),sharey=True,gridspec_kw={'hspace': 0, 'wspace': 0.03})
    for T in p_files_dict:
        p_file = p_files_dict[T]
        H = []
        M = []
        errM = []
        with open(p_file) as csvfile:
                lines = csv.reader(csvfile, delimiter=' ')
                for row in lines:
                    H.append(float(row[0]))
                    M.append(float(row[1]))
                    errM.append(float(row[2]))
        if (float(T) < 2.28):
            axs[0].errorbar(H,M,errM,label="T = "+T[0:4],ls="",marker="+",alpha=0.5)
        else:
            axs[1].errorbar(H,M,errM,label="T = "+T[0:4],ls="",marker="+",alpha=0.5)
        #ax.scatter(H,M,alpha=0.5)
    fig.suptitle(r"Magnetisation vs External Filed, L=36")
    axs[0].annotate("metastable states",
            xy=(0.5,-1300), xycoords='data',
            xytext=(1,-1100), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )
    axs[0].annotate("discontinuities",
            xy=(-0.01,50), xycoords='data',
            xytext=(-1.5,500), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )
    axs[0].annotate("",
            xy=(-0.63,-500), xycoords='data',
            xytext=(-1.1,480), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )
    axs[1].annotate("no remnant\nmagnetisation",
            xy=(-0.01,10), xycoords='data',
            xytext=(-1,300), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )
    for ax in axs:
        # put axis in mid
        ax.spines['left'].set_position('center')
        ax.spines['bottom'].set_position('center')
        ax.spines['left'].set_alpha(0.2)
        # Eliminate upper and right axes
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.set_ylabel("M",rotation=0)
        ax.set_xlabel(r"H / $J/k_B$")
        ax.xaxis.set_label_coords(0.975, 0.475)
        ax.yaxis.set_label_coords(0.475, 1.0)
        ax.legend()

    fig.savefig(folder2+"magn_vs_field.png")
    fig.savefig(texfolder+"magn_vs_field.pdf")

def plot_2():
    """
    energy vs external field at different temperatures
    """
    p_files_dict = {}
    L_ = 36
    dim = 2
    # please change file path when running on different devices
    for file in sorted(os.listdir(folder)):
        if file.endswith(".txt") and not file.endswith(").txt"):
            D = int(((os.path.splitext(os.path.basename(file))[0]).split('_',4)[2])[0])
            L = int((os.path.splitext(os.path.basename(file))[0]).split('_',4)[3])
            T = (os.path.splitext(os.path.basename(file))[0]).split('_',4)[4]
            if (D != dim) or (L != L_):
                continue
            if T not in p_files_dict:
                p_files_dict[T] = os.path.join(folder,file)
    #print(p_files_dict)
    fig, axs = plt.subplots(1,2,figsize=(16,8),sharey=True,gridspec_kw={'hspace': 0, 'wspace': 0.03})
    for T in p_files_dict:
        p_file = p_files_dict[T]
        H = []
        E = []
        errE = []
        with open(p_file) as csvfile:
                lines = csv.reader(csvfile, delimiter=' ')
                for row in lines:
                    H.append(float(row[0]))
                    E.append(float(row[3]))
                    errE.append(float(row[4]))
        if (float(T) < 2.28):
            axs[0].errorbar(H,E,errE,label="T = "+T[0:4],ls="",marker="+",alpha=0.5)
        else:
            axs[1].errorbar(H,E,errE,label="T = "+T[0:4],ls="",marker="+",alpha=0.5)
        #ax.scatter(H,E,alpha=0.5)
    fig.suptitle(r"Energy vs External Filed, L=36")
    axs[0].annotate("metastable\n states",
            xy=(0.4,-2000), xycoords='data',
            xytext=(-1.5,-1000), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )
    axs[0].annotate("",
            xy=(-0.55,-1800), xycoords='data',
            xytext=(-1.1,-1100), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )
    axs[0].annotate("discontinuity",
            xy=(-0.6,-2900), xycoords='data',
            xytext=(-1.5,-2500), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )
    for ax in axs:
        # put axis in mid
        ax.spines['left'].set_position('center')
        ax.spines['left'].set_alpha(0.2)
        #ax.spines['bottom'].set_position('center')
        # Eliminate upper and right axes
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.set_ylabel("E",rotation=0)
        ax.set_xlabel(r"H / $J/k_B$")
        #ax.xaxis.set_label_coords(0.975, 0.475)
        ax.yaxis.set_label_coords(0.475, 1.0)
        ax.legend()
    fig.savefig(folder2+"energy_vs_field.png")
    fig.savefig(texfolder+"energy_vs_field.pdf")

def main():
    plot_1()
    plot_2()
    plt.show()

if (__name__ == '__main__'):
    main()
