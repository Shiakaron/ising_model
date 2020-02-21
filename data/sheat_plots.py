import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import itertools
from scipy.optimize import curve_fit

folder = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\data\\specific_heat"
folder2 = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\pngs\\"
texfolder = "C:\\Users\\savva\\OneDrive - University of Cambridge\\Part2\\Computational Projects\\ising_model\\Report\\texfigures\\"

def gauss(x, mean, sigma, scale):
    return (scale*np.exp(-((x-mean)/sigma)**2))

def plot_1():
    """
    specific heat vs temperature with different L's
    """
    p_files = []
    dim = 2
    for file in sorted(os.listdir(folder)):
        if file.endswith(".txt") and not file.endswith(").txt") and file.startswith("specific_heat_data"):
            p_files.append(os.path.join(folder,file))

    L_list = []
    fig, ax = plt.subplots(figsize=(12,8))
    ax.axvline(2.2692, label="$T_c$", linestyle="--",color="k")
    marker = itertools.cycle(('*', '+', '.', ',', 'o'))
    for p_file in p_files:
        L = (os.path.splitext(os.path.basename(p_file))[0]).split('_',4)[4]
        avgC = []
        errC = []
        T = []
        if (L not in L_list):
            L_list.append(L)
            with open(p_file) as csvfile:
                lines = csv.reader(csvfile, delimiter=' ')
                for row in lines:
                    T.append(float(row[0]))
                    avgC.append(float(row[1]))
                    errC.append(float(row[2]))
            ax.errorbar(T, avgC, errC, ls='',marker = next(marker), label="L = "+str(L))

    ax.set_title("Specific heat vs Temperature")
    ax.set_ylabel(r"C / $k_B$")
    ax.set_xlabel(r"T / $J/k_B$")
    ax.set_yscale("log")
    ax.legend()

    fig.savefig(texfolder+"c_vs_temp.pdf")
    fig.savefig(folder2+"c_vs_temp.png")

def plot_2():
    """
    specific heat vs temperature with different L's
    """
    p_files = []
    dim = 2
    for file in sorted(os.listdir(folder)):
        if file.endswith(".txt") and not file.endswith(").txt") and file.startswith("specific_heat_peak"):
            p_files.append(os.path.join(folder,file))

    L_list = []
    fig, ax = plt.subplots(figsize=(12,8))
    ax.axvline(2.2692, label="$T_c$", linestyle="--",color="k")
    marker = itertools.cycle(('*', '+', '.', ',', 'o'))
    for p_file in p_files:
        L = (os.path.splitext(os.path.basename(p_file))[0]).split('_',4)[4]
        avgC = []
        errC = []
        T = []
        if (L not in L_list):
            L_list.append(L)
            with open(p_file) as csvfile:
                lines = csv.reader(csvfile, delimiter=' ')
                for row in lines:
                    T.append(float(row[0]))
                    avgC.append(float(row[1]))
                    errC.append(float(row[2]))
            ax.errorbar(T, avgC, errC, ls='',marker = next(marker), label="L = "+str(L))
            popt, pcov = curve_fit(gauss, T, avgC, sigma=errC, absolute_sigma=True, maxfev=5000, p0=[2.3, 0.1, 100], bounds=((2.25,-np.inf,-np.inf),(2.4,np.inf,np.inf)))
            # value_list.append(popt[2])
            # error = abs(np.sqrt(abs((np.diag(pcov)))[2]))
            # error_list.append(error)
            x = np.linspace(T[0],T[-1],100)
            ax.plot(x,gauss2(x, *popt), color="c",linewidth=1)
            print(L,popt,np.sqrt(np.diag(pcov)))

    ax.set_title("Specific heat vs Temperature around critical")
    ax.set_ylabel(r"C / $k_B$")
    ax.set_xlabel(r"T / $J/k_B$")
    ax.set_yscale("log")
    ax.legend()

    #fig.savefig(texfolder+"c_vs_temp_peak.pdf")
    #fig.savefig(folder2+"c_vs_temp_peak.png")

def main():
    #plot_1()
    plot_2()
    plt.show()

if (__name__ == '__main__'):
    main()
