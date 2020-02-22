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

def linear(x, nu, a):
    """
    ln(T_c(N) - T_c(inf)) = (1/nu) * ln(L) + ln(a)
    y = (1/nu) * x + ln(a)
    """
    return (1/nu) * x + np.log(a)

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
    and
    peak vs L
    """
    p_files = []
    dim = 2
    for file in sorted(os.listdir(folder)):
        if file.endswith(".txt") and not file.endswith(").txt") and file.startswith("specific_heat_peak"):
            p_files.append(os.path.join(folder,file))

    L_list = []
    ln_L_list = []
    y_list = []
    y_err_list = []
    T_c_inf = 2/np.log(1+np.sqrt(2))
    fig, ax = plt.subplots(figsize=(12,8))
    ax.axvline(T_c_inf, label="$T_c$", linestyle="--",color="k")
    marker = itertools.cycle(('*', '+', '.', ',', 'o'))
    limits = {
        "8":[0,-1],
        "16":[0,-1],
        "24":[3,-3],
        "32":[5,15],
        "36":[0,-1],
        "40":[0,-1],
        "44":[0,-1],
        "48":[5,20]
    }
    for key in limits:
        L = int(key)
        for p_file in p_files:
            if p_file.endswith("_"+key+".txt"):
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
                    left = limits[key][0]
                    right = limits[key][1]
                    T_fit = T[left:right]
                    C_fit = avgC[left:right]
                    err_fit = errC[left:right]
                    popt, pcov = curve_fit(gauss, T_fit, C_fit, sigma=err_fit, absolute_sigma=True, maxfev=5000, p0=[2.3, 0.1, 1000], bounds=((2.25,-5,-np.inf),(2.4,5,np.inf)))
                    # value_list.append(popt[2])
                    # error = abs(np.sqrt(abs((np.diag(pcov)))[2]))
                    # error_list.append(error)
                    x = np.linspace(T[left],T[right],100)
                    ax.plot(x,gauss(x, *popt), color="c",linewidth=1)
                    #print(L,popt,np.sqrt(np.diag(pcov)))
                    ln_L_list.append(np.log(int(L)))
                    y_list.append(np.log(popt[0] - T_c_inf))
                    y_err_list.append(np.sqrt(np.diag(pcov)[0]))
                    print(L, ln_L_list[-1],y_list[-1],y_err_list[-1])

    ax.set_title("Specific heat vs Temperature around critical")
    ax.set_ylabel(r"C / $k_B$")
    ax.set_xlabel(r"T / $J/k_B$")
    ax.set_yscale("log")
    ax.legend()

    #fig.savefig(texfolder+"c_vs_temp_peak.pdf")
    #fig.savefig(folder2+"c_vs_temp_peak.png")

    fig2, ax2 = plt.subplots(figsize=(12,8))
    ax2.errorbar(ln_L_list, y_list, y_err_list,ls="",marker='+')
    popt2, pcov2 = curve_fit(linear, ln_L_list, y_list, sigma=y_err_list, absolute_sigma=True, maxfev=5000, p0=[-1, 1], bounds=((-np.inf,0.000001),(np.inf,np.inf)))
    x2 = np.linspace(ln_L_list[0],ln_L_list[-1],100)
    ax2.plot(x2,linear(x2, *popt2), color="c",linewidth=1)
    print(popt2,np.sqrt(np.diag(pcov2)))

def main():
    #plot_1()
    plot_2()
    plt.show()

if (__name__ == '__main__'):
    main()
