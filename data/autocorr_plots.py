import numpy as np
import matplotlib.pyplot as plt
import csv
import os
from scipy.optimize import curve_fit

folder2 = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\pngs\\"
texfolder = "C:\\Users\\savva\\OneDrive - University of Cambridge\\Part2\\Computational Project\\report\\texfigures\\"

def gauss(x, mean, sigma, scale, offset):
    return (scale*np.exp(-((x-mean)/sigma)**2)+offset)

def gauss2(x, mean, sigma, scale):
    return (scale*np.exp(-((x-mean)/sigma)**2))

def power_fit(x,pow,scale,offset):
    return (scale*np.power(x,pow)+offset)

def plot_1():
    """
    Plot tau_e vs temperature
    """
    p_files = []
    folder = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\data\\autocorrelation_data\\initial_investigation"
    for file in sorted(os.listdir(folder)):
        if file.startswith("autocorr_times") and file.endswith("L.txt"):
            p_files.append(os.path.join(folder,file))

    fig, ax = plt.subplots()
    ax.axvline(2.2692, label="$T_c$", linestyle="--",color="k")
    for p_file in p_files:
        L = (os.path.splitext(os.path.basename(p_file))[0]).split('_',4)[2]
        T_list = []
        tau_list = []
        with open(p_file) as csvfile:
            lines = csv.reader(csvfile, delimiter=' ')
            for row in lines:
                T_list.append(float(row[0]))
                tau_list.append(float(row[1]))
        ax.plot(T_list,tau_list,marker='+',label="L = "+str(L))
    ax.set_title("Time lag vs Temperature")
    ax.set_ylabel(r"$\tau_e$ / sweeps")
    ax.set_xlabel(r"T / $J/k_B$")
    ax.set_yscale("log")
    ax.legend()

    fig.savefig(folder2+"tau_e_vs_temp.png")
    fig.savefig(texfolder+"tau_e_vs_temp.pdf")

def plot_2():
    """
    1. fit gauss/lorentz on tau vs T for different L to get tau_peak
    2. plot of tau_peak vs L
    """
    L_list = []
    value_list = []
    error_list = []
    folder="C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\data\\autocorrelation_data\\peak_investigation"
    p_files = []
    for file in sorted(os.listdir(folder)):
        if file.startswith("autocorr_peak") and file.endswith("L.txt"):
            p_files.append(os.path.join(folder,file))

    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    # fig.suptitle("Time lag vs Temperature")

    for p_file in p_files:
        L = (os.path.splitext(os.path.basename(p_file))[0]).split('_',4)[2]
        T = []
        tau = []
        tau_err = []
        with open(p_file) as csvfile:
            lines = csv.reader(csvfile, delimiter=' ')
            for row in lines:
                T.append(float(row[0]))
                tau.append(float(row[1]))
                tau_err.append(float(row[2]))

        if int(L) in [16,20,24,28,32,36,40]:
            ax1.errorbar(T,tau,tau_err,marker='+',linestyle=" ",label=L)
            popt, pcov = curve_fit(gauss2, T, tau, sigma=tau_err, absolute_sigma=True, maxfev=5000, p0=[2.3, 0.1, 100], bounds=((2.28,-np.inf,-np.inf),(2.32,np.inf,np.inf)))
            value_list.append(popt[2])
            error = abs(np.sqrt(abs((np.diag(pcov)))[2]))
            error_list.append(error)
            x = np.linspace(T[0],T[-1],100)
            ax1.plot(x,gauss2(x, *popt), color="c",linewidth=1)
            L_list.append(int(L))
            # print(L_list[-1],value_list[-1],error_list[-1])

    ax1.legend()
    ax1.set_title(r"Gaussian Fit on Time lag vs Temperature")
    ax1.set_ylabel(r"$\tau_e$ / sweeps")
    ax1.set_xlabel(r"T / $J/k_B$")

    ax2.set_title(r"Power fit on peak Time lag vs Lattice size")
    ax2.set_ylabel(r"Peak $\tau_e$ / sweeps")
    ax2.set_xlabel("L")
    ax2.errorbar(L_list,value_list,error_list,marker='+',linestyle=" ")
    popt2, pcov2 = curve_fit(power_fit, L_list, value_list, sigma=error_list, absolute_sigma=True, maxfev=2000, p0=[1,1,0])
    x2 = np.linspace(L_list[0],L_list[-1],100)
    ax2.plot(x2,power_fit(x2,*popt2), label = "pow = " + '{:.2f}'.format(popt2[0]) + r"$ \pm $" + '{:.2f}'.format(abs(np.sqrt(abs((np.diag(pcov2)))[0]))))
    print("power = ",popt2[0],abs(np.sqrt(abs((np.diag(pcov2)))[0])))


    # popt_plus = popt2 + np.sqrt((np.diag(pcov2)))
    # print(popt2,"+",np.sqrt((np.diag(pcov2))),"=",popt_plus)
    # ax2.plot(x2,power_fit(x2,*popt_plus), label = "pow = " + '{:.2f}'.format(popt_plus[0]),linestyle="--",color="k")
    #
    # popt_minus = popt2 - np.sqrt((np.diag(pcov2)))
    # print(popt2,"+",np.sqrt((np.diag(pcov2))),"=",popt_minus)
    # ax2.plot(x2,power_fit(x2,*popt_minus), label = "pow = " + '{:.2f}'.format(popt_minus[0]),linestyle="--",color="k")

    ax2.legend()
    fig1.savefig(folder2+"tau_e_peak_vs_temp.png")
    fig2.savefig(folder2+"tau_e_peak_vs_L.png")
    fig1.savefig(texfolder+"\\tau_e_peak_vs_temp.pdf")
    fig2.savefig(texfolder+"\\tau_e_peak_vs_L.pdf")


def main():
    plot_1()
    plot_2()
    plt.show()

if (__name__ == '__main__'):
    main()
