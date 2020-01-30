import numpy as np
import matplotlib.pyplot as plt
import csv
import os
from scipy.optimize import curve_fit

def gauss(x, mean, sigma, scale, offset):
    return (scale*np.exp(-((x-mean)/sigma)**2)+offset)

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
    ax.set_ylabel(r"$\tau_e / sweeps$")
    ax.set_xlabel("T / K")
    ax.set_yscale("log")
    ax.legend()
    fig.savefig(folder+"\\tau_e_vs_temp.pdf")
    plt.show()

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
        if file.startswith("autocorr_peak") and not file.endswith(").txt"):
            p_files.append(os.path.join(folder,file))

    fig, axs = plt.subplots(1,2,figsize=(16,6))
    fig.suptitle("Time lag vs Temperature")
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

        axs[0].errorbar(T,tau,tau_err,marker='+',linestyle=" ",label=L)

        # print("T")
        # print(T)
        # print(tau)
        # print(tau_err)
        # left_limit = 40
        # right_limit = len(tau)-left_limit
        if int(L) in [16,24,32]:
            try:
                popt, pcov = curve_fit(gauss, T, tau, sigma=tau_err, absolute_sigma=True, maxfev=5000)
                value_list.append(popt[2]+popt[3])
                error_list.append(np.sqrt((np.diag(pcov)))[2]+np.sqrt((np.diag(pcov)))[3])
                x = np.linspace(T[0],T[-1],100)
                axs[0].plot(x,gauss(x, *popt), label="fit-"+L)
                print(L)
                print("popt :",popt)
                print("pcov^2 :",np.diag(pcov))
                print(value_list[-1],error_list[-1])
            except:
                max_tau = max(tau)
                value_list.append(max_tau)
                error_list.append(tau_err[tau.index(max_tau)])
                print(L,value_list[-1],error_list[-1])
            L_list.append(int(L))


    print(L_list,value_list,error_list)
    axs[1].plot(L_list,value_list,marker='+',linestyle=" ")
    x = np.linspace(10,40,1000)
    axs[1].plot(x,power_fit(x,2.1665,1))

    axs[0].legend()
    for ax in axs:
        ax.set_ylabel(r"$\tau_e$")
        ax.set_xlabel("T")

    fig.savefig(folder+"\\tau_e_vs_temp2.png")

    # fig2, ax2 = plt.subplots()
    # ax2.errorbar(L_list,value_list,error_list,ls='', marker='+', color='Red')

    plt.show()

def main():
    plot_1()
    # plot_2()

if (__name__ == '__main__'):
    main()
