import numpy as np
import matplotlib.pyplot as plt
import csv
import os
from scipy.optimize import curve_fit

def plot_0(L,T,limit):
    """
    
    """
    p_file = None
    folder = "autocorrelation_data"
    filename = "autocorr_data_"+str(L)+"_L_"+str("{0:.6f}".format(T))+"_T.txt"
    for file in sorted(os.listdir(folder)):
        if file.endswith(filename):
            p_file = os.path.join(folder,file)

    M = []
    autoc = []
    t = []
    fig, ax = plt.subplots()
    ax.axhline(np.exp(-1), label="$e^{-1}$", linestyle="--",color="k")
    with open(p_file) as csvfile:
        lines = csv.reader(csvfile, delimiter=" ")
        sweep = 0
        for row in lines:
            M.append(float(row[0]))
            autoc.append(float(row[1]))
            t.append(sweep)
            sweep += 1
            if (sweep==limit):
                break
    ax.plot(t,M,label="M")
    ax.plot(t,autoc,label="Autocorrelation")
    ax.set_title("Autocorrelation and <|M|> vs time \nL = "+str(L)+ ", T = "+str(T))
    ax.set_ylabel(r"<|M|> and A($\tau$)")
    ax.set_xlabel("t/sweeps")
    ax.legend()
    figname = "MnA_vs_t_"+str(L)+"_L_"+str(T)+"_T.png"
    fig.savefig(folder+"\\"+figname)
    plt.show()

def plot_1():
    """
    Plot tau_e vs temperature
    """
    p_files = []
    folder = "autocorrelation_data"
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
    ax.set_ylabel(r"$\tau_e$")
    ax.set_xlabel("T")
    ax.set_yscale("log")
    ax.legend()
    fig.savefig(folder+"\\tau_e_vs_temp.png")
    plt.show()

def gauss(x, mean, sigma, scale, offset):
    return (scale*np.exp(-((x-mean)/sigma)**2)+offset)

def power_fit(x,pow,scale):
    return (scale*np.power(x,pow))

def plot_2():
    """
    1. fit gauss/lorentz on tau vs T for each L to identify peak
    2. plot tau peak vs L
    """
    L_list = []
    value_list = []
    error_list = []

    folder = "autocorrelation_data"
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
    # plot_0(32,2.25,400)
    # plot_1()
    plot_2()

if (__name__ == '__main__'):
    main()
