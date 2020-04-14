import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import itertools
from scipy.optimize import curve_fit

# Define the directories to read data from and to output figures
currentworkingdirectory = os.getcwd()
folder = currentworkingdirectory+"\\data\\chi_data"
folder2 = currentworkingdirectory+"\\pngs\\"

def gauss(x, mean, sigma, scale):
    """
    Gaussian function
    """
    return (scale*np.exp(-((x-mean)/sigma)**2))

def for_T_c_fun(L,T_c_inf,a,nu):
    """
    Finite Size scaling function
    """
    return T_c_inf + a * np.power(L,-1/nu)

def plot_1():
    """
    Plot Susceptibility vs Temperature with different L's to get psuedo-critical Temperatures.
    Then following that we fit the function for_T_c_fun to get the critical Temperature
    """
    # Reading the files in the folder. Please change file path when running on different devices
    p_files_dict = {} # save the file paths in a dictionary with the lattice size L as the key
    L_list = [] # to sort later in order to be able to retrieve file paths in ascending order of lattice size
    dim = 2
    for file in sorted(os.listdir(folder)):
        if file.endswith(".txt") and not file.endswith(").txt"):
            D = int(((os.path.splitext(os.path.basename(file))[0]).split('_',4)[3])[0])
            L = int((os.path.splitext(os.path.basename(file))[0]).split('_',4)[4])
            if (D == dim) and (L not in L_list):
                L_list.append(L)
                p_files_dict[L] = os.path.join(folder,file)
    L_list.sort()

    # define T_c and lists to input pseudo-critical temps (T_c_N) and their errors (T_c_N_err)
    T_c_N_list = []
    T_c_N_err_list = []
    T_c_inf = 2/np.log(1+np.sqrt(2)) # Onsager's solution
    # to adjust fitting limits so that the peak is captured as best as possible from the copmuted data
    limits = {
        8:[0,-1],
        12:[0,-1],
        16:[0,-1],
        20:[8,-5],
        24:[8,-7],
        28:[8,-10],
        32:[7,-2],
        36:[7,-3],
        40:[7,-4],
        44:[7,-4],
        48:[9,-7],
        52:[8,-10],
        56:[8,-10],
        64:[3,-5]
    }
    # iterate over the lattice sizes conveniently defined in the limits dictionary
    for key in limits:
        p_file = p_files_dict[key] # retrieve file path
        # initialise lists
        avgChi = []
        errChi = []
        T = []
        # read the file
        with open(p_file) as csvfile:
            lines = csv.reader(csvfile, delimiter=' ')
            for row in lines:
                T.append(float(row[0]))
                avgChi.append(float(row[1]))
                errChi.append(float(row[2]))
        # retrieve limits previously set and define the fitting regions as [left:right]
        left = limits[key][0]
        right = limits[key][1]
        T_fit = T[left:right]
        Chi_fit = avgChi[left:right]
        errChi_fit = errChi[left:right]
        # fit peak with Gaussian to get value for T_c_N ("mean" parameter in Gaussian) and it's error
        popt, pcov = curve_fit(gauss, T_fit, Chi_fit, sigma=errChi_fit, absolute_sigma=True, maxfev=5000, p0=[2.3, 0.1, 1000], bounds=((2.25,-5,-np.inf),(2.6,5,np.inf)))
        # append the results in the previously initialised lists
        T_c_N_list.append(popt[0])
        T_c_N_err_list.append(np.sqrt(np.diag(pcov)[0]))
        # Plot to check fit and adjust limits in order to get the best possible fit on the peak.
        L_plot = 40 # change this to view each fit individually
        if (key == L_plot):
            fig, ax = plt.subplots(figsize=(12,8))
            ax.axvline(T_c_inf, label="$T_c$", linestyle="--",color="k",alpha=0.5)
            x = np.linspace(T[left],T[right],100)
            ax.errorbar(T, avgChi, errChi, ls='',marker = "+", label="L = "+str(key))
            ax.plot(x,gauss(x, *popt), color="c",linewidth=1)
            ax.set_title("Specific heat vs Temperature around critical")
            ax.set_ylabel(r"$\chi$")
            ax.set_xlabel(r"T / $J/k_B$")
            ax.legend()
    # fit the results with the Finite Size scaling function to get value for T_c and the critical exponent nu
    popt3, pcov3 = curve_fit(for_T_c_fun, L_list, T_c_N_list, sigma=T_c_N_err_list, absolute_sigma=True, maxfev=5000, p0=[2.26,1,1], bounds=((0,-np.inf,0.000001),(np.inf,np.inf,np.inf)))
    # Plot the results
    fig3, ax3 = plt.subplots()
    plt.subplots_adjust(right=1,top=1)
    ax3.errorbar(L_list, T_c_N_list, T_c_N_err_list,ls="",marker='+')
    x3 = np.linspace(L_list[0],L_list[-1],100)
    ax3.plot(x3,for_T_c_fun(x3, *popt3), color="k",linewidth=1)
    ax3.set_ylabel(r"$T_c$(L) (J/$k_B$)")
    ax3.set_xlabel("L")
    ax3.spines['right'].set_color('none')
    ax3.spines['top'].set_color('none')
    # add the resulting T_c and nu copmuted on the plot
    s1 = r"$T_c$ = "+str('{:.3f}'.format(popt3[0]))+r"$\pm$"+str('{:.3f}'.format(np.sqrt(np.diag(pcov3)[0])))
    s2 = r"$\nu$ = "+str('{:.2f}'.format(popt3[2]))+r"$\pm$"+str('{:.2f}'.format(np.sqrt(np.diag(pcov3)[2])))
    s = s1+"\n"+s2
    ax3.text(35,2.4,s)
    # save figure in folder2
    fig3.savefig(folder2+"chi_check_Onsager.png")

def main():
    plot_1()
    plt.show()

if (__name__ == '__main__'):
    main()
