import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
import os
from matplotlib.animation import FuncAnimation

folder = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\data\\configs\\figure"
folder2 = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\pngs\\"
texfolder = "C:\\Users\\savva\\OneDrive - University of Cambridge\\Part2\\Computational Project\\report\\texfigures\\"

L_ = 200
dim = 2
#Temp = 2.28
#Temp = 1.0
Temp = 3.0
configs = []

def get_config():
    for file in sorted(os.listdir(folder)):
        if file.endswith(".txt") and not file.endswith(").txt"):
            split = (os.path.splitext(os.path.basename(file))[0]).split('_',3)
            D = int((split[1])[0])
            L = int(split[2])
            T = float(split[3])
            if (D == dim) and (L == L_) and (Temp == T):
                p_file = os.path.join(folder,file)
                with open(p_file) as csvfile:
                    lines = csv.reader(csvfile, delimiter=' ')
                    for row in lines:
                        config = []
                        for item in row:
                            config.append(int(item))
                        configs.append(config)

def plot_1():
    fig, ax = plt.subplots(figsize=(10,10))
    #fig.patch.set_facecolor('xkcd:mint green')
    #ax.set_title("L = "+str(L_)+", Temp = "+str(Temp))
    ax.axis("off")
    ax.imshow(np.reshape(configs[0],(L_,L_)), cmap=cm.binary, interpolation='nearest')
    plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.invert_yaxis()
    fig.savefig(texfolder+"lattice_"+str(dim)+"D_"+str(L_)+"_"+str(Temp)+".pdf")
    fig.savefig(folder2+"lattice_"+str(dim)+"D_"+str(L_)+"_"+str(Temp)+".png")

def main():
    get_config()
    plot_1()


if (__name__ == '__main__'):
    main()
