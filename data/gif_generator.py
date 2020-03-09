import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
import os
from matplotlib.animation import FuncAnimation

folder = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\data\\gif_configs"
folder2 = "C:\\Users\\savva\\Documents\\GitHub\\ising_model_2.0\\pngs\\"

L_ = 100
dim = 2
Temp = 1.5
configs = []

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


fig, ax = plt.subplots()
#fig.patch.set_facecolor('xkcd:mint green')
ax.set_title("L = "+str(L_)+", Temp = "+str(Temp))
im = ax.imshow(np.reshape(configs[0],(L_,L_)), cmap=cm.binary, interpolation='nearest')
ax.invert_yaxis()

def animate(i):
    ax.set_xlabel(i)
    plot = im.set_array(np.reshape(configs[i],(L_,L_)))
    return plot


def main():
    anim = FuncAnimation(fig, animate, frames=len(configs), interval=50)
    anim.save(folder2+"2d_gif_"+str(L_)+"_"+str(Temp)+".gif", writer='imagemagick', savefig_kwargs={'facecolor':'green'})

if (__name__ == '__main__'):
    main()
