import imp
from Nano_Play.nanite import *
from ..structure.tools import *
from  ..structure import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
mpl.rcParams['axes.linewidth']=2.0
#mpl.rcParams["font.family"] = "Times New Roman"
#mpl.rcParams["mathtext.fontset"] = "cm"
mpl.rcParams["font.style"] = "normal"
mpl.rcParams["font.size"]=10
mpl.rcParams["legend.frameon"] = False
mpl.rcParams["legend.fancybox"] = False

def distance(traj,pairs,width=2.0,fig_ext='',format='.png'):
    """Plots the variation of distances for a list of pair of Atoms in a trajectory"""
    for pair in pairs:
        bn,dist=trajectory.dist_var(traj.frames,pair)
        plt.plot(traj.time,dist,linewidth=width,label=bn)
    plt.xlim(0,max(traj.time))
    plt.xlabel("Time(fs)",fontweight='bold',size=15)
    plt.tick_params(axis='both',which='major',direction='in',length=5,top=True,bottom=True,left=True,right=True,width=2)
    plt.ylabel(r"Distance ($\mathregular{\AA}$)",fontweight='bold',size=15)
    plt.legend()
    plt.savefig('distance'+fig_ext+format,dpi=400.00,bbox_inches = 'tight',pad_inches = 0)
    plt.show()
def mu(traj,atoms,fig_ext='',format='png'):
    for atom in atoms:
        magmom=np.asarray([frame.atoms[atom].mu for frame in traj.frames])
        label='$\mathregular{\mu_{'+traj.frames[0].atoms[atom].ele+'-'+str(traj.frames[0].atoms[atom].ind)+'} }$'
        plt.plot(traj.time,magmom,linewidth = 2.5, label=label)
    plt.xlim(0,max(traj.time))
    plt.xlabel("Time (fs)",fontweight='bold',size=15)
    plt.tick_params(axis='both',which='major',direction='in',length=5,top=True,bottom=True,left=True,right=True,width=2)
    plt.ylabel(r"Mag. Moment",fontweight='bold',size=15)
    plt.legend()
    plt.savefig('mu'+fig_ext+format,dpi=400.00,bbox_inches = 'tight',pad_inches = 0)
    plt.show()