from Nano_Play.nanite import Atom
from . import read 
from . import write
from .tools import *
import numpy 
#import matplotlib.pyplot as plt



def add_atoms(add_atoms):
    """Adds list of atoms and return single list with all Atom """
    new_atoms=[]
    for atoms in add_atoms:
        for atom in atoms:
            new_atoms.append(atom)
    return new_atoms

def add_xyz(xyzs,signature=False,xyzfile='final.xyz'):
    """Adds multiple xyz files and write """    
    atoms_list,new_atoms=[],[]
    for xyz in xyzs:
        atoms_list.append(read.xyz(xyz))
    for atoms in atoms_list:
        for atom in atoms:
            new_atoms.append(atom)
    if signature == False:
        write.xyz(new_atoms,xyzfile=xyzfile)
    else:
        out=open(xyzfile,'w')
        print(len(new_atoms),"\n",file=out)
        for i in range(len(atoms_list)):
            for atom in atoms_list[i]:
                print("%2s(%i)   %15.9f  %15.9f  %15.9f"%(atom.ele,i+1,atom.xyz[0],atom.xyz[1],atom.xyz[2]),file=out)
        out.close()
def allign(atoms,a_axis,r_axis,ai=None,aj=None):
    """Allign Atoms with the cartesian axes and changes the origin"""
    if ai == None and aj == None:
        ai=int(input("Enter the atom index (i=0) of the center of rotation..\t"))
        aj=int(input("Eneter the atom index (i=0) of the rotating atom....\t"))
    elif ai == None:
        ai=int(input("Enter the atom index (i=0) of the center of rotation..\t"))
    elif aj == None:
        aj=int(input("Eneter the atom index (i=0) of the rotating atom....\t"))
    v1=atoms[aj].xyz-atoms[ai].xyz
    print(v1)
    #if r_axis == 'x':
    #    v2=np.asarray([1.0,0.0,0.0])
    #elif r_axis == 'y':
    #    v2=np.asarray([0.0,1.0,0.0])
    #elif r_axis == 'z':
    #    v2=np.asarray([0.0,0.0,1.0])
    if a_axis == 'x':
        v2=np.asarray([1.0,0.0,0.0])
    elif a_axis == 'y':
        v2=np.asarray([0.0,1.0,0.0])
    elif a_axis == 'z':
        v2=np.asarray([0.0,0.0,1.0])
    theta=angle(v1,v2,unit='rad')
    print((theta*180)/math.pi)
    d=input("Enter the direction of rotation; clockwise(c) or anti-clocwise(ac)")
    

    ct=numpy.cos(theta)
    if d=="c": 
        st=-numpy.sin(theta)
    else:
        st=numpy.sin(theta)

    #plt.quiver(,)


    if r_axis=='z':
        r_mat=[[ct,-st,0],[st,ct,0],[0,0,1]]
    elif r_axis=='y':
        r_mat=[[ct,0,-st],[0,1,0],[st,0,ct]]
    elif r_axis=='x':
        r_mat=[[1,0,0],[0,ct,-st],[0,st,ct]]
    #r_mat=[[ct,-st,0],[st,ct,0],[0,0,1]]
    a_mat=[a.xyz-atoms[ai].xyz for a in atoms]
    ad_mat=np.transpose(np.matmul(r_mat,np.transpose(a_mat)))
    new_atoms=[]
    for i in range(len(atoms)):
        new_atoms.append(Atom(i,atoms[i].ele))
        new_atoms[i].xyz=ad_mat[i]
    return new_atoms

