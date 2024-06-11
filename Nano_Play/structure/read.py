import os
import numpy as np
from ..nanite import *
def xyz(xyzfile):
    """The xyz file will be read and array of Atom objects with coordinates will be returned"""
    xyz=open(xyzfile,'r')
    split=xyz.readline().split()
    tot=int(split[0])
    xyz.readline()
    atoms=[]
    for i in range(tot):
        split=xyz.readline().split()
        atoms.append(Atom(i,symbol=split[0]))
        atoms[i].update_position([split[1],split[2],split[3]])
    return atoms
def contcar(contcar):
    """The VASP format coordinate file will be read and the b-matrix(Lattice Parameters)
      and the atom coordinates(xyz and frac) will be saved in atom.xyz and atom.frac as numpy array 
      Return values atoms,b"""
    con=open(contcar,'r')
    for i in range(2):
        con.readline()
    cmpnts=[]
    for i in range (3):
        split=con.readline().split()
        for x in split:
            cmpnts.append(float(x))
    bv=np.asarray(cmpnts)
    b=bv.reshape(3,3)
    split=con.readline().split()
    atomname=[str(x) for x in split]
    split=con.readline().split()
    natoms=[int(x) for x in split]
    atoms,index=[],0
    for i in range(len(natoms)):
        for j in range(natoms[i]):
            atoms.append(Atom(index,symbol=atomname[i]))
            index=index+1
    line=con.readline()
    mob=False
    if 'S' in line:
        line=con.readline()
        mob=True
    if 'D' in line or 'F' in line or 'd' in line or 'f' in line:
        coord=[]
        for i in range(len(atoms)):
            split=con.readline().split()
            for j in range(3):
                coord.append(float(split[j]))
                atoms[i].frac[j]=float(split[j])
                if mob==True:
                    atoms[i].as_mobile([split[-3],split[-2],split[-1]])
        tmp=np.asarray(coord)
        a=tmp.reshape(len(atoms),3)
        con.close()
        xyz=np.transpose(np.matmul(np.transpose(b),np.transpose(a)))
        for i in range(len(atoms)):
            atoms[i].update_position(xyz[i])

    else:
        xyz=[]
        for atom in atoms:
            split=con.readline().split()
            atom.update_position([split[0],split[1],split[2]])
            xyz.append([float(split[0]),float(split[1]),float(split[2])])
            if mob==True:
                atom.as_mobile([split[-3],split[-2],split[-1]])
        b_t_inv=np.linalg.inv(np.transpose(b))
        a=np.transpose(np.matmul(b_t_inv,np.transpose(xyz)))
        for i in range(len(atoms)):
            atoms[i].frac=np.asarray(a[i])         
    con.close()
    return atoms,b
#####Trajectory Formats ##########################
def xyz_traj(xyz_traj):
    """ The trajectory files written in xyz format will be read and array of Frame objects will be returned. """
    xyz=open(xyz_traj,'r')
    line=xyz.readline()
    tot=int(line.split()[0])
    frames,nframe=[],0
    while line.split():
        xyz.readline()
        frames.append(Frame(nframe))
        atoms=[]
        for i in range(tot):
            split=xyz.readline().split()
            atoms.append(Atom(i,symbol=split[0]))
            atoms[-1].update_postion([split[1],split[2],split[3]])
        frames[-1].update_atoms(atoms)
        line=xyz.readline()
        nframe=nframe+1
    return frames
def xdatcar(path='./',file='XDATCAR'):
    xdat=open(os.path.join(path,file),'r')
    trajectory=Trajectory()
    nframe=0
    for i in range(2):
        xdat.readline()
    cmpnts=[]
    for i in range(3):
        split=xdat.readline().split()
        for x in split:
            cmpnts.append(float(x))
    bv=np.asarray(cmpnts)
    b=bv.reshape(3,3)
    split=xdat.readline().split()
    atomname=[str(x) for x in split]
    split=xdat.readline().split()
    natoms=[int(x) for x in split]
    trajectory.elements=atomname
    trajectory.ele_nums=natoms
    frames=[]
    line=xdat.readline()
    nframe=0
    for line in xdat:
        if not line.strip():
            break
        atoms,ind=[],0
        for i in range(len(natoms)):
            for j in range(natoms[i]):
                atoms.append(Atom(ind,atomname[i]))
                ind=ind+1
        nframe=nframe+1
        coord=[]
        frame=Frame(nframe)
        for i in range(len(atoms)):
            split=line.split()
            for j in range(3):
                coord.append(float(split[j]))
            line=xdat.readline()
        tmp=np.asarray(coord)
        a=tmp.reshape(len(atoms),3)
        #line=xdat.readline()
        xyz=np.transpose(np.matmul(np.transpose(b),np.transpose(a)))
        for i in range(len(atoms)):
            atoms[i].frac=a[i]
            atoms[i].update_position([xyz[i][0],xyz[i][1],xyz[i][2]])
        frame.atoms=atoms
        frames.append(frame)
    trajectory.b_mat=b
    trajectory.frames=frames
    return trajectory