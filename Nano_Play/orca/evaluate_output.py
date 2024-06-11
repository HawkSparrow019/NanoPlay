from Nano_Play.nanite import *
import os
from .read import *
from .tools import *
def read1(out_file, path = './'):
    fi=open(os.path.join(path,out_file),'r')
    spin=True
    calculation=None
    while 'END OF INPUT' not in fi.readline():
        continue
    for i in range(3):
        line=fi.readline()
    line=fi.readline()
    atoms=[]
    if "Single Point Calculation" in line:
        for i in range(5):
            fi.readline()
    elif "Geometry Optimization Run" in line:
        while "FINAL ENERGY EVALUATION AT THE STATIONARY POINT" not in fi.readline():
            continue
        for i in range(5):
            fi.readline()
    ###------Creation of atom objects from optimized geometry-----------------
    natom=0
    for line in fi:
        if not line.strip():
            break
        else:
            split=line.split()
            atoms.append(Atom(natom,split[0]))
            atoms[-1].xyz=np.asarray([float(split[1]),float(split[2]),float(split[3])])
            natom=natom+1
    mol=Nanite(name=out_file.split('/')[-1].split('.')[0],spin=spin)
    mol.atoms=atoms
    while "Hamiltonian:" not in fi.readline():
        continue
    print("Reached calculation type")
    line= fi.readline()
    if "DFT" in line:
        calculation='dft'
    if "Hartree-Fock" in line:
        calculation='hf'
    if calculation == 'dft':
        print("it is a dft calculation")
        mol=dft.dft(fi,mol)
    return mol