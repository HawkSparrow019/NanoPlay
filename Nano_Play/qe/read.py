import os
import numpy as np
from ..nanite import *
def pdos(system):
    dos_path=os.path.join(system.path,'plot')
    tot_pdos=open(os.path.join(dos_path,[f for f in os.listdir(dos_path) if f.endswith('pdos_tot')][0]),'r')
    tot_pdos.readline()
    en=[]
    for line in tot_pdos:
        en.append(float(line.split()[0]))
    pdos_files=[f for f in os.listdir(dos_path) if '#' in f]
    atom_files={}
    for f in pdos_files:
        ind=f.split('.')[-1].split('#')[-2].split('(')[0]
        if ind in atom_files:
                atom_files.get(ind).append(f)
        else:
                atom_files.update({ind:[f]})
    atoms=[]
    for atom in atom_files:
        atoms.append(Atom(atom))
        for file in atom_files.get(atom):
            pdos=open(os.path.join(dos_path,file),'r')
            if 's' in file.split('#')[-1]:
                s_up,s_dwn=[],[]
                pdos.readline()
                for line in pdos:
                    s_up.append(float(line.split()[-2]))
                    s_dwn.append(float(line.split()[-1]))
                atoms[-1].orbs.update({'s':Orbital('s')})
                atoms[-1].orbs.get('s').update_dos([s_up,s_dwn])
            elif 'p' in file.split("#")[-1]:
                pdos.readline()
                px_up,px_dwn,py_up,py_dwn,pz_up,pz_dwn=[],[],[],[],[],[]
                for line in pdos:
                    split=line.split()
                    pz_up.append(float(split[3]))
                    pz_dwn.append(float(split[4]))
                    px_up.append(float(split[5]))                
                    px_dwn.append(float(split[6]))
                    py_up.append(float(split[7]))
                    py_dwn.append(float(split[8]))
                atoms[-1].orbs.update({'px':Orbital('px')})
                atoms[-1].orbs.get('px').update_dos([px_up,px_dwn])
                atoms[-1].orbs.update({'py':Orbital('py')})
                atoms[-1].orbs.get('py').update_dos([py_up,py_dwn])
                atoms[-1].orbs.update({'pz':Orbital('pz')})
                atoms[-1].orbs.get('pz').update_dos([pz_up,pz_dwn])
            elif 'd' in file.split("#")[-1]:
                pdos.readline()
                dxz_up,dxz_dwn,dyz_up,dyz_dwn,dz2_up,dz2_dwn,dx2y2_up,dx2y2_dwn,dxy_up,dxy_dwn=[],[],[],[],[],[],[],[],[],[]
                for line in pdos:
                    split=line.split()
                    dz2_up.append(float(split[3]))
                    dz2_dwn.append(float(split[4]))
                    dxz_up.append(float(split[5]))
                    dxz_dwn.append(float(split[6]))
                    dyz_up.append(float(split[7]))
                    dyz_dwn.append(float(split[8]))
                    dx2y2_up.append(float(split[9]))
                    dx2y2_dwn.append(float(split[10]))
                    dxy_up.append(float(split[11]))
                    dxy_dwn.append(float(split[12]))
                atoms[-1].orbs.update({'dxz':Orbital('dxz')})
                atoms[-1].orbs.get('dxz').update_dos([dxz_up,dyz_dwn])
                atoms[-1].orbs.update({'dyz':Orbital('dyz')})
                atoms[-1].orbs.get('dyz').update_dos([dyz_up,dyz_dwn])
                atoms[-1].orbs.update({'dz2':Orbital('dz2')})
                atoms[-1].orbs.get('dz2').update_dos([dz2_up,dz2_dwn])
                atoms[-1].orbs.update({'dx2y2':Orbital('dx2y2')})
                atoms[-1].orbs.get('dx2y2').update_dos([dx2y2_up,dx2y2_dwn])
                atoms[-1].orbs.update({'dxy':Orbital('dxy')})
                atoms[-1].orbs.get('dxy').update_dos([dxy_up,dxy_dwn])
    output=[f for f in os.listdir(system.path) if f.endswith('.out')][0]
    scf_out=open(os.path.join(system.path,output),'r')
    system.atoms=atoms
    ef=0.0
    for line in scf_out:
        if "Fermi energy" in line:
            ef=float(line.split()[-2])
        elif "JOB DONE" in line:
            break
    system.update_es(en,ef)
