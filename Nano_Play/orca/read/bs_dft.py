import os
from Nano_Play.nanite import *
from ..tools import *
def bs_dft(out_file):
    fi=open(out_file,'r')
    spin=True
    while 'END OF INPUT' not in fi.readline():
        continue
    for i in range(3):
        line=fi.readline()
    line=fi.readline()
    atoms,bs_atoms=[],[]
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
    bs_atoms=atoms
    ###---------Updatation of the enrgetics data--------------------------
    while "TOTAL SCF ENERGY" not in fi.readline():
        continue
    line=fi.readline()
    while 'Total Energy' not in line:
        line=fi.readline()
        continue
    split=line.split()
    tot_en=float(split[5])            
    while 'SCF CONVERGENCE' not in fi.readline():
        continue
    #--------This block will determine whether it is a spin polirized calculation or not-----
    while "ENERGY FILE WAS UPDATED" not in fi.readline():
        continue
    fi.readline()
    if "SPIN" not in fi.readline():
        spin=False
    #----Creation of Molecule object with the spin information----------------
    mol=Nanite(name=out_file.split('/')[-1].split('.')[0],spin=spin)
    mol.tot_en=tot_en
    mol.as_molecule()
    if spin==True:
        for i in range(6):
            fi.readline()
        mol.s2=float(fi.readline().split()[5])
    ###----------------This block is for MOs' initialization----------------------------------
    if spin == True:
        amos,bmos=[],[]
        while "SPIN UP ORBITALS" not in fi.readline():
            continue
        fi.readline()
        for line in fi:
            if not line.strip():
                break
            else:
                split=line.split()
                amos.append(MO(int(split[0]),float(split[1]),float(split[3])))
        for i in range(2):
            fi.readline()
        for line in fi:
            if not line.strip():
                break
            else:
                split=line.split()
                bmos.append(MO(int(split[0]),float(split[1]),float(split[3])))
    else:
        mos=[]
        while "ORBITAL ENERGIES" not in fi.readline():
            fi.readline()
        for i in range(3):
            fi.readline()
        for line in fi:
            if not line.strip():
                break
            else:
                split=line.split()
                mos.append(MO(int(split[0]),float(split[1]),float(split[3])))    
    ###----------------This block is to update the population to atoms ----------------------------------
    ###------------------Mulliken Population---------------------------------------
    if spin==True:
        while "MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS" not in fi.readline():
            continue
    else:
        while "MULLIKEN ATOMIC CHARGES" not in fi.readline():
            continue
    fi.readline()
    for i in range(len(atoms)):
        split=fi.readline().split(':')
        if spin==True:
            atoms[i].update_elec_prop()    
            atoms[i].mlkn_chg=float(split[-1].split()[-2])
            atoms[i].mlkn_spin=float(split[-1].split()[-1])
        else:
            atoms[i].mlkn_chg=float(split[-1])
    ####-------Total Population
    mol.mlkn_chg=float(fi.readline().split()[-1])
    if spin == True:
        mol.mlkn_spin=float(fi.readline().split()[5])
    mlkn_chg_block=[]
    for i in range(5):
        fi.readline()
    for line in fi:
        if not line.strip():
            break
        else:
            mlkn_chg_block.append(line)
    update_ml_resolved_population(atoms,mlkn_chg_block,'mlkn_chg')
    if spin == True:
        mlkn_spin_block=[]
        fi.readline()
        for line in fi:
            if not line.strip():
                break
            else:
                mlkn_spin_block.append(line)
        update_ml_resolved_population(atoms,mlkn_spin_block,'mlkn_spin')
    ######------------------------------LOEWDIN Population------------------- 
    if spin==True:
        while "LOEWDIN ATOMIC CHARGES AND SPIN POPULATIONS" not in fi.readline():
            continue
    else:
        while "LOEWDIN ATOMIC CHARGES" not in fi.readline():
            continue
    fi.readline()
    for i in range(len(atoms)):
        split=fi.readline().split(":")
        if spin==True:
            atoms[i].lwdn_chg=float(split[-1].split()[-2])
            atoms[i].lwdn_spin=float(split[-1].split()[-1])
        else:
            atoms[i].lwdn_chg=float(split[-1])
    tlwdnchg=0.0
    tlwdnspin=0.0
    for atom in atoms:
        tlwdnchg=tlwdnchg+atom.lwdn_chg
        if spin==True:
            tlwdnspin=tlwdnspin+atom.lwdn_spin
    mol.lwdn_chg=tlwdnchg
    if spin == True:
        mol.lwdn_spin=tlwdnspin
    while "LOEWDIN REDUCED ORBITAL CHARGES" not in fi.readline():
        continue
    for i in range(2):
        fi.readline()
    lwdn_chg_block=[]
    line=fi.readline()
    while line.strip():
        lwdn_chg_block.append(line)
        line=fi.readline()
    update_ml_resolved_population(atoms,lwdn_chg_block,'lwdn_chg')
    fi.readline()
    if spin==True:
        lwdn_spin_block=[]
        line=fi.readline()
        while line.strip():
            lwdn_spin_block.append(line)
            line=fi.readline()
        update_ml_resolved_population(atoms,lwdn_spin_block,'lwdn_spin')  
     ####----Loewdin reduced MO block----------------
    hirshfeld,reduced_mo=False,False
    for line in fi:
        if "HIRSHFELD ANALYSIS" in line:
            hirshfeld=True
            break
        elif "LOEWDIN REDUCED ORBITAL POPULATIONS PER MO" in line:
            reduced_mo=True
            break
        elif "TIMINGS" in line:
            break
    if reduced_mo == True:
        if spin == True:
            while 'SPIN UP' not in fi.readline():
                continue
            line=fi.readline()
            while 'SPIN DOWN' not in line:
                mo_block=[]
                while line.strip():
                    mo_block.append(line)
                    line=fi.readline()
                if len(mo_block)>3:
                    update_composition(amos,mo_block)
                line=fi.readline()
            line=fi.readline()
            while '********' not in line:
                mo_block=[]
                while line.strip():
                    mo_block.append(line)
                    line=fi.readline()
                if len(mo_block)>3:
                    update_composition(bmos,mo_block)
                line=fi.readline()    
     ####----Hirshfeld population block----------------
    if reduced_mo == True:
        for line in fi:
            if "HIRSHFELD ANALYSIS" in line: 
                hirshfeld=True
                break
            elif "BROKEN SYMMETRY" in line:
                break
    if hirshfeld==True:
        if spin==True:
            for i in range(6):
                fi.readline()
            for atom in atoms:
                split=fi.readline().split()
                atom.hf_chg=float(split[-2])
                atom.hf_spin=float(split[-1])
        else:
            for i in range(6):
                fi.readline()
            for atom in atoms:
                split=fi.readline().split()
                atom.hf_chg=float(split[-2])
    mol.atoms=atoms      
    if mol.spin == True:
        mol.a_mos,mol.b_mos=amos,bmos
    else:
        mol.mos=mos

##-----------------------------------------Reading BS state data--------------
    
    bs_state=Nanite(name='BS',spin=spin)
    ###---------Updatation of the enrgetics data--------------------------
    for line in fi:
        if "TOTAL SCF ENERGY" in line:
            break
    line=fi.readline()
    while 'Total Energy' not in line:
        line=fi.readline()
        continue
    split=line.split()
    tot_en=float(split[5])            
    while 'SCF CONVERGENCE' not in fi.readline():
        continue
    #--------This block will determine whether it is a spin polirized calculation or not-----
    while "ENERGY FILE WAS UPDATED" not in fi.readline():
        continue
    fi.readline()
    if "SPIN" not in fi.readline():
        spin=False
    #----Creation of Molecule object with the spin information----------------
    bs_state.tot_en=tot_en
    bs_state.as_molecule()
    if spin==True:
        for i in range(6):
            fi.readline()
        bs_state.s2=float(fi.readline().split()[5])
    ###----------------This block is for MOs' initialization----------------------------------
    if spin == True:
        amos,bmos=[],[]
        while "SPIN UP ORBITALS" not in fi.readline():
            continue
        fi.readline()
        for line in fi:
            if not line.strip():
                break
            else:
                split=line.split()
                amos.append(MO(int(split[0]),float(split[1]),float(split[3])))
        for i in range(2):
            fi.readline()
        for line in fi:
            if not line.strip():
                break
            else:
                split=line.split()
                bmos.append(MO(int(split[0]),float(split[1]),float(split[3])))
    else:
        mos=[]
        while "ORBITAL ENERGIES" not in fi.readline():
            fi.readline()
        for i in range(3):
            fi.readline()
        for line in fi:
            if not line.strip():
                break
            else:
                split=line.split()
                mos.append(MO(int(split[0]),float(split[1]),float(split[3])))    
    ###----------------This block is to update the population to atoms ----------------------------------
    ###------------------Mulliken Population---------------------------------------
    if spin==True:
        while "MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS" not in fi.readline():
            continue
    else:
        while "MULLIKEN ATOMIC CHARGES" not in fi.readline():
            continue
    fi.readline()
    for i in range(len(bs_atoms)):
        split=fi.readline().split(':')
        if spin==True:
            bs_atoms[i].update_elec_prop()    
            bs_atoms[i].mlkn_chg=float(split[-1].split()[-2])
            bs_atoms[i].mlkn_spin=float(split[-1].split()[-1])
        else:
            bs_atoms[i].mlkn_chg=float(split[-1])
    ####-------Total Population
    bs_state.mlkn_chg=float(fi.readline().split()[-1])
    if spin == True:
        bs_state.mlkn_spin=float(fi.readline().split()[5])
    mlkn_chg_block=[]
    for i in range(5):
        fi.readline()
    for line in fi:
        if not line.strip():
            break
        else:
            mlkn_chg_block.append(line)
    update_ml_resolved_population(bs_atoms,mlkn_chg_block,'mlkn_chg')
    if spin == True:
        mlkn_spin_block=[]
        fi.readline()
        for line in fi:
            if not line.strip():
                break
            else:
                mlkn_spin_block.append(line)
        update_ml_resolved_population(bs_atoms,mlkn_spin_block,'mlkn_spin')
    ######------------------------------LOEWDIN Population------------------- 
    if spin==True:
        while "LOEWDIN ATOMIC CHARGES AND SPIN POPULATIONS" not in fi.readline():
            continue
    else:
        while "LOEWDIN ATOMIC CHARGES" not in fi.readline():
            continue
    fi.readline()
    for i in range(len(bs_atoms)):
        split=fi.readline().split(":")
        if spin==True:
            bs_atoms[i].lwdn_chg=float(split[-1].split()[-2])
            bs_atoms[i].lwdn_spin=float(split[-1].split()[-1])
        else:
            bs_atoms[i].lwdn_chg=float(split[-1])
    tlwdnchg=0.0
    tlwdnspin=0.0
    for atom in bs_atoms:
        tlwdnchg=tlwdnchg+atom.lwdn_chg
        if spin==True:
            tlwdnspin=tlwdnspin+atom.lwdn_spin
    bs_state.lwdn_chg=tlwdnchg
    if spin == True:
        bs_state.lwdn_spin=tlwdnspin
    while "LOEWDIN REDUCED ORBITAL CHARGES" not in fi.readline():
        continue
    for i in range(2):
        fi.readline()
    lwdn_chg_block=[]
    line=fi.readline()
    while line.strip():
        lwdn_chg_block.append(line)
        line=fi.readline()
    update_ml_resolved_population(bs_atoms,lwdn_chg_block,'lwdn_chg')
    fi.readline()
    if spin==True:
        lwdn_spin_block=[]
        line=fi.readline()
        while line.strip():
            lwdn_spin_block.append(line)
            line=fi.readline()
        update_ml_resolved_population(bs_atoms,lwdn_spin_block,'lwdn_spin')  
     ####----Loewdin reduced MO block----------------
    hirshfeld,reduced_mo=False,False
    for line in fi:
        if "HIRSHFELD ANALYSIS" in line:
            hirshfeld=True
            break
        elif "LOEWDIN REDUCED ORBITAL POPULATIONS PER MO" in line:
            reduced_mo=True
            break
        elif "TIMINGS" in line:
            break
    if reduced_mo == True:
        if spin == True:
            while 'SPIN UP' not in fi.readline():
                continue
            line=fi.readline()
            while 'SPIN DOWN' not in line:
                mo_block=[]
                while line.strip():
                    mo_block.append(line)
                    line=fi.readline()
                if len(mo_block)>3:
                    update_composition(amos,mo_block)
                line=fi.readline()
            line=fi.readline()
            while '********' not in line:
                mo_block=[]
                while line.strip():
                    mo_block.append(line)
                    line=fi.readline()
                if len(mo_block)>3:
                    update_composition(bmos,mo_block)
                line=fi.readline()    
     ####----Hirshfeld population block----------------
    if reduced_mo == True:
        for line in fi:
            if "HIRSHFELD ANALYSIS" in line: 
                hirshfeld=True
                break
            elif "TIMINGS" in line:
                break
    if hirshfeld==True:
        if spin==True:
            for i in range(6):
                fi.readline()
            for atom in bs_atoms:
                split=fi.readline().split()
                atom.hf_chg=float(split[-2])
                atom.hf_spin=float(split[-1])
        else:
            for i in range(6):
                fi.readline()
            for atom in bs_atoms:
                split=fi.readline().split()
                atom.hf_chg=float(split[-2])
    bs_state.atoms=bs_atoms      
    if mol.spin == True:
        bs_state.a_mos,bs_state.b_mos=amos,bmos
    else:
        bs_state.mos=mos
    return mol,bs_state

    