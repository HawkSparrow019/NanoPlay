from ..nanite import *
from .tools import *

def dft(out_file):
    fi=open(out_file,'r')
    spin=True
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
    for i in range(11):
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
            elif "TIMINGS" in line:
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
    fi.close()
    return mol
def stddft(out_file,out_path='./',nanite=None):
    if nanite == None:
        nanite=Nanite(out_path)
    fi=open(out_file,'r')
    spin=True
    while "ORBITAL ENERGIES" not in fi.readline():
        continue
    fi.readline()
    if not fi.readline().strip():
        nanite.spin=False
    nanite.as_molecule()
    ####---------This block will read the MOs------------------------
    if nanite.spin==True:
        fi.readline()
        for line in fi:
            if not line.strip():
                break
            else:
                split=line.split()
                nanite.a_mos.append(MO(int(split[0]),float(split[1]),float(split[3])))
        for i in range(2):
            fi.readline()
        for line in fi:
            if "------" in line:
                break
            else:
                split=line.split()
                nanite.b_mos.append(MO(int(split[0]),float(split[1]),float(split[3])))
        ahomo=find_fmo(nanite.a_mos)
        bhomo=find_fmo(nanite.b_mos)
    else:
        fi.readline()
        for line in fi:
            if "----" in line:
                break
            else:
                split=line.split()
                nanite.mos.append(MO(int(split[0]),float(split[1]),float(split[3])))
        homo=find_fmo(nanite.mos)
        print(homo)
    ##--This block will identify the ordered HOMO-----
    while "ordered frontier orbitals" not in fi.readline():
        continue
    fi.readline()
    if nanite.spin==True:
        line=fi.readline()
        while line.strip():
            a_ofo=int(line.split()[0])
            line=fi.readline()
        scale=ahomo-a_ofo
    else:
        line=fi.readline()
        while line.strip():
            ofo=int(line.split()[0])
            line=fi.readline()
        scale=homo-ofo
    ###----This block will be update the transitions----------------------
    while "excitation energies, transition moments and amplitudes" not in fi.readline():
        continue
    while "state" not in fi.readline():
        continue
    states=[]
    line=fi.readline()
    while line.strip():
        states.append(update_TDDFT_state(line,nanite.spin,scale))
        line=fi.readline()
    while "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" not in fi.readline():
        continue
    for i in range(4):
        fi.readline()
    for state in states:
        state.intsty=float(fi.readline().split()[4])
    nanite.update_uvvis(states)
    fi.close()
    return nanite
def casscf(out_file,out_path='./',nanite=None):
    if nanite == None:
        nanite=Nanite(out_path)
    fi=open(out_file,'r')
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
    ###---------Updatation of the enrgetics data--------------------------
    while "CASSCF RESULTS" not in fi.readline():
        continue
    