import os
from tkinter.messagebox import YES

from numpy import asarray
from Nano_Play.nanite import *
from ..tools import *

def tddft(out_file,out_path='./',nanite=None):
    if nanite == None:
        nanite=Nanite(out_path)
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
    
            
    ### Initializing the MOs ############        
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
                nanite.a_mos.append(MO(int(split[0]),float(split[1]),float(split[2])))
        for i in range(2):
            fi.readline()
        for line in fi:
            if not line.strip():
                break
            else:
                split=line.split()
                nanite.b_mos.append(MO(int(split[0]),float(split[1]),float(split[2])))
        ahomo=find_fmo(nanite.a_mos)
        bhomo=find_fmo(nanite.b_mos)
        nanite.mos=[nanite.a_mos,nanite.b_mos]
    else:
        fi.readline()
        for line in fi:
            if not line.strip():
                break
            else:
                split=line.split()
                nanite.mos.append(MO(int(split[0]),float(split[1]),float(split[2])))
        homo=find_fmo(nanite.mos)
    
    ###---TD-DFT ---
    abs,sf,uv_vis,xas=False,False,False,False
    for line in fi:
        if "ORCA TD-DFT/TDA CALCULATION" in line:
            break
        elif "ORCA TERMINATED NORMALLY" in line:
            print("Not a TD-DFT output")
            exit
    ###----For getting data of the TD-DFT transitions------------------
  
    for line in fi:
        if "TD-DFT/TDA EXCITED STATES"  in line:
            abs=True
            break
        elif "SF-TDA EXCITED STATES" in line:
            sf=True
            break
    
    ##--Spin Flip TDDFT ---
        
    if sf==True:
        while "SPIN-FLIP GROUND STATE" not in fi.readline():
            continue
        line=fi.readline()
        states=[]
        while line.strip():
            split=line.split()
            states.append(State(int(split[1][:-1]),float(split[3])))
            states[-1].s2=float(split[-1])
            line=fi.readline()
            while line.strip():
                update_sftddft_transition(states[-1],line,nanite.mos)
                line=fi.readline()
            line=fi.readline()
        nanite.states=states

    ### Normal TDDFT spectra        
    if abs == True:
        for i in range(5):
            line=fi.readline()
        states=[]
        while line.strip():
            split=line.split()
            states.append(State(int(split[1][:-1]),float(split[3])))
            if '<S**2>' in line:
                states[-1].s2=float(split[-1])
            line=fi.readline()
            while line.strip():
                update_tddft_transition(states[-1],line,nanite.mos)
                line=fi.readline()
            line=fi.readline()
        nanite.states=states
    
        if nanite.states[0].lam < 1:
            uv_vis = True
        elif nanite.states[0].lam > 10:
            xas= True
    if uv_vis == True:
        while "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS" not in fi.readline():
            continue
        for i in range(4):
            fi.readline()
        for state in nanite.states:
            split=fi.readline().split()
            state.lam=float(split[2])
            state.intsty=float(split[3])
        mint=max([x.intsty for x in nanite.states])
        lam,intsty=[],[]
        for state in nanite.states:
            lam.append(state.lam)
            intsty.append(state.intsty)
            state.norm_intsty=state.intsty/mint
        nanite.uvvis=Spectrum(lam,intsty,'UVVIS')
        fi.close()

        



    print(abs,sf,uv_vis,xas)
    return nanite




def XAS(outfile):
    """"""
    fi=open(outfile,'r')
    values={}
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
    while "ENERGY FILE WAS UPDATED" not in fi.readline():
        continue
    fi.readline()
    spin=True
    if "SPIN" not in fi.readline():
        spin=False
    #----Creation of Molecule object with the spin information----------------
    mol=Nanite(name=outfile.split('/')[-1].split('.')[0],spin=spin)
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
    
     ###-----------------For Detecting the reduced orbital per MO-------------------------------
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

    if mol.spin == True:
        mol.a_mos,mol.b_mos=amos,bmos
        mos=[mol.a_mos,mol.b_mos]
    else:
        mol.mos=mos  
    ###----For getting data of the TD-DFT transitions------------------
    while "TD-DFT/TDA EXCITED STATES" not in fi.readline():
        continue
    for i in range(5):
        line=fi.readline()
    states=[]
    while line.strip():
        split=line.split()
        states.append(State(int(split[1][:-1]),float(split[5])))
        if '<S**2>' in line:
            states[-1].s2=float(split[-1])
        line=fi.readline()
        while line.strip():
            update_tddft_transition(states[-1],line,mos)
            line=fi.readline()
        line=fi.readline()
    mol.states=states
    freq,intsty=[7108],[0.0]
    for line in fi:
        if "(Exact Formulation, Velocity)" in line or "COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM (Origin Independent, Length Representation)" in line:
            break
    for i in range(4):
        fi.readline()
    for state  in mol.states:
        line=fi.readline()
        if not line.strip():
            break
        else:
            split=line.split()
            if len(split) > 5:
                freq.append((float(split[1])/8065.544718579)+14.208)
                intsty.append(float(split[6]))
                state.intsty=float(split[6])
    freq.append(7150.0)
    intsty.append(0)
    mol.xas=Spectrum(freq,intsty,'XAS')
    return mol



def IR(outfile):
    """This will read the vibrational frequencies and return the Nanite with ir Spectrum"""
    fi=open(outfile,'r')
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
    ###Creation of nanite object###
    mol=Nanite(name=outfile.split('/')[-1].split('.')[0])
    while "VIBRATIONAL FREQUENCIES" not in fi.readline():
        continue    
    for i in range(4):
        fi.readline()
    imag=False
    for line in fi:
        if not line.strip():
            break
        else:
            split=line.split()
            if float(split[1])<0.0:
                print("imaginary mode found at \n", float(split[1]),"  cm-1\n\n")
                imag=True
    if imag == False:            
        print("No imaginary frequency has been found. Cheers!\n\n")
    freq,intsty=[0.0],[0.0]
    while "IR SPECTRUM" not in line:
        line=fi.readline()
    while "---------------------" not in fi.readline():
        continue
    for line in fi:
        if not line.strip():
            break
        else:
            split=line.split()
            freq.append(float(split[1]))
            intsty.append(float(split[2]))
    freq.append(4500)
    intsty.append(0.0)
    ir=Spectrum(freq,intsty,type='ir')
    mol.spectra(ir)
    return mol