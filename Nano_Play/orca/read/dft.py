import os
from Nano_Play.nanite import *
from ..tools import *

def dft(out_file):
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
    conv,noiter,spin=False,False,True
    while "ORCA SCF" not in fi.readline():
        continue
    for line in fi:
        if "TOTAL SCF ENERGY" in line:
            conv = True
            break
        elif "skipping the SCF" in line:
            noiter=True
            spin,tot_en=True,0.00 ## By default spin will be true for the noiter case!
            break
    if conv== True: 
        while 'Total Energy' not in line:
            line=fi.readline()
            continue
        split=line.split()
        tot_en=float(split[-4])            
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
    if spin==True and noiter == False:
        line=fi.readline()
        while "Expectation value of <S**2>" not in line:
            line=fi.readline()
        mol.s2=float(line.split()[5])
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
                amos.append(MO(int(split[0]),float(split[1]),float(split[2])))
        for i in range(2):
            fi.readline()
        for line in fi:
            if not line.strip():
                break
            else:
                split=line.split()
                bmos.append(MO(int(split[0]),float(split[1]),float(split[2])))
    else:
        mos=[]
        while "ORBITAL ENERGIES" not in fi.readline():
            continue
        for i in range(3):
            fi.readline()
        for line in fi:
            if not line.strip():
                break
            else:
                split=line.split()
                mos.append(MO(int(split[0]),float(split[1]),float(split[2])))    
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
    mol.atoms=atoms
    ####----Loewdin reduced MO block----------------
    mbo,hirshfeld,reduced_mo,freq,uno=False,False,False,False,False
    for line in fi:
        if "HIRSHFELD ANALYSIS" in line:
            hirshfeld=True
            break
        elif "LOEWDIN REDUCED ORBITAL POPULATIONS PER MO" in line:
            reduced_mo=True
            break
        elif "Mayer bond orders" in line:
            mbo=True
            break
        elif "VIBRATIONAL FREQUENCIES" in line:
            freq= True
            break
        elif "UHF Natural Orbitals" in line:
            uno=True
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
        else:
            
            for i in range(3):
                line=fi.readline()
            while '********' not in line:
                mo_block=[]
                while line.strip():
                    mo_block.append(line)
                    line=fi.readline()
                if len(mo_block)>3:
                    update_composition(mos,mo_block)
                line=fi.readline()

        for line in fi:
            if "HIRSHFELD ANALYSIS" in line: 
                hirshfeld=True
                break
            elif "Mayer bond orders" in line:
                mbo=True
                break
            elif "***UHF Natural Orbitals" in line:
                uno=True
                break
            elif "VIBRATIONAL FREQUENCIES" in line:
                freq=True
                break
            elif "ORCA TERMINATED NORMALLY" in line:
                break    
    #####---Mayer Bond order----------
    if mbo ==True:
        block,line=[],fi.readline()
        while line.strip():
            block.append(line)
            line=fi.readline()
        mol.block=block
        update_mbo_data(block,mol)
        for line in fi:
            if "***UHF Natural Orbitals" in line:
                uno=True
                break
            if "HIRSHFELD ANALYSIS" in line: 
                hirshfeld=True
                break
            elif "VIBRATIONAL FREQUENCIES" in line:
                freq=True
                break
            elif "ORCA TERMINATED NORMALLY" in line:
                break 
    ####----Hirshfeld population block----------------
    if hirshfeld==True:
        if spin==True:
            for i in range(6):
                fi.readline()
            for atom in mol.atoms:
                split=fi.readline().split()
                atom.hf_chg=float(split[-2])
                atom.hf_spin=float(split[-1])
        else:
            for i in range(6):
                fi.readline()
            for atom in atoms:
                split=fi.readline().split()
                atom.hf_chg=float(split[-2])
        for line in fi:
            if "***UHF Natural Orbitals" in line:
                uno=True
                break
            elif 'VIBRATIONAL FREQUENCIES' in line:
                freq=True
                break
            elif "ORCA TERMINATED NORMALLY" in line:
                break
    ####----UNO orbitals------
    if uno == True:
        ###UNSO##
        while "LOEWDIN REDUCED ORBITAL POPULATIONS PER UNSO" not in fi.readline():
            continue
        unso=[MO(x,0,0) for x in range(len(amos))]
        for i in range(3):
            line=fi.readline()
        while line.strip():
            mo_block=[]
            while line.strip():
                mo_block.append(line)
                line=fi.readline()
            if len(mo_block)>3:
                update_composition(unso,mo_block,en_occ=False)
            line=fi.readline()    
        mol.unso=unso
        ####UNO##
        while "LOEWDIN REDUCED ORBITAL POPULATIONS PER UNO" not in fi.readline():
            continue
        un=[MO(x,0,0) for x in range(len(amos))]
        for i in range(3):
            line=fi.readline()
        while "QR-MO GENERATION" not in line:
            mo_block=[]
            while line.strip():
                mo_block.append(line)
                line=fi.readline()
            if len(mo_block)>3:
                update_composition(un,mo_block,en_occ=False)
            line=fi.readline()    
        mol.uno=un
        ##QRO###
        while "Orbital Energies of Quasi-Restricted MO's" not in fi.readline():
            continue
        qros=[]
        for line in fi:
            if "-------" in line :
                break
            else:
                split=line.split()
                ind=split[0].split('(')[0]
                en=float(split[5])*0.037
                occ=split[1].split(')')[0]
                qros.append(MO(ind,occ,en))
        mol.qro=qros
        for line in fi:
            if 'VIBRATIONAL FREQUENCIES' in line:
                freq=True
                break
            elif "ORCA TERMINATED NORMALLY" in line:
                break



    #mol.atoms=atoms      
    if mol.spin == True:
        mol.a_mos,mol.b_mos=amos,bmos
    else:
        mol.mos=mos
    
    mol.freq=False
    if freq == True:
        mol.freq=True
        nus=[]
        for i in range(4):
            fi.readline()
        for line in fi:
            if line.strip():
                split=line.split()
                if float(split[1])>0.0:
                    nus.append(float(split[1]))
                elif float(split[1])<0.0:
                    print("WARNING: negative frequency found at {} cm-1".format(float(split[1])))
            else:
                break
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
        mol.nus=nus
        line=fi.readline()
        ##while "ENTHAPLY" not in line:
         #   print(line)
         #   line=fi.readline()  
            #continue 
        for line in fi:
            if 'Total Enthalpy' in line:
                mol.h=float(line.split()[-2])
                break
        #while "ENTROPY" not in fi.readline():
        #    continue 
        for line in fi:
            if 'Final entropy term' in line:
                mol.s=float(line.split()[-4])
                break
        #while "GIBBS FREE ENTHALPY" not in fi.readline():
        #    continue 
        for line in fi:
            if 'Final Gibbs free enthalpy ' in line:
                mol.g=float(line.split()[-2])
                break
    fi.close()
    #print(reduced_mo,mbo,uno,hirshfeld,freq)
    return mol