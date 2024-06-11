from Nano_Play.nanite import *
from ..tools import *

def casscf(out_file,out_path='./',nanite=None,orca_version = 5):
    if nanite == None:
        nanite=Nanite(name=out_file.split('/')[-1].split('.')[0])
    fi=open(out_file,'r')
    nanite.as_molecule(mr=True)
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
    
    ## Determining the active space ---

    while "ORCA-CASSCF" not in fi.readline():
        continue
    while "Determined orbital ranges" not in fi.readline():
        continue
    for i in range(2):
        line=fi.readline()
    acti,actf=int(line.split()[1]),int(line.split()[3])


    ###---------Updatation of the enrgetics data--------------------------
    while "CASSCF RESULTS" not in fi.readline():
        continue
    for i in range(3):
        line=fi.readline()
    nanite.tot_en=float(line.split()[-2])
    for i in range(6):
        fi.readline()
    mos=[]
    for line in fi:
        if not line.strip():
            break
        else:
            split=line.split()
            mos.append(MO(int(split[0]),float(split[1]),float(split[2])))
    ### Initialization of the CASSCF roots ###
    while "--------------------" not in fi.readline():
        continue
    #line=fi.readline()
    for line in fi:
        if "SA-CASSCF TRANSITION ENERGIES" in line:
            break
        else:
            multi=int(line.split()[6])
            block=[]
            for i in range(3):
                line=fi.readline()
            while line.strip():
                block.append(line)
                line=fi.readline()
        nanite.multis.update({multi:update_roots(block,multi)})
        for i in range(2):
            fi.readline()    
    for multi in nanite.multis:
        roots=nanite.multis.get(multi)
        for root in roots:
            nanite.roots.append(root)
    nanite.roots.sort(key=lambda x:x.en)
    
    ####  Evaluating the type of output ###
    reducd,nevpt2,soc=False,False,False    
    for line in fi:
        if "LOEWDIN ORBITAL-COMPOSITIONS"  in line:
            reducd=True
            break
        elif "Mayer bond orders" in line:
            mbo=True
            break
        elif "ORCA TERMINATED NORMALLY" in line:
            break
    ### Updating MO compositions ####
    if reducd == True:
        for i in range(3):
            line=fi.readline()
        while '----------------' not in line:
            mo_block=[]
            while line.strip():
                mo_block.append(line)
                line=fi.readline()
            if len(mo_block)>3:
                update_composition(mos,mo_block)
            line=fi.readline()    
    if nanite.spin==True:
        while "MULLIKEN ATOMIC CHARGES AND SPIN DENSITIES" not in fi.readline():
            continue
    else:
        while "MULLIKEN ATOMIC CHARGES" not in fi.readline():
            continue
    fi.readline()
    for i in range(len(atoms)):
        split=fi.readline().split(':')
        if nanite.spin==True:
            atoms[i].update_elec_prop()    
            atoms[i].mlkn_chg=float(split[-1].split()[-2])
            atoms[i].mlkn_spin=float(split[-1].split()[-1])
        else:
            atoms[i].mlkn_chg=float(split[-1])
    ####-------Total Population
    nanite.mlkn_chg=float(fi.readline().split()[-1])
    if nanite.spin == True:
        nanite.mlkn_spin=float(fi.readline().split()[5])
    mlkn_chg_block=[]
    for i in range(5):
        fi.readline()
    for line in fi:
        if not line.strip():
            break
        else:
            mlkn_chg_block.append(line)
    update_ml_resolved_population(atoms,mlkn_chg_block,'mlkn_chg')
    if nanite.spin == True:
        mlkn_spin_block=[]
        fi.readline()
        for line in fi:
            if not line.strip():
                break
            else:
                mlkn_spin_block.append(line)
        update_ml_resolved_population(atoms,mlkn_spin_block,'mlkn_spin')
    
    
        ######------------------------------LOEWDIN Population------------------- 
    if nanite.spin==True:
        while "LOEWDIN ATOMIC CHARGES AND SPIN DENSITIES" not in fi.readline():
            continue
    else:
        while "LOEWDIN ATOMIC CHARGES" not in fi.readline():
            continue
    fi.readline()
    for i in range(len(atoms)):
        split=fi.readline().split(":")
        if nanite.spin==True:
            atoms[i].lwdn_chg=float(split[-1].split()[-2])
            atoms[i].lwdn_spin=float(split[-1].split()[-1])
        else:
            atoms[i].lwdn_chg=float(split[-1])
    tlwdnchg=0.0
    tlwdnspin=0.0
    for atom in atoms:
        tlwdnchg=tlwdnchg+atom.lwdn_chg
        if nanite.spin==True:
            tlwdnspin=tlwdnspin+atom.lwdn_spin
    nanite.lwdn_chg=tlwdnchg
    if nanite.spin == True:
        nanite.lwdn_spin=tlwdnspin
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
    if nanite.spin==True:
        lwdn_spin_block=[]
        line=fi.readline()
        while line.strip():
            lwdn_spin_block.append(line)
            line=fi.readline()
        update_ml_resolved_population(atoms,lwdn_spin_block,'lwdn_spin')  
    #### MBO anlysis###
    nanite.atoms=atoms
    while "Mayer bond orders" not in fi.readline():
        continue
    block,line=[],fi.readline()
    while line.strip():
        block.append(line)
        line=fi.readline()
    update_mbo_data(block,nanite)
###    

    for line in fi:
        if "< NEVPT2  >" in line:
            nevpt2=True
            break
        elif "CASSCF RELATIVISTIC PROPERTIES" in line :
            soc=True
            break
        elif "ORCA TERMINATED NORMALLY" in line:
            break
    if nevpt2 == True:
        while "NEVPT2 TOTAL ENERGIES" not in fi.readline():
            continue
        while "STATE   ROOT MULT" not in fi.readline():
            continue
        for multi in nanite.multis:
            for root in nanite.multis.get(multi):
                energy=fi.readline().split()[3]
                root.cas_en=root.en
                root.en=float(energy)
        for line in fi:
            if "CASSCF RELATIVISTIC PROPERTIES" in line :
                soc=True
                break
            elif "ORCA TERMINATED NORMALLY" in line:
                break
    nanite.cas_roots=nanite.roots
    nanite.roots=[]
    for multi in nanite.multis:
        roots=nanite.multis.get(multi)
        for root in roots:
            nanite.roots.append(root)
    nanite.roots.sort(key=lambda x:x.en)
    nanite.soc,nanite.gtensor=False,False
    if soc == True:
        nanite.soc=True
        gtensor=False
        if orca_version == 4:
            while  "EFFECTIVE HAMILTONIAN SPIN-ORBIT COUPLING CONTRIBUTION" not in fi.readline():
                continue
        else:
            while "EFFECTIVE HAMILTONIAN SOC CONTRIBUTION" not in fi.readline():
                continue
        while "Eigenvectors:" not in fi.readline():
            continue
        dT=[]
        for i in range(3):
            split=fi.readline().split()
            dT.append([float(x) for x in split])
        nanite.dT=np.asarray(dT)
        for line in fi:
            if "D   =" in line:
                nanite.D = float(line.split()[-2])
                break
        nanite.EbD=float(fi.readline().split()[-1])
        for i in range(3):
            fi.readline()
        for line in fi:
            if line.strip():
                split=line.split()
            else:
                break
            for root in nanite.roots:
                if root.mlt == int(split[1]) and root.ind == int(split[2]):
                    root.D=float(split[3])
                    root.E=float(split[4])        
        for line in fi:
            if "ELECTRONIC G-MATRIX FROM EFFECTIVE HAMILTONIAN" in line:
                gtensor=True
                break
            elif "ORCA TERMINATED NORMALLY" in line:
                break
        if gtensor== True:
            nanite.gtensor=True
            gT=[]
            while "g-factors:" not in fi.readline():
                continue
            nanite.g = np.asarray([float(x) for x in fi.readline().split()[:3]])
            while "Orientation:" not in fi.readline():
                continue
            for i in range(3):
                gT.append(np.asarray([float(x) for x in fi.readline().split()[1:]]))
            nanite.gT=np.asarray(gT)
        
        
        if nevpt2 ==True:
            while "QDPT WITH NEVPT2 DIAGONAL ENERGIES" not in fi.readline():
                continue
            while "EFFECTIVE HAMILTONIAN SPIN-ORBIT COUPLING CONTRIBUTION" not in fi.readline():
                continue
            while "Eigenvectors:" not in fi.readline():
                continue            
            dT=[]
            for i in range(3):
                split=fi.readline().split()
                dT.append([float(x) for x in split])
            nanite.dT=np.asarray(dT)
            for line in fi:
                if "D   =" in line:
                    nanite.D = float(line.split()[-2])
                    break
            nanite.EbD=float(fi.readline().split()[-1])
            for i in range(3):
                fi.readline()
            nanite.roots[0].D=0.0
            nanite.roots[0].E=0.0
            for line in fi:
                if line.strip():
                    split=line.split()
                else:
                    break
                for root in nanite.roots:
                    if root.mlt == int(split[1]) and root.ind == int(split[2]):
                            root.D=float(split[3])
                            root.E=float(split[4])

            for line in fi:
                if "ELECTRONIC G-MATRIX FROM EFFECTIVE HAMILTONIAN" in line:
                    gtensor=True
                    break
                elif "ORCA TERMINATED NORMALLY" in line:
                    break
            if gtensor== True:
                nanite.gtensor=True
                gT=[]
                while "g-factors:" not in fi.readline():
                    continue
                nanite.g = np.asarray([float(x) for x in fi.readline().split()[:3]])
                while "Orientation:" not in fi.readline():
                    continue
                for i in range(3):
                    gT.append(np.asarray([float(x) for x in fi.readline().split()[1:]]))
                nanite.gT=np.asarray(gT)
    print(reducd,nevpt2,soc)
    nanite.mos=mos
    nanite.act_mos=mos[acti:actf+1]
    #nanite.update_bonds()
    return nanite