from Nano_Play.nanite import *
from ..tools import *

def rohf(out_file,out_path='./',nanite=None):
    if nanite == None:
        nanite=Nanite(name=out_file.split('/')[-1].split('.')[0])
    fi=open(out_file,'r')
    nanite.as_molecule()
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
    while "ORBITAL ENERGIES" not in fi.readline():
        continue
    for i in range(3):
        line=fi.readline()
    #nanite.tot_en=float(line.split()[-2])
    #for i in range(6):
    #    fi.readline()
    mos=[]
    for line in fi:
        if not line.strip():
            break
        else:
            split=line.split()
            mos.append(MO(int(split[0]),float(split[1]),float(split[2])))
    ### Updating MO compositions ####    
    while "LOEWDIN REDUCED ORBITAL POPULATIONS PER MO" not in fi.readline():
        continue
    for i in range(3):
        line=fi.readline()
    while '**********' not in line:
        mo_block=[]
        while line.strip():
            mo_block.append(line)
            line=fi.readline()
        if len(mo_block)>3:
            update_composition(mos,mo_block)
        line=fi.readline()    

#    if nanite.spin==True:
#        while "MULLIKEN ATOMIC CHARGES AND SPIN DENSITIES" not in fi.readline():
#            continue
#    else:
#        while "MULLIKEN ATOMIC CHARGES" not in fi.readline():
#            continue
#    fi.readline()
#    for i in range(len(atoms)):
#        split=fi.readline().split(':')
#        if nanite.spin==True:
#            atoms[i].update_elec_prop()    
#            atoms[i].mlkn_chg=float(split[-1].split()[-2])
#            atoms[i].mlkn_spin=float(split[-1].split()[-1])
#        else:
#            atoms[i].mlkn_chg=float(split[-1])
#    ####-------Total Population
#    nanite.mlkn_chg=float(fi.readline().split()[-1])
#    if nanite.spin == True:
#        nanite.mlkn_spin=float(fi.readline().split()[5])
#    mlkn_chg_block=[]
#    for i in range(5):
#        fi.readline()
#    for line in fi:
#        if not line.strip():
#            break
#        else:
#            mlkn_chg_block.append(line)
#    update_ml_resolved_population(atoms,mlkn_chg_block,'mlkn_chg')
#    if nanite.spin == True:
#        mlkn_spin_block=[]
#        fi.readline()
#        for line in fi:
#            if not line.strip():
#                break
#            else:
#                mlkn_spin_block.append(line)
#        update_ml_resolved_population(atoms,mlkn_spin_block,'mlkn_spin')
#    
#    
    
    nanite.mos=mos
    nanite.atoms=atoms
    return nanite