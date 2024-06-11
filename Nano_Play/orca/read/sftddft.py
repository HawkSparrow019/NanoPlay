from Nano_Play.nanite import *
from ..tools import *

def sftddft(out_file,out_path='./',nanite=None):
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
            if not line.strip():
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
    
    

    mos=[nanite.a_mos,nanite.b_mos]
    while "SF-TDA EXCITED STATES" not in fi.readline():
        continue
    while "SPIN-FLIP GROUND STATE" not in fi.readline():
        continue
    line=fi.readline()
    states=[]
    while line.strip():
        split=line.split()
        states.append(State(int(split[1][:-1]),float(split[5])))
        states[-1].s2=float(split[-1])
        line=fi.readline()
        while line.strip():
            update_sftddft_transition(states[-1],line,mos)
            line=fi.readline()
        line=fi.readline()
    nanite.states=states
    return nanite
    
    