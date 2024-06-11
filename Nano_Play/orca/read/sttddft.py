from Nano_Play.nanite import *
from ..tools import *

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