import os
import numpy as np
from   ..nanite import *
from  ..structure.read import contcar as read_contcar
from ..structure.tools import *
def pdos(path='./',nanite=None):
    '''This will read the VASP DOSCAR file and return a Nanite object with the orbital resolved DOS
    orbital resolved DOS are saved in a dictionary with {'orb':Orbital}
    '''
    ret=False
    if nanite == None:
        nanite=Nanite(path=path)
        nanite.as_nano()
        ret = True
    file=os.path.join(nanite.path,'DOSCAR')
    fi=open(file,'r')
    split=fi.readline().split()
    natoms=int(split[1])
    for i in range (4):
        fi.readline()
    split=fi.readline().strip().split()
    if len(split) == 4:
       dospt=int(split[1][12:])
       efermi=float(split[2])    
    else: 
        dospt=int(split[2])
        efermi=float(split[3])
    en,atoms,up_dos,dwn_dos,dos=[],[],[],[],[]
    for i in range (dospt):
        split=fi.readline().split()
        if i == 0:
            if len(split) == 3:
                nanite.spin = False
        en.append(float(split[0]))
        if nanite.spin == True:
            up_dos.append(float(split[1]))
            dwn_dos.append(float(split[2]))
        else:
            dos.append(float(split[1]))
    if len(nanite.atoms) == 0:
        for i in range(natoms):
            nanite.atoms.append(Atom(i+1))
    for i in range(natoms):
        fi.readline()
        dospts=[]    
        for j in range(dospt):
            split=fi.readline().split()
            dospts.append([float(split[k])for k in range(len(split))])
        lm_proj=[]
        for dpt in dospts:
            tmp=[]
            for x in dpt[1:]:
                tmp.append(float(x))
            lm_proj.append(tmp)
        lm_proj=np.transpose(np.asarray(lm_proj))
        count,orbs=0,['s','py','pz','px','dxy','dyz','dz2','dxz','dx2y2']
        #print(len(lm_proj))
        if len(lm_proj) == 18:
            for orb in orbs:
                nanite.atoms[i].orbs.update({orb:Orbital(orb)})
                orbital=nanite.atoms[i].orbs.get(orb)
                orbital.update_dos([lm_proj[count],lm_proj[count+1]])
                count=count+2
        elif len(lm_proj) == 36:
            nanite.ncl=True
            for orb in orbs:
                nanite.atoms[i].orbs.update({orb:Orbital(orb)})
                orbital=nanite.atoms[i].orbs.get(orb)
                orbital.update_dos(lm_proj[count:count+4])
                count=count+4
        else:
            nanite.spin=False
            for orb in orbs:
                nanite.atoms[i].orbs.update({orb:Orbital(orb)})
                orbital=nanite.atoms[i].orbs.get(orb)
                orbital.update_dos(lm_proj[count])
                count=count+1
    fi.close()
    if nanite.spin == True:
        nanite.update_es(en,efermi,[up_dos,dwn_dos])
    else:
        nanite.update_es(en,efermi,dos)
    for atom in nanite.atoms:
        orbs,tmp=atom.orbs,[]
        for o in orbs:
            orb=orbs.get(o)
            if nanite.spin == True and nanite.ncl == False:    
                tmp.append(max(orb.up_dos))
                tmp.append(max(orb.down_dos))
            else:
                tmp.append(max(orb.dos))
        nf=max(tmp)
        orbs.update({'p':orbs.get('px')+orbs.get('py')+orbs.get('pz')})    
        orbs.get('p').sym='$Total-p$'
        orbs.update({'d':orbs.get('dxy')+orbs.get('dxz')+orbs.get('dyz')+orbs.get('dz2')+orbs.get('dx2y2')})
        orbs.get('d').sym='$Total-d$'
        tot=orbs.get('s')+orbs.get('p')+orbs.get('d')
        if nanite.spin == True and nanite.ncl== False:
            nf=max([max(tot.up_dos),max(tot.down_dos)])
        else:
            nf=max(tot.dos)
        orbs.update({'total':tot})
        atom.nf=nf
        for o in orbs:
            orbs.get(o).occupancy(nanite.e_m_e_f)
            orbs.get(o).normalise(nf)
        #orbs.get('p').normalise(max([max(orbs.get('p').up_dos),max(orbs.get('p').down_dos)]))
        #orbs.get('d').normalise(max([max(orbs.get('d').up_dos),max(orbs.get('d').down_dos)]))
    if ret == True:
        return nanite
def scf(path='./',nanite=None,nsw=500,name=None):
    if nanite == None:
        ret=True
        if name != None:
            nanite=Nanite(path=path,name=name)
        else:
            nanite=Nanite(path=path)
    else:
        ret=False
    nanite.conv=True
    if "OSZICAR" in os.listdir(nanite.path):
        osz=open(os.path.join(nanite.path,'OSZICAR'),'r')
        final_step=int(osz.readlines()[-2].split()[1])
        if final_step >= nsw:
            nanite.conv=False
            #print("WARNING: It seems that the convergence is not achieved properly, check your output in\n",nanite.path)
    out=open(os.path.join(nanite.path,'OUTCAR'),'r')
    con=os.path.join(nanite.path,'CONTCAR')

    nanite.ncl=False
    while "Startparameter for this run:" not in out.readline():
        continue
    for line in out:
        if not line.strip():
            break
        elif "ISPIN" in line: 
            if '1' in line:
                nanite.spin=False
        elif "LNONCOLLINEAR =      T" in line:# and 'T' in line:
            nanite.ncl=True
   
    #print(ncl,nanite.spin)   
    nanite.as_nano(ncl=nanite.ncl)
    nanite.atoms,nanite.cell=read_contcar(con)
    for line in out:
        if "e  e" in line:
            nanite.tot_en=float(line.split()[4])
            break
    if nanite.ncl == True:
        count = 0
        while "Spin-Orbit-Coupling matrix elements" not in out.readline():
            continue
        for atom in nanite.atoms:
            #for line in out.readline():
            soc_mat=[]
            line=out.readline()
            while 'Ion' not in line:
                line=out.readline()
                continue
            atom.esoc=float(line.split()[-1])            
            for i in range(3):
                out.readline()
                l_mat=[]
                for j in range(2*(i+1)+1):
                    l_mat.append(np.asarray([float(x) for x in out.readline().split()]))
                soc_mat.append(np.asarray(l_mat))
            atom.soc_mat=soc_mat
    for atom in nanite.atoms:
        atom.update_elec_prop()
    if nanite.spin==True:
        while 'magnetization (x)' not in out.readline():
            continue
        for i in range(3):
            out.readline()
        for atom in nanite.atoms:
            split=out.readline().split()
            atom.update_mu([split[1],split[2],split[3],split[4]])
    elif nanite.ncl == True:
        mu=[[[],[],[],[]] for i in range(len(nanite.atoms))]
        while 'magnetization (x)' not in out.readline():
            continue
        for i in range(3):
            out.readline()
        for i in range(len(nanite.atoms)):
            split=out.readline().split()
            mu[i][0].append(float(split[1]))
            mu[i][1].append(float(split[2]))
            mu[i][2].append(float(split[3]))
            mu[i][3].append(float(split[4]))
        out.readline()
        nanite.mu[0]=float(out.readline().split()[-1])
        while 'magnetization (y)' not in out.readline():
            continue
        for i in range(3):
            out.readline()

        for i in range(len(nanite.atoms)):
            split=out.readline().split()
            mu[i][0].append(float(split[1]))
            mu[i][1].append(float(split[2]))
            mu[i][2].append(float(split[3]))
            mu[i][3].append(float(split[4]))
        out.readline()
        nanite.mu[1]=float(out.readline().split()[-1])
        while 'magnetization (z)' not in out.readline():
            continue
        for i in range(3):
            out.readline()
        for i in range(len(nanite.atoms)):
            split=out.readline().split()
            mu[i][0].append(float(split[1]))
            mu[i][1].append(float(split[2]))
            mu[i][2].append(float(split[3]))
            mu[i][3].append(float(split[4]))
        out.readline()
        nanite.mu[2]=float(out.readline().split()[-1])

        
        #mu=np.transpose(mu)
        for i,atom in enumerate(nanite.atoms):
            atom.update_ncl_magnetization(mu[i])
    if "ACF.dat" in os.listdir(nanite.path):
        acf=open(os.path.join(nanite.path,'ACF.dat'),'r')
        nelect=0.0
        for i in range(2):
            acf.readline()
        for atom in nanite.atoms:
            atom.bader_pop=float(acf.readline().split()[4])
            nelect=nelect+atom.bader_pop
        nanite.nelect=nelect
    con=open(os.path.join(nanite.path,'CONTCAR'),'r')
    for i in range(5):
        con.readline()
    element=[str(x) for x in con.readline().split()]
    ele_num=[int(x) for x in con.readline().split()]
    ions,index=[],0
    for i in range(len(ele_num)):
        ions.append(Ion(element[i]))
        for j in range(ele_num[i]):
            ions[i].atoms.append(nanite.atoms[index])
            index=index+1
    for ion in ions:
        if nanite.spin==True:
            ion.mu()
        ion.bader_population()
    nanite.ions=ions
    nanite.bader_population()
    find_layers(nanite)
    out.close()
    if ret == True:
        return nanite
def oszicar(traj,path='./'):
    osz=open(os.path.join(path,'OSZICAR'),'r')
    for line in osz:
        if "T=" in line:
            split=line.split()
            traj.frames[int(split[0])-1].temp=float(split[2])
            traj.frames[int(split[0])-1].en=float(split[4])
            traj.frames[int(split[0])-1].free_en=float(split[6])
            traj.frames[int(split[0])-1].magmom=float(split[-1])
    traj.update()      
def outcar_md(traj,path='./'):
    out=open(os.path.join(path,'OUTCAR'),'r')
    spin=True
    for frame in traj.frames:
        line=out.readline()
        while "magnetization (x)" not in line:
            if "General timing and accounting informations for this job" in line:
                spin=False
                break
            else:
                line=out.readline()
                #print(line)
        if spin==True:
            for i in range(3): 
                out.readline()
            for atom in frame.atoms:
                split=out.readline().split()
                atom.update_mu([split[1],split[2],split[3],split[4]])
    out.close()


