import scipy as sp
from ..nanite import *
import os


###### Reading tools ####
def update_ml_resolved_population(atoms,block,which):
    """Updates the ml-resolved population to the atom objects."""
    set=False
    count=0
    split=block[0].split()
    split=split[2:]
    block_counter=1
    while count < len(atoms):
        atom=atoms[count]
        orbs=['s','pz','px','py','dz2','dxz','dyz','dx2y2','dxy','f0','f1','f-1','f2','f-2','f3','f-3','g0',
              'g1','g-1','g2','g-2','g3','g-3','g4','g-4','h0','h1','h-1','h2','h-2','h3','h-3','h4','h-4',
              'h5','h-5']
        for orb in orbs:
            if len(block) < block_counter+1:
                count=count+1
                break
            if len(split) <= 6:
                if orb not in atom.orbs:
                    atom.orbs.update({orb:Orbital(orb)})
                    atom.orbs.get(orb).as_atomic_orb()
                orbital=atom.orbs.get(orb)
                if which == 'lwdn_chg':
                    if len(split) == 3:
                        orbital.lwdn_chg=float(split[-1])
                    else:
                        orbital.lwdn_chg=float(split[-4])
                elif which == 'lwdn_spin':
                    if len(split) == 3:
                        orbital.lwdn_spin=float(split[-1])
                    else:
                        orbital.lwdn_spin=float(split[-4])
                elif which == 'mlkn_chg':
                    if len(split) == 3:
                        orbital.mlkn_chg=float(split[-1])
                    else:
                        orbital.mlkn_chg=float(split[-4])
                elif which == 'mlkn_spin':
                    if len(split) == 3:
                        orbital.mlkn_spin=float(split[-1])
                    else:
                        orbital.mlkn_spin=float(split[-4])
                split=block[block_counter].split()
                block_counter=block_counter+1
            else:
                split=split[2:]
                count=count+1
                break        
def update_composition(mos,data,en_occ=True):
    mo_ind=[int(x) for x in data[0].split()]
    if en_occ == False:
        en=[float(x) for x in data[1].split()]    
        occ=[float(x) for x in data[2].split()]
        for i in range(len(mo_ind)):
            mos[mo_ind[i]].en=en[i]
            mos[mo_ind[i]].occ=occ[i]
    for line in data[4:]:
        split=line.split()
        basis=split[0]+'-'+split[1]+'-'+split[2]
        compo=[float(x) for x in split[3:]]
        for i in range(len(mo_ind)):
            mos[mo_ind[i]].c.update({basis:compo[i]})

    for x in mo_ind:
        mos[x].find_major_contributors()

def update_TDDFT_state(line,spin,scale):
    """For sTDDFT aprroach."""
    split=line.split()        
    state=State(split[0],split[2])
    state.transitions.append(s_Transition(line.split("(")[1].split()[:-1],line.split("(")[0].split()[-1],spin,scale))
    state.transitions.append(s_Transition(line.split("(")[2].split()[:-1],line.split("(")[1].split()[-1],spin,scale))
    state.transitions.append(s_Transition(line.split("(")[3].split(),line.split("(")[2].split()[-1],spin,scale))
    return state
def update_sftddft_transition(state,line,mos):
    """In SFTDDFT method the initial orbital is always alpha and the final is always beta""" 
    split=line.split()
    init_ind=int(split[0][:-1])
    final_ind=int(split[2][:-1])
    compo=float(split[4])
    state.transitions.append(Transition(mos[0][init_ind],mos[1][final_ind],final_type='b',compo=compo))
def update_tddft_transition(state,line,mos):
    """Updates trasitions corresponding to a TDDFT State""" 
    split=line.split()
    if len(mos) == 2:
        mos={'a':mos[0],'b':mos[1]}
    else:
        mos={'a':mos}
    init_ind=int(split[0][:-1])
    final_ind=int(split[2][:-1])
    init_type=split[0][-1]
    final_type=split[2][-1]
    compo=float(split[4])
    state.transitions.append(Transition(mos.get(init_type)[init_ind],mos.get(final_type)[final_ind],init_type=init_type,final_type=final_type,compo=compo))
    
##### Molecular Orbital analysis Tools ###############
def find_fmo(mos):
    for mo in mos:
        if mo.occ < 0.50:
            lumo=mo.ind
            break
    homo=lumo-1
    return homo
def del_en(mos):
    homo=find_fmo(mos)
    del_en=mos[homo+1].en-mos[homo].en
    return del_en
def find_degenerate_mos(mos,tol=0.05):
    dg={mos[0].en:[mos[0]]}
    for mo in mos[1:]:
        prev=False    
        for x in dg:
            if mo.en-tol < x  and mo.en+tol > x :
                dg.get(x).append(mo)
                prev=True
                append=False
            else:
                append=True
        if append == True and prev != True:
            dg.update({mo.en:[mo]})                
    return dg 

def find_degenerate_states(mos,tol=0.05):
    """This will return a directory of State like objects with a common attribute Object.en
       The tolerance values are in """
    dg={mos[0].en:[mos[0]]}
    for mo in mos[1:]:
        prev=False    
        for x in dg:
            if mo.en-tol < x  and mo.en+tol > x :
                dg.get(x).append(mo)
                prev=True
                append=False
            else:
                append=True
        if append == True and prev != True:
            dg.update({mo.en:[mo]})                
    return dg 

### spectra analysis tools
def normalise(spectra,upper=None,lower=None):
    """Normailse the Spectra between the limits of lower and upper and returns the normalised Spectrum"""
    en_tmp,intsty_tmp=[],[]
    if upper == None and lower != None:
        for i in range(len(spectra.en)):
            if spectra.en[i] >= lower:
                en_tmp.append(spectra.en[i])
                intsty_tmp.append(spectra.intsty[i])
    elif upper != None and lower == None:
        for i in range(len(spectra.en)):
            if spectra.en[i] <= upper:
                en_tmp.append(spectra.en[i])
                intsty_tmp.append(spectra.intsty[i])
    elif upper != None and lower != None:
        for i in range(len(spectra.en)):
            if spectra.en[i] <= upper and spectra.en[i] >= lower :
                en_tmp.append(spectra.en[i])
                intsty_tmp.append(spectra.intsty[i])
    else:
        en_tmp=spectra.en
        intsty_tmp=spectra.intsty
    return Spectrum(en_tmp,np.asarray(intsty_tmp)/max(intsty_tmp),spectra.type)
def broden(spectra,sigma,upper=None,lower=None):
    """Broaden the Spectra between the limits of lower and upper and returns brodened Spectrum"""
    ranges={'UVVIS':[100.00,800.00],'IR':[400.00,4000.00],'ir':[400.00,4000.00],'XAS':[7100,7150]}
    en_tmp,intsty_tmp=[],[]
    if upper == None and lower != None:
        upper=ranges.get(spectra.type)[1]
    elif upper != None and lower == None:
        lower=ranges.get(spectra.type)[0]
    elif upper == None and lower == None:            
        lower=ranges.get(spectra.type)[0]
        upper=ranges.get(spectra.type)[1]
    for i in range(len(spectra.en)):
        if spectra.en[i] <= upper and  spectra.en[i] >= lower:
            en_tmp.append(spectra.en[i])
            intsty_tmp.append(spectra.intsty[i])
    intsty_tmp.append(0.00)
    print(en_tmp)
    x=np.linspace(upper,lower,2000, endpoint=True)
    gE=[]
    for Ei in x:
        tot=0
        for Ej,os in zip(en_tmp,intsty_tmp):
            tot+=os*np.exp(-((((Ej-Ei)/sigma)**2)))
        gE.append(tot)
    return Spectrum(x,gE,spectra.type)

### Population analysis ###
def orb_pop(atom,l,type='lwdn_chg'):
    """Returns the total population for a l value."""
    pop=0
    for ao in atom.orbs:
        if l in ao:
            if type == 'lwdn_chg':
                pop=pop+atom.orbs.get(ao).lwdn_chg
            #print(atom.ind,orb,atom.orbs.get(ao).lwdn_spin)
            if type == 'lwdn_spin':
                pop=pop+atom.orbs.get(ao).lwdn_spin
            if type == 'mlkn_chg':
                pop=pop+atom.orbs.get(ao).lwdn_spin
            if type == 'mlkn_spin':
                pop=pop+atom.orbs.get(ao).lwdn_spin
    return pop
def update_mbo_data(block,mol):
    """Updates the mayer bond order for the available bonds of the Nanite object"""
    #print(block)
    for line in block:
        #print(line)
        split=line[:-1].split('B(')[1:]
        #print(split)
        for bond in split:
            s1=bond.split(') :')
            s2=s1[0].split(',')
            mol.bonds.append(Bond(atom1=mol.atoms[int(s2[0].split('-')[0])],atom2=mol.atoms[int(s2[1].split('-')[0])]))#,name=s2[0].strip()+'-'+s2[1].strip()))
            mol.bonds[-1].bo=float(s1[1])


### Multi-reference ###
def update_roots(block,mlt):
    roots,i=[],0
    for line in block:
        if 'ROOT' in line:
            roots.append(Root(i,line.split()[3],mlt))
            i=i+1
        else:
            split=line.split()
            roots[-1].cfgs.append(Cfg(split[0],split[-1]))
    return roots

