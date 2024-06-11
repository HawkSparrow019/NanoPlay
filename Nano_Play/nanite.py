from matplotlib.colors import Normalize
import numpy as np
from scipy.integrate import simps as intgr

class Nanite:
    def __init__(self,path='./',name='Nanite',spin=True):
        """Nanite object will be created using the path of the output directory and its spin polarization """
        self.name=name
        self.path=path
        self.spin=spin
        self.atoms=[]
        self.tot_en=0.0
        self.bonds=[]
    def update_trajectory(self,trajectory):
        self.traj=trajectory
    def as_molecule(self,mr=False):
        if self.spin == True:
            self.s2=0.0
            self.mlkn_chg=0.0
            self.mlkn_spin=0.0
            self.lwdn_chg=0.0
            self.lwdn_spin=0.0
            self.hf_chg=0.0
            self.hf_spin=0.0
            self.a_mos=[]
            self.b_mos=[]
        else:
            self.mlkn_chg=0.0
            self.lwdn_chg=0.0
            self.hf_chg=0.0
            self.mos=[]
        if mr == True:
            self.mr=True
            self.multis={} #{multiplicity:Roots}
            self.roots=[]
            self.act_mos=[]
    def as_nano(self,ncl=False):
        self.cell=np.asarray([np.zeros(3),np.zeros(3),np.zeros(3)])
        self.layers=[]
        self.ions=[]
        self.charge=0.0
        self.nelect=0.0
        self.inv=False
        self.ncl=ncl
        if self.ncl == False:
            self.mu=0.0
        else:
            self.mu=np.zeros(3)
        
    def update_es(self,en,e_f,dos):
        self.en=np.asarray(en)
        self.e_f=float(e_f)
        self.e_m_e_f=self.en-self.e_f
        if self.spin == True:
            self.up_dos = np.asarray(dos[0])
            self.dwn_dos = np.asarray(dos[1])
        else:
            self.dos=np.asarray(dos)
    def bader_population(self):
        self.bader_pop=0
        if self.atoms[0].bader_pop != None:
            for atom in self.atoms:
                self.bader_pop=self.bader_pop+atom.bader_pop
    def spectra(self,spectrum):
        if spectrum.type == 'ir':
            self.ir=spectrum
        elif spectrum.type == 'raman':
            self.raman=spectrum
        elif spectrum.type == 'XAS' or spectrum.type == 'xas':
            self.xas=spectrum        
    def update_bonds(self):
        rad={'X':2.00,'H':0.56,'C':0.9,'N':1.0,'O':1.0,'F':1.5,'Cl':1.75,'Br':1.85,'I':1.98,'Ti':1.6,'Mn':1.35,'Fe':1.35,'Co':1.30,'Au':1.4,'Cu':1.45,'V':1.55,'Cr':1.55,'Mo':1.45,'W':1.46}
        for i in range(len(self.atoms)):
            for j in range(i+1,len(self.atoms)):
                bond=Bond(self.atoms[i],self.atoms[j])
                if bond.bd < rad.get(self.atoms[i].ele)+rad.get(self.atoms[j].ele):
                    self.bonds.append(bond)
    def update_uvvis(self,states):
        """This function must be removed in future for simplicity"""
        self.states=states
        lam,intsty=[],[]
        for state in self.states:
            lam.append(state.lam)
            intsty.append(state.intsty)
        self.uvvis=Spectrum(lam,intsty,'UVVIS')
        norm=np.asarray(self.uvvis.intsty)/max(self.uvvis.intsty)
        for i  in range(len(self.states)):
            self.states[i].norm_intsty=norm[i]
class Ion:
    """"""
    def __init__(self,element):
        self.ele=element
        self.natoms=0
        self.atoms=[]
    def mu(self):
        self.mu=0
        for atom in self.atoms:
            self.mu=self.mu+atom.mu
    def bader_population(self):
        self.bader_pop,self.bader_chg=0,0
        if self.atoms[0].bader_pop != None:
            for atom in self.atoms:
                self.bader_pop=self.bader_pop+atom.bader_pop
                #self.bader_chg=self.bader_chg+atom.bader_chg
class Layer:
    def __init__(self,name):
        self.name=name
        self.atoms=[]
        self.pos=0.0
class Atom:
    """Atom object."""
    def __init__(self,index=0,symbol='X'):
        mass={'X':0.0,'H':1.0,'C':12.00,'N':14.00,'Fe':55.85,'Au':197.00}
        self.ind=int(index)
        self.ele=symbol
        self.xyz=np.zeros(3)
        self.frac=np.zeros(3)
        self.mass=mass.get(symbol)
        self.orbs={}
        self.mob=None
    def update_elec_prop(self,chg=None,magmom=None,bader_pop=None,
       lwdn_chg=None,lwdn_spin=None,mlkn_chg=None,mlkn_spin=None,hf_chg=None,hf_spin=None):
        """This will update the electronic properties of a atom in a molecule"""
        self.chg=chg
        self.magmom=magmom
        self.lwdn_chg=lwdn_chg
        self.lwdn_spin=lwdn_spin
        self.mlkn_chg=mlkn_chg
        self.mlkn_spin=mlkn_spin
        self.hf_chg=hf_chg
        self.hf_spin=hf_spin
        self.bader_pop=bader_pop
    def update_mu(self,mu):
        self.s_mu=float(mu[0])
        self.p_mu=float(mu[1])
        self.d_mu=float(mu[2])
        self.mu=float(mu[3])
    def update_ncl_magnetization(self,mu):
        self.s_mu=np.asarray(mu[0])
        self.p_mu=np.asarray(mu[1])
        self.d_mu=np.asarray(mu[2])
        self.mu=np.asarray(mu[3])

    def update_position(self,pos):
        self.xyz[0]=float(pos[0])
        self.xyz[1]=float(pos[1])
        self.xyz[2]=float(pos[2])
    def as_mobile(self,mob):
        self.mob=mob
class Bond:
    """A Bond object is beween two atoms it will have the properties like bond length, bond order etc."""
    def  __init__(self,atom1:Atom,atom2:Atom,name=None):
        rad={'X':2.00,'H':0.56,'C':0.9,'N':1.0,'O':1.0,'F':1.5,'Cl':1.75,'Br':1.85,'I':1.98,'Ti':1.6,'Mn':1.35,'Fe':1.35,'Co':1.30,'Au':1.4,'Cu':1.45,'V':1.55,'Cr':1.55,'Mo':1.45,'W':1.46}
        if name == None:
            self.name=str(atom1.ind)+'-'+atom1.ele+'-'+str(atom2.ind)+'-'+atom2.ele
            self.bd=np.linalg.norm(atom1.xyz-atom2.xyz)
        else:
            self.name = name
            self.bd=np.linalg.norm(atom1.xyz-atom2.xyz)
        self.bo=0.0
class Orbital:
    def __init__(self,name,spin=True,ncl=False):
        """This will create an orbital object with the name of the orbital"""
        syms={'s':'$\mathregular{ s }$','py':'$\mathregular{p_{y} }$','pz':'$\mathregular{p_{z} }$','px':'$\mathregular{p_{x} }$',
          'dxy':"$\mathregular{d_{xy} }$",'dyz':"$\mathregular{d_{yz} }$",'dz2':"$\mathregular{d_{z^2} }$",
          'dxz':"$\mathregular{d_{xz} }$",'dx2y2':"$\mathregular{d_{x^2-y^2} }$",'dx2-y2':"$\mathregular{d_{x^2-y^2} }$"}
        self.name=str(name)
        if name not in syms:                   
            self.sym=name
        else:
            self.sym=syms.get(name)
        self.spin = spin
        self.ncl = ncl
    def as_atomic_orb(self):    
        self.lwdn_chg=0.0
        self.lwdn_spin=0.0
        self.mlkn_chg=0.0
        self.mlkn_spin=0.0
        self.hf_chg=0.0
        self.hf_spin=0.0
    def update_dos(self,dos):
        if len(dos) == 2:
            self.up_dos=np.asarray(dos[0])
            self.down_dos=-np.asarray(dos[1])
            self.spin=True
        elif len(dos) == 4:
            self.dos=np.asarray(dos[0])
            self.xdos=np.asarray(dos[1])
            self.ydos=np.asarray(dos[2])
            self.zdos=np.asarray(dos[3])
            self.ncl =True
        else:
            self.spin=False
            self.dos=np.asarray(dos)
    def __add__(self,orb):
        """Adds two orbital and form a hybrid orbital"""
        new_orb=Orbital(self.name+'+'+orb.name,spin=orb.spin,ncl=orb.ncl)
        new_orb.sym=self.sym+'+'+orb.sym
        if orb.spin == True and orb.ncl == False:
            new_orb.update_dos([(self.up_dos+orb.up_dos),-(self.down_dos+orb.down_dos)])
        elif orb.spin == True and orb.ncl == True:
            new_orb.update_dos([self.dos+orb.dos,self.xdos+orb.xdos,self.ydos+orb.ydos,self.zdos+orb.zdos])
        else:
            new_orb.update_dos((self.dos+orb.dos))
        return new_orb
    def occupancy(self,e_m_e_f):
        """Calculate the occupancy for an orbital"""
        delta=e_m_e_f[1]-e_m_e_f[0]
        yup,ydwn,y,yx,yy,yz=[],[],[],[],[],[]
        for q in range(len(e_m_e_f)):
            if e_m_e_f[q] <= 0.0:
                if self.spin ==True and self.ncl == False:
                    yup.append(self.up_dos[q])
                    ydwn.append(self.down_dos[q])
                elif self.spin ==True and self.ncl == True:
                    yx.append(self.xdos[q])
                    yy.append(self.ydos[q])
                    yz.append(self.zdos[q])                
                else:
                    y.append(self.dos[q])
        if self.spin==True and self.ncl == False:
            self.up_occ=intgr(yup,dx=delta)
            self.down_occ=-intgr(ydwn,dx=delta)
            self.occ=self.up_occ+self.down_occ
        elif self.spin == True and self.ncl == True:
            self.m_x=intgr(yx,dx=delta)
            self.m_y=intgr(yy,dx=delta)
            self.m_z=intgr(yz,dx=delta)
            self.occ=self.m_x+self.m_y+self.m_z
        else:
            self.occ=intgr(y,dx=delta)
    def normalise(self,nf):
        if self.spin == True and self.ncl == False:
            self.nup_dos=self.up_dos/nf
            self.ndown_dos=self.down_dos/nf
        elif self.spin == True and self.ncl == True:
            self.ndos=self.dos/nf
            self.nxdos=self.xdos/nf
            self.nydos=self.ydos/nf
            self.nzdos=self.zdos/nf
        else:
            self.ndos=self.dos/nf

class MO:
    """ The MO object has to initialised with index of the MO,it's energy and occupancy. 
    The type(alpha or beta) is not going to help here"""
    def __init__(self,index,occupancy,energy):
        """Defination of the MO"""
        self.ind=int(index)
        self.occ=float(occupancy)
        self.en=float(energy)
        self.c={}
    def find_major_contributors(self,cutoff=2.0):
        """Finds the major atomic orbital contributor to the MO"""
        self.rc={}
        orbs=[x for x in self.c]
        compo_mat=[self.c.get(x) for x in self.c]
        self.mao=orbs[compo_mat.index(max(compo_mat))]
        self.maoc=max(compo_mat)
        for x in self.c:
            if self.c.get(x)>=cutoff:
                self.rc.update({x:self.c.get(x)})
class Frame:
    """This will create a frame object. """
    def __init__(self,n_frame):
        self.nf=int(n_frame)
        self.magmom=0.00
        self.temp=0.00
        self.en=0.00
        self.free_en=0.00
        self.atoms=[]
class Trajectory:
    """Default timestep 1.0 fs. PBC is stored in b_mat"""
    def __init__(self,timestep=1.0,b_mat=np.asarray(np.zeros(3))):
        #print("The default timestep for the trajectory is set to 1 fs.")
        self.timestep=float(timestep)
        self.b_mat=b_mat
        self.frames=[]
        self.elements=[]
        self.ele_nums=[]
    def update(self):  
        self.time=np.asarray([x.nf for x in self.frames])*self.timestep
        self.temp=np.asarray([x.temp for x in self.frames])
        self.en=np.asarray([x.en for x in self.frames])
        self.free_en=np.asarray([x.free_en for x in self.frames])
        self.magmom=np.asarray([x.magmom for x in self.frames])
class State:
    """The state object will contain information about a peak in the spectra. 
     This contains wavelength of the spectral peak, intensity and 
     the corresponding transitions and their contribution."""
    def __init__(self,index,lam):
        self.ind=int(index)
        self.lam=float(lam)
        self.transitions=[]
        self.intsty=0.0
        self.norm_intsty=0.0
class s_Transition:
    def __init__(self,trans,contri,spin,scale):
        if spin==True:
            spin_mo={'a':'a','b':'b'}
            self.init_ind=int(trans[0][:-3])+scale
            self.init_type=spin_mo.get(trans[0][-3])
            self.final_ind=int(trans[1][:-2])+scale
            self.final_type=spin_mo.get(trans[1][-2])
            self.contri=float(contri)*100
        else:
            self.init_ind=int(trans[0][:-2])+scale
            self.final_ind=int(trans[1][:-1])+scale
            self.contri=float(contri)*100
class Transition:
    """"""
    def __init__(self,init_mo,final_mo,spin=True,init_type='a',final_type='a',compo=0.0):
        if spin == True:
            self.init=init_mo
            self.final=final_mo
            self.init_type=init_type
            self.final_type=final_type
            self.compo=compo
        else:
            self.init=init_mo
            self.final=final_mo
            self.compo=compo
class Spectrum:
    """The base class for a spectrum. The normalization and broadeing has to be done"""
    def __init__(self,en,intsty,type):
        self.type=type
        self.en=np.asarray(en)
        self.intsty=np.asarray(intsty)
class Root:
    """Root object of a multireference calculation"""
    def __init__(self,index,energy,mlt):
        self.ind=int(index)
        self.en=float(energy)
        self.mlt=int(mlt)
        self.cfgs=[]
        self.D=0.0
        self.E=0.0
class Cfg:
    """Electronic configuration that constitutes a root for a MR object"""
    def __init__(self,weight,configuration):
        self.w=float(weight)
        self.c=[int(x) for x in configuration]
        self.mlt=self.c.count(1)+1