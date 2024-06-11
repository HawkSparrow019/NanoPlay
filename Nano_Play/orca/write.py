from ..nanite import *
from ..structure import write as wt 
def summary(mol):
    print(mol.tot_en,mol.s2)
def reorient_gT(nanite):
    newatoms=[]
    for atom in nanite.atoms:
        newatoms.append(Atom(atom.ind,atom.ele))
        newatoms[-1].xyz=np.matmul(atom.xyz,nanite.gT)
    wt.xyz(newatoms,xyzfile='reorient_gT.xyz')
    nanite.rntd_atoms=newatoms
def gT_orientation(nanite):
    trans_gT=np.transpose(nanite.gT)
    gx,gy,gz=nanite.g[0]*trans_gT[0],nanite.g[1]*trans_gT[1],nanite.g[2]*trans_gT[2]    
    v0=nanite.atoms[0].xyz
    wt.vmd_vectors(v0,gx,'red',filename='g_vectors.tcl')
    wt.vmd_vectors(v0,gy,'blue',filename='g_vectors.tcl')
    wt.vmd_vectors(v0,gz,'green3',filename='g_vectors.tcl')
    for v in [gx,gy,gz]:
        nanite.atoms.append(Atom(index=len(nanite.atoms)))
        nanite.atoms[-1].xyz=v+nanite.atoms[0].xyz
def mag_axes(nanite,origin=0,cartesian=True):
    trans_gT=np.transpose(nanite.gT)
    trans_dT=np.transpose(nanite.dT)
    gx,gy,gz=nanite.g[0]*trans_gT[0],nanite.g[1]*trans_gT[1],nanite.g[2]*trans_gT[2]    
    dx,dy,dz=(2*nanite.g[0])*trans_dT[0],(2*nanite.g[1])*trans_dT[1],(2*nanite.g[2])*trans_dT[2] 
    v0=nanite.atoms[origin].xyz
    wt.vmd_vectors(v0,gx,'red',filename='g_vectors.tcl')
    wt.vmd_vectors(v0,gy,'blue',filename='g_vectors.tcl')
    wt.vmd_vectors(v0,gz,'green3',filename='g_vectors.tcl')
    wt.vmd_vectors(v0,dx,'red',drawing_object='dashed_line',filename='g_vectors.tcl')
    wt.vmd_vectors(v0,dy,'blue',drawing_object='dashed_line',filename='g_vectors.tcl')
    wt.vmd_vectors(v0,dz,'green3',drawing_object='dashed_line',filename='g_vectors.tcl')
    if cartesian== True:
        wt.vmd_vectors(v0,np.asarray([1,0,0]),'red',drawing_object='line',filename='g_vectors.tcl')
        wt.vmd_vectors(v0,np.asarray([0,1,0]),'blue',drawing_object='line',filename='g_vectors.tcl')
        wt.vmd_vectors(v0,np.asarray([0,0,1]),'green3',drawing_object='line',filename='g_vectors.tcl')

    #for v in [gx,gy,gz]:
    #    nanite.atoms.append(Atom(index=len(nanite.atoms)))
    #    nanite.atoms[-1].xyz=v+nanite.atoms[0].xyz
