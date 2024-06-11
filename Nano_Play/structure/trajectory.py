from .tools import *
def dist_var(frames,atoms):
    """Returns numpy array of the distances between two atoms over Frames"""
    dist=[]
    bn=frames[0].atoms[atoms[0]].ele+'-'+str(atoms[0])+'--'+frames[0].atoms[atoms[1]].ele+'-'+str(atoms[1])
    for frame in frames:
        dist.append(bond_dist(frame.atoms[atoms[0]],frame.atoms[atoms[1]]))
    return bn,np.asarray(dist)
def angle_var(frames,atom1,atom2,atom3):
    """Returns numpy array of the angles between two atoms over Frames"""
    angle_var=[]
    for frame in frames:
        angle_var.append(bond_ang(frame.atoms[atom1],frame.atoms[atom2],frame.atoms[atom3]))
    return np.asarray(angle_var) 
def rmsd(frames,ref_frame=0,atoms=[]):
    """This function will calculate the rmsd of given set of Frames"""
    ref=frames[ref_frame]
    rmsd=[]
    if len(atoms) == 0:
        atoms=[x.ind for x in frames[0].atoms]
    for frame in frames:
        sd=0.0
        for i in range(len(frame.atoms)):
            if frame.atoms[i].ind in atoms:
                dr=frame.atoms[i].xyz-ref.atoms[i].xyz
                sd=dr[0]**2+dr[1]**2+dr[2]**2+sd
        rmsd.append(math.sqrt(sd/len(atoms)))
    return rmsd    
#def pd(frames):
    
