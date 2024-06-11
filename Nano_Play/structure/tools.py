from ..nanite import *
import numpy as np
import math

def angle(v1,v2,unit='deg'):
    """Return the angle between two vectors v1 and v2"""
    cos_theta = np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
    #print(cos_theta)
    sin_theta = np.linalg.norm(np.cross(v1, v2))
    theta=np.arccos(cos_theta)
    #print(theta*180/3.1416 ) 
    #theta = np.arctan2(sin_theta, cos_theta)
    if unit == 'deg':
        return 180.0 * theta / np.pi
    elif unit == 'rad':
        return theta 

def bond_ang(atom1,atom2,atom3,unit='deg'):
    """Return angle between three Atoms"""
    rij = atom1.xyz - atom2.xyz
    rkj = atom3.xyz - atom2.xyz
    return angle(rij,rkj,unit)
def bond_dist(atom1,atom2):
    """Returns distance between two Atoms"""
    dist=np.linalg.norm(atom1.xyz-atom2.xyz)
    return dist
def com(atoms):
    """This function returns the Center of Mass of a set of Atoms""" 
    sum_mi_ri,sum_mi=np.zeros(3),0.0
    for atom in atoms:
        sum_mi_ri=sum_mi_ri+(atom.mass*atom.xyz)
        sum_mi=sum_mi+atom.mass
    return sum_mi_ri/sum_mi
def plane(v1,v2,v3):
    """Returns the numpy array containing the coefficients of the equation 
    ax+by+cz+d=0 passing through the v1,v2,v3 points"""
    v21=v2-v1
    v31=v3-v1
    n=np.cross(v21,v31)
    d=-(np.matmul(n,np.transpose(v1)))
    return np.append(n,d)
def plane_dist(p,v):
    """Returns the distance between a point and a plane"""
    vdash=np.append(v,1)
    return (np.matmul(p,np.transpose(vdash)))/math.sqrt(((p[0]**2)+(p[1]**2)+(p[2]**2)))




def search_bond(nanite,ind=None,ele1=None,ele2=None,print_info=True):
    """This will search all the chemical bonds and print the possible bonds between """
    if ele1 == None:
        ele1=input("Enter the element you want to search bond with :")
    if ele2 == None:
        ele2= input("Enter the element you want to search bond with :")     
    bonds=[]
    for x in nanite.bonds:
        split=x.name.split('-')
        if (split[1]== ele1 and split[3] == ele2) or (split[1]== ele2 and split[3] == ele1) :
            if print_info == True:
                print("%s \t   %8.5f"%(x.name,x.bd))
            bonds.append(x)
    return bonds

def find_layers(nanite,element=None):
    """Search the atoms and update the atomic layers into the nanite object"""
    cutoff=1.5/np.linalg.norm(nanite.cell[2])
    layers=[Layer('L0')]
    layers[0].pos=nanite.atoms[0].frac[2]
    #for layer in layers:
    layered=[]
    for atom in nanite.atoms:
        if atom.ind not in layered:    
            found=False
            diff=[abs(x.pos-atom.frac[2]) for x in layers]
            for i in range(len(layers)):
                if diff[i] <=cutoff:
                    layers[i].atoms.append(atom)
                    found=True
                    break
            if found !=True:
                layers.append(Layer('L'+str(len(layers))))
                layers[-1].pos=atom.frac[2]
                layers[-1].atoms.append(atom)
                layered.append(atom.ind)        
    nanite.layers = layers       

