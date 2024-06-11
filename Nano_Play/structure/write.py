import os
import numpy as np
def xyz(atoms,xyzfile='mol.xyz',path='./'):
    """mol.xyz will be writeen for an Atom object array atoms"""
    fo=open(os.path.join(path,xyzfile),'w')
    print(len(atoms),"\n",file=fo)
    for i in range(len(atoms)):
        print("%2s   %15.9f  %15.9f  %15.9f"%(atoms[i].ele,atoms[i].xyz[0],atoms[i].xyz[1],atoms[i].xyz[2]),file=fo)
    fo.close()
def poscar(atoms,cell,poscar='POSCAR',path='./'):
    pos=open(os.path.join(path,poscar),'w')
    print('Name\n1.0000000',file=pos)
    for vec in cell:
        print("        %20.15f   %20.15f   %20.15f"%tuple(vec),file=pos)
    ele=[]
    for atom in atoms:
        if atom.ele not in ele:
            ele.append(atom.ele)
    count,ions=[],[]
    for e in ele:
        c,ion=0,[]
        for atom in atoms:
            if atom.ele == e:
                c=c+1
                ion.append(atom)
        ions.append(ion)    
        count.append(c)
    new_atoms=[]
    for ion in ions:
        for atom in ion:
            new_atoms.append(atom)
    line1,line2='  ','  '
    for i in range(len(ele)):
        line1=line1+'  '+str(ele[i])
        line2=line2+'  '+str(count[i])
    print(line1,'\n',line2,file=pos)
    if new_atoms[0].mob != None:
        print("Selective Dynamics",file=pos)
    print('Cartesian',file=pos)
    for atom in new_atoms:
        if atom.mob == None:
            print("%20.15f    %20.15f   %20.15f"%tuple(atom.xyz),file=pos)
        else:
            print("%20.15f    %20.15f   %20.15f  %s %s %s"%(tuple(atom.xyz)+tuple(atom.mob)),file=pos)

    pos.close()
def xyz_traj(traj,out_file='traj.xyz',path='./'):
    xyz=open(os.path.join(path,out_file),'w')
    for frame in traj.frames:
        print(len(frame.atoms),"\n",file=xyz)
        for atom in frame.atoms:
            print("%s   %12.8f   %12.8f  %12.8f"%(atom.ele,atom.xyz[0],atom.xyz[1],atom.xyz[2]),file=xyz)
def xdatcar(traj,out_file='XDATCAR',path='./'):
    xdat=open(os.path.join(path,out_file),'w')
    print("traj\n      1",file=xdat)
    for vec in traj.b_mat:
        print("    %16.11f   %16.11f  %16.11f"%tuple(vec),file=xdat)
        ele=[]
    line1,line2='',''
    for i in range(len(traj.elements)):
        line1=line1+'  '+str(traj.elements[i])
        line2=line2+'  '+str(traj.ele_nums[i])
    print(line1,'\n',line2,file=xdat)
    dc=1
    for frame in traj.frames:
        print("Direct configuration=     ",dc,file=xdat)
        for atom in frame.atoms:
            print("%14.8f    %14.8f   %14.8f"%tuple(atom.frac),file=xdat)
        dc=dc+1   
def vmd_vectors(v0:np.array,v1:np.array,color:str,drawing_object='arrow',point=False,path='./',filename='vectors.tcl',res=500,radius=0.08,cone=1/5,cone_radius=0.12):
    if filename in os.listdir(path):
        fo=open(os.path.join(path,filename),'a')
    else:
        fo=open(os.path.join(path,filename),'w')
        print("axes location Off;",file=fo)
        print("display height 3.1;\ncolor Display Background white;\ndisplay projection Orthographic;\ndisplay depthcue off;",file=fo)
        print("mol color Element;\nmol selection all;\nmol representation CPK 0.5 0.3 500.0 500.0;\nmol material Glossy;\nmol addrep 0;\nmol delrep 0 0;",file=fo)
    if point==True:
        v=v1-v0
    else:
        v=v1
    l=np.linalg.norm(v)
    uv=v/l
    if drawing_object == 'arrow': 
        v2=v0+((l*(1-cone))*uv)
        v3=v2+((l*cone)*uv)
        #v0=v0-(l*uv)
        print("draw color {};".format(color),file=fo)
        print("draw cylinder  {{ {} {} {} }} {{  {} {} {} }} resolution {} radius {};".format(v0[0],v0[1],v0[2],v2[0],v2[1],v2[2],res,radius),file=fo)
        print("draw cone  {{ {} {} {} }} {{  {} {} {} }} resolution {} radius {}".format(v2[0],v2[1],v2[2],v3[0],v3[1],v3[2],res,cone_radius),file=fo)
    elif drawing_object == 'dashed_line':
        v2=v0+v
        #v0=v0-v
        print("draw color {};".format(color),file=fo)
        print("draw line  {{ {} {} {} }} {{  {} {} {} }} width 4 style dashed;".format(v0[0],v0[1],v0[2],v2[0],v2[1],v2[2]),file=fo)
    elif drawing_object == 'line':
        v2=v0+v
        v0=v0-v
        print("draw color {};".format(color),file=fo)
        print("draw line  {{ {} {} {} }} {{  {} {} {} }} width 4;".format(v0[0],v0[1],v0[2],v2[0],v2[1],v2[2]),file=fo)

    fo.close()