import os
import numpy as np
from   ..nanite import *
from  ..structure.read import contcar as read_contcar
def magmom(incar,ions):
    mom={'Co':2.0,'Mn':5.0,'Fe':4.0}
    line='   MAGMOM = '
    for ion in ions:
        if ion.element in mom:
            if ion.afm==True:
                line=line+str(ion.natoms)+'*'+'-'+str(mom.get(ion.element))+' '
            else:
                line=line+str(ion.natoms)+'*'+str(mom.get(ion.element))+' '

        else:
            line=line+str(ion.natoms)+'*0.0 '
    print(line,file=incar)
def ncl_magmom(incar,ions,direction):
    mom={'Co':2.0,'Mn':5.0,'Fe':4.0}
    line='   MAGMOM = '
    for ion in ions:
        if ion.element in mom:
            if direction == 'x':
                mu=str(mom.get(ion.element)) +' 0.0 0.0 '
            elif direction == 'y':
                mu='0.0 '+str(mom.get(ion.element)) +' 0.0 '
            elif direction == 'z':
                mu=' 0.0 0.0 '+str(mom.get(ion.element)) 
            for i in range(ion.natoms):
                    line=line+mu+' '
        else:
            for i in range(3):
                line=line+str(ion.natoms)+'*0.0 '
    print(line,file=incar)
    #return line
def hub_u(incar,ions):
    ldau={'Mn':4.00,'Fe':4.00,'Zn':8.50}
    ldaj={'Mn':1.00,'Fe':1.00,'Zn':1.00}
    ldaul,ldauu,ldauj='LDAUL = ','LDAUU = ','LDAUJ = '
    for ion in ions:
        if ion.element in ldau:
            ldaul=ldaul+'   2'
            ldauu=ldauu+'  '+str(ldau.get(ion.element))
            ldauj=ldauj+'  '+str(ldaj.get(ion.element))
        else:
            ldaul=ldaul+'  -2'
            ldauu=ldauu+'  0.00'
            ldauj=ldauj+'  0.00'
    print("LDAU = .TRUE. \nLDAUTYPE = 2 \n"+ldaul+"\n"+ldauu+"\n"+ldauj+"\nLDAUPRINT = 2",file=incar)
    #return "LDAU = .TRUE. \nLDAUTYPE = 2 \n"+ldaul+"\n"+ldauu+"\n"+ldauj+"\nLDAUPRINT = 2"
def write_general_input(path,specs):
    ions=specs.get('ions')
    incar=open(path+"INCAR",'w')
    print("System=  %s"%specs.get('sys'),file=incar)
    print("\n  NPAR     =   6 \n  LPLANE   =   .TRUE. \n  LSCALU   =   .FALSE. \n  NSIM     =   4 ! blocked algorithm update, four bands at a time",file=incar)
    print("#------Charge density, Restart and Wavefuncion----------",file=incar)
    print("\n#   ISTART   =  1 \n  # INIWAV   =   0\n#   ICHARG   =  11\n   NGXF    =  200\n   NGYF    =  200\n   NGZF    =  210",file=incar)
    print("\n#----DOS and band-----   ",file=incar)
    print("\n   NEDOS    =   3000\n   ISMEAR   =   1\n   SIGMA    =   0.1",file=incar)
    print("\n#----Tollerences ------",file=incar)
    print("\n   EDIFF    =   1E-06\n   EDIFFG   =  -1E-02",file=incar)
    print("\n#------Calculations Specifications----",file=incar)
    print("\n   IBRION   =   2\n   NSW      =   0\n   NELM     =   500\n   NELMIN   =   6\n   POTIM    =   0.5\n   ISIF     =   2 \n   PREC     =   high",file=incar)
    print("\n#--------  Charge and Magnetic Specifications  --------",file=incar)
    print("\n   ISPIN    =  2 \n #NUPDOWN  =   1.0 \n #NELECT = 41",file=incar)
    if specs.get('ncl') == True:
        print("\n#------- Non Colllinear Magnetism------(Adjust restart settings)--\n",file=incar)
        specs.update({'direction':input("Enter ncl direction : ")})
        print(ncl_magmom(ions,specs.get('direction')),file=incar)
        print("   NBANDS   =   \n   ISYM     =    0\n   LNONCOLLINEAR=.TRUE.\n   LSORBIT = .TRUE.\n   LMAXMIX = 4\n   GGA_COMPAT = .FALSE. \n   LORBMOM=.TRUE.\n   LMAXPAW=0\n ",file=incar)
    else:
        magmom(incar,ions)
    print("\n#---------Algorithms ---------------------------\n",file=incar)
    print("\n   IALGO    =  48 \n   LDIAG    =   .TRUE.\n   NELMDL   =   0\n   NOSTOPCAR = .TRUE.\n\n   LORBIT   =   12\n#   RWIGS    =   0.067 0.153 0.161 0.053 0.056 \n   BMIX_MAG =   0.600",file=incar)
    print("\n#---------------   for large cells ----------",file=incar)
    print("\n   LREAL = A     ! evaluate projection operators in real space",file=incar)
    print("\n#---------------  big memory Outputs --------------------------------",file=incar)
    print("\n   LCHARG   =  .TRUE.\n   LELF     =  .FALSE.\n   LWAVE    =  .TRUE.",file=incar)
    print("\n#----------Functionals and Dispersion corrections------------",file=incar)
    print("\n   GGA = PE \n   IVDW = 11",file=incar)
    if specs.get('u') == True:
        print("\n#--------LDA+U module------------------\n",file=incar)
        incar.write(hub_u(incar,ions))
        
