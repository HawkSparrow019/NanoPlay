import os,sys
from Nano_Play.nanite import *


def write(nanite,resouces,methods,blocks,coordinate=None):
    inp=open(os.path.join(nanite.path,nanite.name+'.inp'),'w')
    print(resouces,file=inp)
    print(methods,file=inp)
    for block in blocks:
        print(block,file=inp)
    if coordinate == None:
        pass

