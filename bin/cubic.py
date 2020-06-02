import sys
import random
from config import Config 
from ion import Ion 

if __name__=="__main__":
    atype=["A"]
    nn=8
    d=1./(nn)
    nion=4*(nn**3)
    print(nion)
    print(nn)
    itype=[nion]
    cell=nn*5.4
    u=[cell,0,0]
    v=[0,cell,0]
    w=[0,0,cell]
    ions=[]	
    conf = Config(ions=ions,u=u,v=v,w=w,system='gen_cubic',\
                  nion=nion,ntype=len(atype),\
                  types=atype, natmpertype = itype, coord_format='Direct')
    
    for ia in range(conf.nion):
        conf.ions.append (  Ion ( index_ion=ia ) )
    conf.typeinfo_init()
    # random structure in the box
    ia=0
    for i in range(nn):
        for j in range(nn):
            for k in range(nn):
                x = i*d
                y = j*d
                z = k*d
                conf.move(ia,pos=[x,y,z])
                ia+=1
                x = (i+0.5)*d
                y = (j+0.5)*d
                z = k*d
                conf.move(ia,pos=[x,y,z])
                ia+=1
                x = (i+0.5)*d
                y = j*d
                z = (k+0.5)*d
                conf.move(ia,pos=[x,y,z])
                ia+=1
                x = i*d
                y = (j+0.5)*d
                z = (k+0.5)*d
                conf.move(ia,pos=[x,y,z])
                ia+=1
    conf.write_POSFF('POSFF')
    print("POSFF generated")
