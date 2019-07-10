#!/usr/bin/env python2
DO_profile=True

import sys
sys.path.append('../../modules/')
import numpy as np
import get_data
from  berry import calcAHC
import functools
from parallel import eval_integral_BZ
from time import time
import tetra_update as tetra


def main():
    seedname="Fe"


    NKFFT=np.array([int(sys.argv[1])]*3)
    NKdiv=np.array([int(sys.argv[2])]*3)
    Efermi=np.linspace(12,13,51)

    t_start=time()
#    Data=get_data.Data(tb_file='Fe_tb.dat',getAA=True)
    Data=get_data.Data(seedname,getAA=True)
    
    t_read=time()
    eval_func=functools.partial(  calcAHC, Efermi=Efermi,tetra=True)
    AHC_all=eval_integral_BZ(eval_func,Data,NKdiv,NKFFT=NKFFT,parallel=False,nproc=4,bcc=True)
    t_calc=time()

    open(seedname+"_w19_ahc_fermi_scan.dat","w").write(
       "    ".join("{0:^15s}".format(s) for s in ["EF",]+
         [a+b for a in ["","J0_","J1_","J2_"] for b in ["x","y","z"]])+"\n"+
      "\n".join(
       "    ".join("{0:15.6f}".format(x) for x in [ef]+[x for X in ahc[(3,0,1,2),:] for x in X]) 
                      for ef,ahc in zip (Efermi,AHC_all) )
       +"\n")  

    print ("Time for reading data    (s) : ",t_read-t_start )
    print ("Time for calculating AHC (s) : ",t_calc-t_read  )
    print ("Total time               (s) : ",t_calc-t_start )

if __name__ == '__main__':
    if DO_profile:
        import cProfile
        cProfile.run('main()')
    else:
        main()



