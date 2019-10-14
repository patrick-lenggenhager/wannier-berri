#!/usr/bin/env python3
DO_profile=True

import sys
sys.path.append('/Users/stepan/ARBEIT/wannier19/modules/')
import numpy as np
import get_data
import gyrotropic as gyro
import berry
import functools
from parallel import eval_integral_BZ

import cProfile

def main():
    seedname="Fe"
#NKFFT=    NKFFT=np.array([int(sys.argv[1])]*3)
    NKFFT=np.array([16,16,16])
    try:
        NKdiv=np.array([int(sys.argv[1])]*3)
    except:
        NKdiv=np.array([1]*3)
        
    Efermi=np.linspace(6.7,7.0,1001)
    Data=get_data.Data(seedname,getAA=True)

    eval_func=functools.partial(  gyro.calcAHC, Efermi=Efermi ,degen_thresh=0.01)
    AHC_all=eval_integral_BZ(eval_func,Data,NKdiv,NKFFT=NKFFT,parallel=False,nproc=8)
    print ("AHC= {}".format(AHC_all))

    eval_func=functools.partial(  berry.calcAHC, Efermi=Efermi )
    AHC_all=eval_integral_BZ(eval_func,Data,NKdiv,NKFFT=NKFFT,parallel=True,nproc=8)
    print ("AHC= {}".format(AHC_all))

#    exit()
#    AHC_all=np.zeros((len(Efermi),3))    
## now write the result
#    open(seedname+"_w19_ahc_fermi_scan_NK={0}-{1}-{2}.dat".format(*(tuple(NKFFT*NKdiv))),"w").write(
#       "    ".join("{0:^15s}".format(s) for s in ["EF","x","y","z"])+"\n"+
#      "\n".join(
#       "    ".join("{0:15.6f}".format(x) for x in [ef,ahc[0],ahc[1],ahc[2]]) 
#                      for ef,ahc in zip (Efermi,AHC_all) )
#       +"\n")  
          


cProfile.run('main()')



