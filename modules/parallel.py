#------------------------------------------------------------#
# This file is distributed as part of the Wannier19 code     #
# under the terms of the GNU General Public License. See the #
# file `LICENSE' in the root directory of the Wannier19      #
# distribution, or http://www.gnu.org/copyleft/gpl.txt       #
#                                                            #
# The Wannier19 code is hosted on GitHub:                    #
# https://github.com/stepan-tsirkin/wannier19                #
#                     written by                             #
#           Stepan Tsirkin, University ofZurich              #
#                                                            #
#------------------------------------------------------------#



import  multiprocessing 
import functools
import numpy as np
from data_dk import Data_dk

dict_xyz={"x":0,"y":1,"z":2}




def eval_integral_BZ(func,Data,NKdiv=np.ones(3,dtype=int),parallel=False,
               nproc=1,NKFFT=None,components_1d="xyz",bcc=False):
    """This function evaluates in parallel or serial an integral over the Brillouin zone 
of a function func, which whould receive only one argument of type Data_dk, and return whatever object,
for which the '+' operation is defined.

the user has to provide 2 grids of K-points - FFT grid anf NKdiv

The parallelisation is done by NKdiv

As a result, the integration will be performed ove NKFFT x NKdiv
"""
    # trnslate xyz to indices
    components_1d=np.array([dict_xyz[c] for c in components_1d])

    NKFFT=Data.NKFFT if NKFFT is None else NKFFT
    dk1=1./(NKFFT*NKdiv)
    
    if bcc:
        from tetra_aux import get_bcc_shift
        dk2=get_bcc_shift(Data.reclattice)/(NKFFT*NKdiv)
        dk_list=[(dk1*np.array([x,y,z]),"bcc1") for x in range(NKdiv[0]) 
            for y in range(NKdiv[1]) for z in range(NKdiv[2]) ]
        dk_list=dk_list+[(dk[0]+dk2,"bcc2") for  dk in dk_list]
        paralfunc=functools.partial(
            _eval_func_dk2, func=func,Data=Data,NKFFT=NKFFT,components_1d=components_1d,NKdiv=NKdiv )
        print dk_list
        if parallel:
            p=multiprocessing.Pool(nproc)
            return sum(p.map(paralfunc,dk_list))/len(dk_list)
        else:
            return sum(paralfunc(dk) for dk in dk_list)/len(dk_list)
    else:
        dk_list=[dk1*np.array([x,y,z]) for x in range(NKdiv[0]) 
            for y in range(NKdiv[1]) for z in range(NKdiv[2]) ]
        paralfunc=functools.partial(
            _eval_func_dk, func=func,Data=Data,NKFFT=NKFFT,components_1d=components_1d,NKdiv=NKdiv )
        if parallel:
            p=multiprocessing.Pool(nproc)
            return sum(p.map(paralfunc,dk_list))/len(dk_list)
        else:
            return sum(paralfunc(dk) for dk in dk_list)/len(dk_list)


def _eval_func_dk(dk,func,Data,NKFFT,components_1d,NKdiv):
    data_dk=Data_dk(Data,dk[0],NKFFT=NKFFT,components_1d=components_1d,NKdiv=NKdiv)
    return func(data_dk)


def _eval_func_dk2(dk,func,Data,NKFFT,components_1d,NKdiv):
    data_dk=Data_dk(Data,dk[0],NKFFT=NKFFT,components_1d=components_1d,NKdiv=NKdiv,tetra_type=dk[1])
    return func(data_dk)

