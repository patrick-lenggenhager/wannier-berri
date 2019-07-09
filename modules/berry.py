#------------------------------------------------------------#
# This file is distributed as part of the Wannier19 code     #
# under the terms of the GNU General Public License. See the #
# file `LICENSE' in the root directory of the Wannier19      #
# distribution, or http://www.gnu.org/copyleft/gpl.txt       #
#                                                            #
# this file initially was  an adapted translation of         #
# the corresponding Fortran90 code from  Wannier 90 project  #
#                                                            #
# with significant modifications for better performance      #
#   it is nor a lot different                                #
#                                                            #
# The Wannier19 code is hosted on GitHub:                    #
# https://github.com/stepan-tsirkin/wannier19                #
#                                                            #
# The webpage of the Wannier90 code is www.wannier.org       #
# The Wannier90 code is hosted on GitHub:                    #
# https://github.com/wannier-developers/wannier90            #
#------------------------------------------------------------#
#                                                            #
#  Translated to python and adapted for wannier19 project by #
#           Stepan Tsirkin, University of Zurich             #
#                                                            #
#------------------------------------------------------------#

import numpy as np
from scipy import constants as constants
from collections import Iterable

fac_ahc = -1.0e8*constants.elementary_charge**2/constants.hbar


def eval_J0(A,occ):
    return A[occ].sum(axis=0)

def eval_J12(B,unoccocc):
    return -2*B[unoccocc].sum(axis=0)

def get_occ(E_K,Efermi):
    return (E_K< Efermi)
    
    
def calcAHC(data,Efermi=None,occ_old=None, evalJ0=True,evalJ1=True,evalJ2=True,smear=None,tetra=False):


    if not( smear is None) or tetra:
        return calcAHC_frac(data,Efermi=Efermi,occ_old=occ_old, evalJ0=evalJ0,evalJ1=evalJ1,evalJ2=evalJ2,smear=smear,tetra=tetra)


#version with integer (actually boolean) occupation numbers

    if occ_old is None: 
        occ_old=np.zeros((data.NKFFT_tot,data.num_wann),dtype=bool)
    
    ncomp=data.ncomp1d


    if isinstance(Efermi, Iterable):
        nFermi=len(Efermi)
        AHC=np.zeros( ( nFermi,4,ncomp) ,dtype=float )
        for iFermi in range(nFermi):
            AHC[iFermi]=calcAHC(data,Efermi=Efermi[iFermi],occ_old=occ_old, evalJ0=evalJ0,evalJ1=evalJ1,evalJ2=evalJ2)
        return np.cumsum(AHC,axis=0)
    
    # now code for a single Fermi level:
    AHC=np.zeros((4,ncomp))

    occ_new=get_occ(data.E_K,Efermi)
    unocc_new=np.logical_not(occ_new)
    unocc_old=np.logical_not(occ_old)

    changed_occ= occ_old!=occ_new
    selectK=np.where(np.any(changed_occ,axis=1))[0]
    selectB=np.where(np.any(changed_occ,axis=0))[0]
    if len(selectB)>0:
        Bmin=selectB.min()
        Bmax=selectB.max()+1
    else:
        return AHC

    occ_old_selk=occ_old[selectK]
    occ_new_selk=occ_new[selectK]
    unocc_old_selk=unocc_old[selectK]
    unocc_new_selk=unocc_new[selectK]
    delocc=occ_new_selk[:,Bmin:Bmax]!=occ_old_selk[:,Bmin:Bmax]
    unoccocc_plus=unocc_new_selk[:,:,None]*delocc[:,None,:]
    unoccocc_minus=delocc[:,:,None]*occ_old_selk[:,None,:]

    if evalJ0:
        AHC[0]= eval_J0(data.OOmegaUU_K[selectK,Bmin:Bmax], delocc)
    if evalJ1:
        B=data.delHH_dE_AA_K[selectK]
        AHC[1]=eval_J12(B[:,:,Bmin:Bmax],unoccocc_plus)-eval_J12(B[:,Bmin:Bmax,:],unoccocc_minus)
    if evalJ2:
        B=data.delHH_dE_SQ_K[selectK]
        AHC[2]=eval_J12(B[:,:,Bmin:Bmax],unoccocc_plus)-eval_J12(B[:,Bmin:Bmax,:],unoccocc_minus)
    AHC[3,:]=AHC[:3,:].sum(axis=0)

    occ_old[:,:]=occ_new[:,:]
    return AHC*fac_ahc/(data.NKFFT_tot*data.cell_volume)



def eval_J0_frac(A,occ):
    return np.sum(A*occ[:,:,None],axis=(0,1))

def eval_J12_frac(B,unoccocc):
    return -2*np.sum(B*unoccocc[:,:,:,None],axis=(0,1,2))

    
    
#version with frctional  occupation numbers (Fermi-Dirac or tetrahedra)
def calcAHC_frac(data,Efermi=None,occ_old=None, evalJ0=True,evalJ1=True,evalJ2=True,smear=0,diff_occ_threshold=1e-5,tetra=False):

    if occ_old is None: 
        occ_old=np.zeros((data.NKFFT_tot,data.num_wann),dtype=float)

    ncomp=data.ncomp1d


    if isinstance(Efermi, Iterable):
        nFermi=len(Efermi)
        AHC=np.zeros( ( nFermi,4,ncomp) ,dtype=float )
        for iFermi in range(nFermi):
            AHC[iFermi]=calcAHC_frac(data,Efermi=Efermi[iFermi],occ_old=occ_old, evalJ0=evalJ0,evalJ1=evalJ1,evalJ2=evalJ2,smear=smear,tetra=tetra)
        return np.cumsum(AHC,axis=0)
    
    # now code for a single Fermi level:
    AHC=np.zeros((4,ncomp))
    
    occ_new=data.get_occ(Efermi,tetra=tetra,smear=smear)
    
#    print ("sum of occupations for Ef={0} : {1}".format(Efermi,occ_new.sum()/data.NKFFT_tot ))

#    unocc_new=1.-occ_new
#    unocc_old=1.-occ_old
    changed_occ= np.abs(occ_old-occ_new)>diff_occ_threshold
    selectK=np.where(np.any(changed_occ,axis=1))[0]
    selectB=np.where(np.any(changed_occ,axis=0))[0]
    Bmin=selectB.min()
    Bmax=selectB.max()+1
    
    occ_old_selk=occ_old[selectK]
    occ_new_selk=occ_new[selectK]
    
    unocc_old_selk=1.-occ_old_selk
#    unocc_new_selk=1.-occ_new_selk

    delocc=occ_new_selk[:,Bmin:Bmax]-occ_old_selk[:,Bmin:Bmax]
        
    if evalJ1 or evalJ2:
        unoccocc_plus =unocc_old_selk[:,:,None]*delocc       [:,None,:]
        unoccocc_minus=delocc[:,:,None]        *occ_new_selk [:,None,:] 

    if evalJ0:
        AHC[0]= eval_J0_frac(data.OOmegaUU_K[selectK,Bmin:Bmax], delocc)
    if evalJ1:
        B=data.delHH_dE_AA_K[selectK]
        AHC[1]=eval_J12_frac(B[:,:,Bmin:Bmax],unoccocc_plus)-eval_J12_frac(B[:,Bmin:Bmax,:],unoccocc_minus)
    if evalJ2:
        B=data.delHH_dE_SQ_K[selectK]
        AHC[2]=eval_J12_frac(B[:,:,Bmin:Bmax],unoccocc_plus)-eval_J12_frac(B[:,Bmin:Bmax,:],unoccocc_minus)
    AHC[3,:]=AHC[:3,:].sum(axis=0)

    occ_old[:,:]=occ_new[:,:]
    return AHC*fac_ahc/(data.NKFFT_tot*data.cell_volume)
