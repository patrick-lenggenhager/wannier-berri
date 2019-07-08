#
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
#------------------------------------------------------------

import numpy as np
import wan_ham as wham
import lazy_property
from get_data import Data
import tetra_vec as tetra

class Data_dk(Data):
    def __init__(self,data,dk=None,AA=None,BB=None,CC=None,SS=None,NKFFT=None,components_1d=(0,1,2),NKdiv=None):
        self.spinors=data.spinors
        self.iRvec=data.iRvec
        self.real_lattice=data.real_lattice
        self.NKdiv=NKdiv    
        self.NKFFT=data.NKFFT if NKFFT is None else NKFFT
        self.num_wann=data.num_wann
        # components_1d should be  a tuple of indices  0,1,2
        self.components_1d=components_1d
        self.alpha=np.array([wham.alpha[c] for c in self.components_1d])
        self.beta =np.array([wham.beta [c] for c in self.components_1d])

        if dk is not None:
            expdk=np.exp(2j*np.pi*self.iRvec.dot(dk))
        else:
            expdk=np.ones(self.nRvec)


        self.HH_R=data.HH_R[:,:,:]*expdk[None,None,:]
        
        if AA in (None,True):
            try:
                self.AA_R=data.AA_R[:,:,:,:]*expdk[None,None,:,None]
            except AttributeError:
                if AA : raise AttributeError("AA_R is not defined")

    @lazy_property.LazyProperty
    def ncomp1d(self):
       return len(self.components_1d)


    @lazy_property.LazyProperty
    def _get_eig_deleig(self):
#        print "running _get_eig_deleig"
        self._E_K,self._delE_K, self._UU_K, self._HH_K, self._delHH_K =   wham.get_eig_deleig(self.NKFFT,self.HH_R,self.iRvec,self.cRvec)

    @lazy_property.LazyProperty
    def NKFFT_tot(self):
        return np.prod(self.NKFFT)
    
    def _get_E_K_shifted(self,shift):
    # shift is given in terms of "small" k-grid as integres
    # TODO : check if the shifted grid belongs to the original or one of calculated neighbour grids
    # So far, let's recalculate them anyway
        dk=np.array(shift)/(self.NKFFT*self.NKdiv)
        return Data_dk(self,dk=dk,AA=False,BB=False,CC=False,SS=False,NKFFT=self.NKFFT,NKdiv=self.NKdiv).E_K_only

    @lazy_property.LazyProperty
    def E_K_neighbours(self):
        ekneigh=dict()
        for shift in tetra.NEIGHBOURS:
            ekneigh[shift]=self._get_E_K_shifted(shift)
#        ekneigh_K=[
#           dict( {shift:ekneigh[shift][ik] for shift in tetra.NEIGHBOURS})
#              for ik in range(self.NKFFT_tot)
#                   ]
        return ekneigh

    def get_occ_tetra(self,Ef):
        print ("evaluating occ tetra for Ef=",Ef)
        return tetra.get_occ(self.E_K,self.E_K_neighbours,Ef)
        return occ
            

    @lazy_property.LazyProperty
    def E_K(self):
        self._get_eig_deleig
        return self._E_K

    @lazy_property.LazyProperty
    def E_K_only(self):
        return wham.get_eig(self.NKFFT,self.HH_R,self.iRvec)


    @lazy_property.LazyProperty
    def delE_K(self):
        self._get_eig_deleig
        return self._delE_K

    @lazy_property.LazyProperty
    def UU_K(self):
        self._get_eig_deleig
        return self._UU_K


    @lazy_property.LazyProperty
    def UUC_K(self):
        return self.UU_K.conj()


    @lazy_property.LazyProperty
    def HH_K(self):
        self._get_eig_deleig
        return self._HH_K

    @lazy_property.LazyProperty
    def delHH_K(self):
        self._get_eig_deleig
        return self._delHH_K

    @lazy_property.LazyProperty
    def delHH_dE_K(self):
            _delHH_K_=np.einsum("kml,kmna,knp->klpa",self.UU_K.conj(),self.delHH_K,self.UU_K)
#          The following is probalby faster:
#            _delHH_K_=np.array([    uu.dot(uuc.dot(delhh))
#                       for uu,delhh,uuc in 
#                             zip(self.UU_K.transpose(0,2,1),self.delHH_K,self.UUC_K.transpose(0,2,1))])
            dEig_threshold=1e-14
            dEig=self.E_K[:,:,None]-self.E_K[:,None,:]
            select=abs(dEig)<dEig_threshold
            dEig[select]=dEig_threshold
            _delHH_K_[select]=0
            return 1j*_delHH_K_/dEig[:,:,:,None]

    @lazy_property.LazyProperty
    def delHH_dE_SQ_K(self):
         return  (self.delHH_dE_K[:,:,:,self.beta]
                         *self.delHH_dE_K[:,:,:,self.alpha].transpose((0,2,1,3))).imag

    @lazy_property.LazyProperty
    def delHH_dE_AA_K(self):
         return ( (self.delHH_dE_K[:,:,:,self.alpha]*self.AAUU_K.transpose((0,2,1,3))[:,:,:,self.beta]).imag+
               (self.delHH_dE_K.transpose((0,2,1,3))[:,:,:,self.beta]*self.AAUU_K[:,:,:,self.alpha]).imag  )

    @lazy_property.LazyProperty
    def UUU_K(self):
        return self.UU_K[:,:,None,:]*self.UUC_K[:,None,:,:]
    
    @lazy_property.LazyProperty
    def AAUU_K(self):
#        print "running get_AAUU_K.."
        _AA_K=wham.fourier_R_to_k_hermitian( self.AA_R,self.iRvec,self.NKFFT)
        return np.einsum("kml,kmna,knp->klpa",self.UUC_K,_AA_K,self.UU_K)

    @lazy_property.LazyProperty
    def OOmegaUU_K(self):
#        print "running get_OOmegaUU_K.."
        _OOmega_K =  wham.fourier_R_to_k_hermitian( -1j*(
                        self.AA_R[:,:,:,self.alpha]*self.cRvec[None,None,:,self.beta ] - 
                        self.AA_R[:,:,:,self.beta ]*self.cRvec[None,None,:,self.alpha])   , self.iRvec, self.NKFFT )
        return np.einsum("knmi,kmna->kia",self.UUU_K,_OOmega_K).real


unused="""


    def get_OOmega_K(self):
        try:
            return self._OOmega_K
        except AttributeError:
            print "running get_OOmega_K.."
            self._OOmega_K=    -1j* wham.fourier_R_to_k( 
                        self.AA_R[:,:,:,self.alpha]*self.cRvec[None,None,:,self.beta ] - 
                        self.AA_R[:,:,:,self.beta ]*self.cRvec[None,None,:,self.alpha]   , self.iRvec, self.NKFFT )
             
            return self._OOmega_K


    def get_AA_K(self):
        try:
            return self._AA_K
        except AttributeError:
            print "running get_AA_K.."
            self._AA_K=wham.fourier_R_to_k( self.AA_R,self.iRvec,self.NKFFT)
            return self._AA_K


"""