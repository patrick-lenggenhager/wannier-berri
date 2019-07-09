#------------------------------------------------------------#
# This file is distributed as part of the Wannier19 code     #
# under the terms of the GNU General Public License. See the #
# file `LICENSE' in the root directory of the Wannier19      #
# distribution, or http://www.gnu.org/copyleft/gpl.txt       #
#                                                            #
# The Wannier19 code is hosted on GitHub:                    #
# https://github.com/stepan-tsirkin/wannier19                #
#                                                            #
#------------------------------------------------------------#
#                                                            #
# written by                                                 #
#           Stepan Tsirkin, University of Zurich             #
#                                                            #
#------------------------------------------------------------#

#
# vectorized version
#

## a version that takes into account the previous calculated Ef

import numpy as np
from lazy_property import LazyProperty
import tetra_aux as taux
from tetra_tmp import Tetrahedra_tmp


class Tetrahedra():

    def __init__(self,data,latt_type="",reclattice=None):
    ### Some variables, defining my thtrahera scheme.

        if  latt_type=="bcc1":
            self.__TETRA_NEIGHBOURS,self._TETRA_WEIGHTS=taux.construct_tetra_bcc1(reclattice)
        elif latt_type=="bcc2":
            self.__TETRA_NEIGHBOURS,self._TETRA_WEIGHTS=taux.construct_tetra_bcc2(reclattice)
        else:
            self.__TETRA_NEIGHBOURS,self._TETRA_WEIGHTS=taux.construct_tetra_generic()
        

        self.__NEIGHBOURS=set(tuple(p) for t in self.__TETRA_NEIGHBOURS for p in t)
        
        
    
        E_neigh=dict()
        for shift in self.__NEIGHBOURS:
            E_neigh[shift]=data._get_E_K_shifted(shift)
        self._Etetra=np.zeros( data.E_K_only.shape+(self.Ntetra,4),dtype=float)
        for it,tetra in enumerate(self.__TETRA_NEIGHBOURS):
#            print ("tetra({0})=\n{1}\n".format(it,tetra))
            self._Etetra[:,:,it,0]=data.E_K_only
            for j in range(3):
                self._Etetra[:,:,it,j+1]=E_neigh[tuple(tetra[j])]

#        print ("Etetra: \n"+"\n".join("ik={0}\n".format(ik)+
#                "\n".join("ib={0}\n".format(ib)+
#                   str(EtetraB)
#                     for ib,EtetraB in enumerate(EtetraK)     )
#                     for ik,EtetraK in enumerate(self._Etetra[:3]) )+"\n\n\n") 

        self._ivertex=self._Etetra.argsort(axis=-1).argmin(axis=-1)
        self._E0=self._Etetra[:,:,:,0]
        self._Etetra=np.sort(self._Etetra,axis=-1)
        self._Etetra_min=self._Etetra[:,:,:,0].min(axis=-1)
        self._Etetra_max=self._Etetra[:,:,:,3].max(axis=-1)

        self.Ef_old=-np.Inf
        

        self._sel_iv_1= (self._ivertex==0)
        self._sel_iv_2= (self._ivertex==1)
        self._sel_iv_3= (self._ivertex==2)
        self._sel_iv_4= (self._ivertex==3)
        self._sel_iv_123=np.logical_or(np.logical_or(self._sel_iv_1,self._sel_iv_2),self._sel_iv_3)
        self._sel_iv_234=np.logical_or(np.logical_or(self._sel_iv_2,self._sel_iv_3),self._sel_iv_4)


        self._Ef_old=-np.Inf
        self._occ_old=np.zeros(self._Etetra_min.shape,dtype=float)

           
    @LazyProperty
    def Ntetra(self):  
        return len(self.__TETRA_NEIGHBOURS)


    def get_occ(self,Ef_new):
        Tetrahedra_tmp(self,Ef_new).update_occ(self._occ_old)
        self._average_degen(self._occ_old)
        self._Ef_old=Ef_new
        return self._occ_old
        
    @LazyProperty
    def _degen_groups(self):
        degengroups=[]
        for ik,E in enumerate(self._E0):
            borders=np.hstack( ( [0], np.where( (E[1:]-E[:-1])>1e-5)[0]+1, [len(E)]) )
            degengroups.append([ (b1,b2) for b1,b2 in zip(borders,borders[1:]) if b2-b1>1])
        return degengroups


    def _average_degen(self,weights):
     # make sure that degenerate bands have same weights
        for ik,groups  in enumerate(self._degen_groups):
           for b1,b2 in groups:
               weights[ik,b1:b2]=weights[ik,b1:b2].mean()


    