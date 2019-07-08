#------------------------------------------------------------#
# This file is distributed as part of the Wannier19 code     #
# under the terms of the GNU General Public License. See the #
# file `LICENSE' in the root directory of the Wannier19      #
# distribution, or http://www.gnu.org/copyleft/gpl.txt       #
#                                                            #
# part of this file is based on                              #
# the corresponding Fortran90 code from                      #
#                                 Quantum Espresso  project  #
#                                                            #                                                            #
# The Wannier19 code is hosted on GitHub:                    #
# https://github.com/stepan-tsirkin/wannier19                #
#                                                            #
# The webpage of the QuantumEspresso  code is                #
#            https://www.quantum-espresso.org/               #
#------------------------------------------------------------#
#                                                            #
#  Translated to python and adapted for wannier19 project by #
#           Stepan Tsirkin, University ofZurich              #
#                                                            #
#------------------------------------------------------------#

#
# vectorized version
#

#! Copyright message from Quantum Espresso:
#! Copyright (C) 2016 Quantum ESPRESSO Foundation
#! This file is distributed under the terms of the
#! GNU General Public License. See the file `License'
#! in the root directory of the present distribution,
#! or http://www.gnu.org/copyleft/gpl.txt .
#!


# TODO : transform this module into a class, 
#            and save quantities that are  needed for different Efermi


import numpy as np
from lazy_property import LazyProperty


class Tetrahedra():

    def __init__(self,data):
    ### Some variables, defining my thtrahera scheme.

        self.__points=np.array([  [ 0, 0, 0],
                              [ 0, 0, 1],
                              [ 0, 1, 0],
                              [ 0, 1, 1],
                              [ 1, 0, 0],
                              [ 1, 0, 1],
                              [ 1, 1, 0],
                              [ 1, 1, 1] ])
        
        self.__tetra_1box=np.array([ [1,2,3,6],
                                     [2,3,4,6],
                                     [1,3,5,6],
                                     [3,4,6,8],
                                     [3,6,7,8],
                                     [3,5,6,7] ])-1
        
        self.__TETRA_NEIGHBOURS=[]
        

        for ix in 0,-1:
          for iy in 0,-1:
            for iz in 0,-1:
              shift=np.array([ix,iy,iz],dtype=int)
              for tetra in self.__tetra_1box:
                tet=[tuple(self.__points[ip]+shift) for ip in tetra]
                if tuple( (0,0,0) ) in tet:
                   self.__TETRA_NEIGHBOURS.append(np.array([t for t in tet if t!=(0,0,0)]))
        
            # now the class definition
        self.__NEIGHBOURS=set(tuple(p) for t in self.__TETRA_NEIGHBOURS for p in t)
        
        
    
        E_neigh=dict()
        for shift in self.__NEIGHBOURS:
            E_neigh[shift]=data._get_E_K_shifted(shift)
        self._Etetra=np.zeros( data.E_K_only.shape+(self.Ntetra,4),dtype=float)
        for it,tetra in enumerate(self.__TETRA_NEIGHBOURS):
            self._Etetra[:,:,it,0]=data.E_K_only
            for j in range(3):
                self._Etetra[:,:,it,j]=E_neigh[tuple(tetra[j])]

        self._ivertex=self._Etetra.argsort(axis=-1).argmin(axis=-1)
        self._E0=self._Etetra[:,:,:,0]
        self._Etetra=np.sort(self._Etetra,axis=-1)


        self._sel_iv_1= (self._ivertex==0)
        self._sel_iv_2= (self._ivertex==1)
        self._sel_iv_3= (self._ivertex==2)
        self._sel_iv_4= (self._ivertex==3)
        self._sel_iv_123=np.logical_or(np.logical_or(self._sel_iv_1,self._sel_iv_2),self._sel_iv_3)
        self._sel_iv_234=np.logical_or(np.logical_or(self._sel_iv_2,self._sel_iv_3),self._sel_iv_4)



           
    @LazyProperty
    def Ntetra(self):  
        return len(self.__TETRA_NEIGHBOURS)

    def _get_e01234(self,select):
        return (self._E0[select],)+tuple([E for E in self._Etetra[select].T])


    def _get_occ_range_43_(self,ef,e0,e1,e2,e3,e4,iv=None):
        c4 =  (e4 - ef) **3 / ((e4 - e1) * (e4 - e2) * (e4 - e3))
        dosef = 0.3 * (e4 - ef) **2 /((e4 - e1)* (e4 - e2) * (e4 - e3))*(e1 + e2 + e3 + e4 - 4. * e0 )
        if iv<=3: 
            return  1.  - c4 * (e4 - ef) / (e4 - e0) + dosef 
        else:
            return  1.  - c4 * (4. - (e4 - ef) * (1. / (e4 - e1) + 1. / (e4 - e2)  + 1. / (e4 - e3) ) ) + dosef 

    def _get_occ_range_32_(self,ef,e0,e1,e2,e3,e4,iv=None):
        c1 = (ef - e1) **2 / ((e4 - e1) * (e3 - e1))
        c2 = (ef - e1) * (ef - e2) * (e3 - ef)  / ((e4 - e1) * (e3 - e2) * (e3 - e1))
        c3 = (ef - e2) **2 * (e4 - ef) /( (e4 - e2)  * (e3 - e2) * (e4 - e1))
        dosef = 0.1 / (e3 - e1) / (e4 - e1) * (3. * 
                   (e2 - e1) + 6. * (ef - e2) - 3. * (e3 - e1 + e4 - e2) 
                   * (ef - e2) **2 / (e3 - e2) / (e4 - e2) )* (e1 + e2 +  e3 + e4 - 4. * e0 ) 
        if   iv==0:
            return  c1 + (c1 + c2) * (e3 - ef) / (e3 - e1) + (c1 + c2 + c3) * (e4 - ef) / (e4 - e1) + dosef
        elif iv==1:
            return  c1 + c2 + c3 + (c2 + c3)  * (e3 - ef) / (e3 - e2) + c3 * (e4 - ef) / (e4 - e2) + dosef 
        elif iv==2:
            return  (c1 + c2) * (ef - e1) / (e3 - e1) + (c2 + c3) * (ef - e2) / (e3 - e2) + dosef 
        elif iv==3:
            return  (c1 + c2 + c3) * (ef - e1)  / (e4 - e1) + c3 * (ef - e2) / (e4 - e2) + dosef 

    def _get_occ_range_21_(self,ef,e0,e1,e2,e3,e4,iv=None):
        c4 = (ef - e1) **3 / (e2 - e1) / (e3 - e1) / (e4 - e1)
        dosef = 0.3 * (ef - e1) **2 / (e2 - e1) / (e3 - e1) / (e4 - e1) * (e1 + e2 + e3 + e4 - 4. * e0 ) 
        if   iv==0:
            return   c4 * (4. - (ef - e1) * (1. / (e2 - e1) + 1. / (e3 - e1) + 1. / (e4 - e1) ) )   + dosef 
        elif iv in (1,2,3):
            return   c4 * (ef - e1) / (e0 - e1)  + dosef
    
    def _update_weights_(self,fun_get_occ,selectIV,iv,ef,weights):
        e0,e1,e2,e3,e4=self._get_e01234(selectIV)
        weights[selectIV]=fun_get_occ(ef,e0,e1,e2,e3,e4,iv)
  
    
    def get_occ(self,ef):
# TODO : save the quantities, which do not depend of Ef (or maybe not needed)
        weights=np.zeros(self._Etetra.shape[:-1])

##  First case :  Ef>=E4
        select= (ef>=self._Etetra[:,:,:,3])
        weights[select] = 1.
##  second case :  E4>Ef>=E3
        select= (ef<self._Etetra[:,:,:,3])*(ef>=self._Etetra[:,:,:,2])
        self._update_weights_( self._get_occ_range_43_ , select*self._sel_iv_123 , 1 , ef, weights)
        self._update_weights_( self._get_occ_range_43_ , select*self._sel_iv_4   , 4 , ef, weights)
##  third  case :  E3>Ef>=E2
        select= (ef<self._Etetra[:,:,:,2])*(ef>=self._Etetra[:,:,:,1])
        self._update_weights_( self._get_occ_range_32_ , select*self._sel_iv_1   , 1 , ef, weights)
        self._update_weights_( self._get_occ_range_32_ , select*self._sel_iv_2   , 2 , ef, weights)
        self._update_weights_( self._get_occ_range_32_ , select*self._sel_iv_3   , 3 , ef, weights)
        self._update_weights_( self._get_occ_range_32_ , select*self._sel_iv_4   , 4 , ef, weights)
##  fourth  case :  E3>Ef>=E2
        select= (ef<self._Etetra[:,:,:,1])*(ef>=self._Etetra[:,:,:,0])
        self._update_weights_( self._get_occ_range_21_ , select*self._sel_iv_1   , 1 , ef, weights)
        self._update_weights_( self._get_occ_range_21_ , select*self._sel_iv_234 , 4 , ef, weights)        

        weights=weights.sum(axis=-1)/24.
        self._average_degen(weights)
        return weights
       

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
    