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

#! Copyright message from Quantum Espresso:
#! Copyright (C) 2016 Quantum ESPRESSO Foundation
#! This file is distributed under the terms of the
#! GNU General Public License. See the file `License'
#! in the root directory of the present distribution,
#! or http://www.gnu.org/copyleft/gpl.txt .
#!

## a version that takes into account the previous calculated Ef


import numpy as np
from lazy_property import LazyProperty

###########
## an Auxillary class to complute the change in occupation numbers

class Tetrahedra_tmp():

    def __init__(self,tetra,Ef_new):
    
        self._TETRA_WEIGHTS=tetra._TETRA_WEIGHTS
        self._Ef_old=tetra._Ef_old
        self._Ef_new=Ef_new
        self._select_update = ( (min(self._Ef_new,self._Ef_old)<tetra._Etetra_max) 
                    *  ( max(self._Ef_new,self._Ef_old)>tetra._Etetra_min) ) 
#        print (" Between EF {0} and {1} update {2} bands".format(self._Ef_old,self._Ef_new,self._select_update.sum()/24))
        self._Etetra     = tetra._Etetra     [self._select_update]
        self._E0         = tetra._E0         [self._select_update]
        self._sel_iv_1   = tetra._sel_iv_1   [self._select_update]
        self._sel_iv_2   = tetra._sel_iv_2   [self._select_update]
        self._sel_iv_3   = tetra._sel_iv_3   [self._select_update]
        self._sel_iv_4   = tetra._sel_iv_4   [self._select_update]
        self._sel_iv_123 = tetra._sel_iv_123 [self._select_update]
        self._sel_iv_234 = tetra._sel_iv_234 [self._select_update]

    def _get_e01234(self,select):
        E = self._Etetra[select].T
        return (self._E0[select],E[0],E[1],E[2],E[3],)#+tuple([E for E in self._Etetra[select].T])


    def _get_occ_range_43_(self,ef,e0,e1,e2,e3,e4,iv=None):
        c4 =  (e4 - ef) **2 / ((e4 - e1) * (e4 - e2) * (e4 - e3))
        dosef = 0.3 * c4 *(e1 + e2 + e3 + e4 - 4. * e0 )
        c4*= (e4 - ef)
        if iv<=3: 
            return  1.  - c4 * (e4 - ef) / (e4 - e0) + dosef 
        elif iv==4:
            return  1.  - c4 * (4. - (e4 - ef) * (1. / (e4 - e1) + 1. / (e4 - e2)  + 1. / (e4 - e3) ) ) + dosef 
        else:
            raise RuntimeError("incorrect value of iv={0} in _get_occ_range_43_".format(iv))

    def _get_occ_range_32_(self,ef,e0,e1,e2,e3,e4,iv=None):
        c1 = (ef - e1) **2 / ((e4 - e1) * (e3 - e1))
        c2 = (ef - e1) * (ef - e2) * (e3 - ef)  / ((e4 - e1) * (e3 - e2) * (e3 - e1))
        c3 = (ef - e2) **2 * (e4 - ef) /( (e4 - e2)  * (e3 - e2) * (e4 - e1))
        dosef = 0.1 / (e3 - e1) / (e4 - e1) * (3. * 
                   (e2 - e1) + 6. * (ef - e2) - 3. * (e3 - e1 + e4 - e2) 
                   * (ef - e2) **2 / (e3 - e2) / (e4 - e2) )* (e1 + e2 +  e3 + e4 - 4. * e0 ) 
        if   iv==1:
            return  c1 + (c1 + c2) * (e3 - ef) / (e3 - e1) + (c1 + c2 + c3) * (e4 - ef) / (e4 - e1) + dosef
        elif iv==2:
            return  c1 + c2 + c3 + (c2 + c3)  * (e3 - ef) / (e3 - e2) + c3 * (e4 - ef) / (e4 - e2) + dosef 
        elif iv==3:
            return  (c1 + c2) * (ef - e1) / (e3 - e1) + (c2 + c3) * (ef - e2) / (e3 - e2) + dosef 
        elif iv==4:
            return  (c1 + c2 + c3) * (ef - e1)  / (e4 - e1) + c3 * (ef - e2) / (e4 - e2) + dosef 
        else:
            raise RuntimeError("incorrect value of iv={0} in _get_occ_range_32_".format(iv))

    def _get_occ_range_21_(self,ef,e0,e1,e2,e3,e4,iv=None):
        c4 = (ef - e1) **2 / (e2 - e1) / (e3 - e1) / (e4 - e1)
        dosef = 0.3  *c4* (e1 + e2 + e3 + e4 - 4. * e0 ) 
        c4 *= (ef - e1)
        if   iv==1:
            return   c4 * (4. - (ef - e1) * (1. / (e2 - e1) + 1. / (e3 - e1) + 1. / (e4 - e1) ) )   + dosef 
        elif iv in (2,3,4):
            return   c4 * (ef - e1) / (e0 - e1)  + dosef
        else:
            raise RuntimeError("incorrect value of iv={0} in _get_occ_range_21_".format(iv))
    
    def _update_weights_(self,fun_get_occ,selectIV,iv,ef,weights):
        e0,e1,e2,e3,e4=self._get_e01234(selectIV)
        weights[selectIV]=fun_get_occ(ef,e0,e1,e2,e3,e4,iv)


    def _get_occ(self,Ef):
        weights=np.zeros(self._Etetra.shape[0:2])
        select3= (Ef>=self._Etetra[:,:,3])
        select2= (Ef>=self._Etetra[:,:,2])
        select1= (Ef>=self._Etetra[:,:,1])
        select0= (Ef>=self._Etetra[:,:,0])

        ##  First case :  Ef>=E4
        select=select3
        range4=select.sum()
        weights[select] = 1.
        ##  second case :  E4>Ef>=E3
        select= np.logical_not(select3)*select2
        range34=select.sum()
        self._update_weights_( self._get_occ_range_43_ , select*self._sel_iv_123 , 1 , Ef, weights)
        self._update_weights_( self._get_occ_range_43_ , select*self._sel_iv_4   , 4 , Ef, weights)
        ##  third  case :  E3>Ef>=E2
        select= np.logical_not(select2)*select1
        range23=select.sum()
        self._update_weights_( self._get_occ_range_32_ , select*self._sel_iv_1   , 1 , Ef, weights)
        self._update_weights_( self._get_occ_range_32_ , select*self._sel_iv_2   , 2 , Ef, weights)
        self._update_weights_( self._get_occ_range_32_ , select*self._sel_iv_3   , 3 , Ef, weights)
        self._update_weights_( self._get_occ_range_32_ , select*self._sel_iv_4   , 4 , Ef, weights)
        ##  fourth  case :  E2>Ef>=E1
        select= np.logical_not(select1)*select0
        range12=select.sum()
        self._update_weights_( self._get_occ_range_21_ , select*self._sel_iv_1   , 1 , Ef, weights)
        self._update_weights_( self._get_occ_range_21_ , select*self._sel_iv_234 , 4 , Ef, weights)
        weights=weights.dot(self._TETRA_WEIGHTS)
#        print ("number of tetra in ranges \n 12 : {0}\n 23 : {1}\n 34 : {2}\n 4- : {3}".format(
#                   range12,range23,range34,range4))

        return weights
        
        
    def update_occ(self,occ_old):
        occ_old[self._select_update]+= self._get_occ(self._Ef_new)-self._get_occ(self._Ef_old)  
        
