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

__points=np.array([ [ 0, 0, 0],
                  [ 0, 0, 1],
                  [ 0, 1, 0],
                  [ 0, 1, 1],
                  [ 1, 0, 0],
                  [ 1, 0, 1],
                  [ 1, 1, 0],
                  [ 1, 1, 1] ])

__tetra_1box=np.array([ [1,2,3,6],
                      [2,3,4,6],
                      [1,3,5,6],
                      [3,4,6,8],
                      [3,6,7,8],
                      [3,5,6,7] ])-1

__TETRA_NEIGHBOURS=[]

for ix in 0,-1:
  for iy in 0,-1:
    for iz in 0,-1:
      shift=np.array([ix,iy,iz],dtype=int)
      for tetra in __tetra_1box:
        tet=[tuple(__points[ip]+shift) for ip in tetra]
        if tuple( (0,0,0) ) in tet:
           __TETRA_NEIGHBOURS.append(np.array([t for t in tet if t!=(0,0,0)]))
           
NEIGHBOURS=set(tuple(p) for t in __TETRA_NEIGHBOURS for p in t)

def Ntetra():  # should be 24
    return len(__TETRA_NEIGHBOURS)

def weights_all(Etetra,ef):
#    Etetra is an array NKFFT x nband x Ntetra x 4
    # energies will be sorted, remember which is at the corner of interest
    # do it for one Ef here
    assert Etetra.shape[-1]==4
    assert Etetra.shape[-2]==Ntetra()
    
    weights=np.zeros(Etetra.shape[:-1])
    
    
    ivertex=Etetra.argsort(axis=-1).argmin(axis=-1)
    E0=Etetra[:,:,:,0]
    Etetra=np.sort(Etetra,axis=-1)

##  First case :  Ef>=E4
    
    select= ef>=Etetra[:,:,:,3]
    weights[select]+=1.

##  second case :  E4>Ef>=E3
    select= (ef<Etetra[:,:,:,3])*(ef>=Etetra[:,:,:,2])
    E1234=Etetra[select]
    E0_sel=E0[select]
    iV=ivertex[select]
    weights_select=np.zeros(E1234.shape[0],dtype=float)

    sel_iv= (iV<=2)
    e1=E1234[sel_iv][:,0]
    e2=E1234[sel_iv][:,1]
    e3=E1234[sel_iv][:,2]
    e4=E1234[sel_iv][:,3]
    e0=E0_sel[sel_iv]
    c4 =  (e4 - ef) **3 / ((e4 - e1) * (e4 - e2) * (e4 - e3))
    dosef = 0.3 * (e4 - ef) **2 /((e4 - e1)* (e4 - e2) * (e4 - e3))*(e1 + e2 + e3 + e4 - 4. * e0 )
    weights_select[sel_iv]  = 1.  - c4 * (e4 - ef) / (e4 - e0) + dosef 


    sel_iv= (iV==3)
    e1=E1234[sel_iv][:,0]
    e2=E1234[sel_iv][:,1]
    e3=E1234[sel_iv][:,2]
    e4=E1234[sel_iv][:,3]
    e0=E0_sel[sel_iv]
    c4 =  (e4 - ef) **3 / ((e4 - e1) * (e4 - e2) * (e4 - e3))
    dosef = 0.3 * (e4 - ef) **2 /((e4 - e1)* (e4 - e2) * (e4 - e3))*(e1 + e2 + e3 + e4 - 4. * e0 )
    weights_select[sel_iv]  = 1.  - c4 * (4. - (e4 - ef) * (1. / (e4 - e1) + 1. / (e4 - e2)  + 1. / (e4 - e3) ) ) + dosef 
    
    weights[select]+=weights_select
##  Third and Forth case to be implemented
   
    return weights.sum(axis=-1)/24.
   
    

    


def average_degen(Eall,weights):
    # make sure that degenerate bands have same weights
    for ik,E in enumerate(Eall):
        borders=np.hstack( ( [0], np.where( (E[1:]-E[:-1])>1e-5)[0]+1, [len(E)]) )
        degengroups=[ (b1,b2) for b1,b2 in zip(borders,borders[1:]) if b2-b1>1]
        for b1,b2 in degengroups:
           weights[ik,b1:b2]=weights[ik,b1:b2].mean()

def weights_all_bands_1tetra(Etetra,Ef):
    weights=np.array([weights_1band(etetra,ef) for etetra in Etetra])

def get_occ(E,E_neigh,Ef):
#  E_neigh is a dict (i,j,k):E 
# where i,j,k = -1,0,+1 - are coordinates of a k-point, relative to the reference point
    ntetra=Ntetra()
    Etetra=np.zeros( E.shape+(ntetra,4),dtype=float)
    for it,tetra in enumerate(__TETRA_NEIGHBOURS):
        Etetra[:,:,it,0]=E
        for j in range(3):
            Etetra[:,:,it,j]=E_neigh[tuple(tetra[j])]
    
    weights=weights_all(Etetra,Ef)
    average_degen(E,weights)
    return weights

    