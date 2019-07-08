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


# returns occupation factor of a band at a k-point in a corner of a tetrahedron
def weights_1band(etetra,Ef):
    # energies will be sorted, remember which is at the corner of interest
    e0=etetea[0]
    ivertex=np.sum(e0>etetra[1:])
    e1,e2,e3,e4=np.sort(etetra)
  
#    weights[Ef>=etetra[3]]=1.
    if Ef>=e4:
        return 1.
    elif Ef>=e3:
        c4 =  (e4 - ef) **3 / ((e4 - e1) * (e4 - e2) * (e4 - e3))
        dosef = 0.3 * (e4 - ef) **2 /((e4 - e1)* (e4 - e2) * (e4 - e3))*(e1 + e2 + e3 + e4 - 4. * e0 )
        if   ivertex in (0,1,2) :
            return 1.  - c4 * (e4 - ef) / (e4 - e0) + dosef 
        elif   ivertex==3: 
            return 1.  - c4 * (4. - (e4 - ef) * (1. / (e4 - e1) + 1. / (e4 - e2) 
                   + 1. / (e4 - e3) ) ) + dosef 
    elif Ef>=e2:
        c1 = (ef - e1) **2 / ((e4 - e1) * (e3 - e1))
        c2 = (ef - e1) * (ef - e2) * (e3 - ef)  / ((e4 - e1) * (e3 - e2) * (e3 - e1))
        c3 = (ef - e2) **2 * (e4 - ef) /( (e4 - e2)  * (e3 - e2) * (e4 - e1))
        dosef = 0.1 / (e3 - e1) / (e4 - e1) * (3. * 
                   (e2 - e1) + 6. * (ef - e2) - 3. * (e3 - e1 + e4 - e2) 
                   * (ef - e2) **2 / (e3 - e2) / (e4 - e2) )* (e1 + e2 +  e3 + e4 - 4. * e0 ) 
        if   ivertex==0:
            return  c1 + (c1 + c2) * (e3 - ef) / (e3 - e1) + (c1 + c2 + c3) * (e4 - ef) / (e4 - e1) + dosef
        elif ivertex==1:
            return  c1 + c2 + c3 + (c2 + c3)  * (e3 - ef) / (e3 - e2) + c3 * (e4 - ef) / (e4 - e2) + dosef 
        elif ivertex==2:
            return  (c1 + c2) * (ef - e1) / (e3 - e1) + (c2 + c3) * (ef - e2) / (e3 - e2) + dosef 
        elif ivertex==3:
            return  (c1 + c2 + c3) * (ef - e1)  / (e4 - e1) + c3 * (ef - e2) / (e4 - e2) + dosef 
    elif Ef>=e1:
        c4 = (ef - e1) **3 / (e2 - e1) / (e3 - e1) / (e4 - e1)
        dosef = 0.3 * (ef - e1) **2 / (e2 - e1) / (e3 - e1) / (e4 - e1) * (e1 + e2 + e3 + e4 - 4. * e0 ) 
        if   ivertex==0:
            return   c4 * (4. - (ef - e1) * (1. / (e2 - e1) + 1. / (e3 - e1) + 1. / (e4 - e1) ) )   + dosef 
        elif ivertex in (1,2,3):
            return   c4 * (ef - e1) / (e0 - e1)  + dosef     
    else:
        return  0


def average_degen(E,weights):
    # make sure that degenerate bands have same weights
    borders=np.hstack( ( [0], np.where( (E[1:]-E[:-1])>1e-5)[0]+1, [len(E)]) )
    degengroups=[ (b1,b2) for b1,b2 in zip(borders,borders[1:]) if b2-b1>1]
    for b1,b2 in degengroups:
       weights[b1:b2]=weights[b1:b2].mean()

def weights_all_bands_1tetra(Etetra,Ef):
    weights=np.array([weights_1band(etetra,ef) for etetra in Etetra])

def occ_factors(E,E_neigh,Ef):
#  E_neigh is a dict (i,j,k):E 
# where i,j,k = -1,0,+1 - are coordinates of a k-point, relative to the reference point
    num_wann=E.shape[0]
    assert( num_wann.shape[:3]==(3,3,3))
    occ=np.zeros(num_wann)
    weights=np.zeros(num_wann)
    Etetra=np.zeros( (num_wann,4),dtype=float)
    Etetra[:,0]=E
    for tetra in __TETRA_NEIGHBOURS:
       for i,p in enumerate(tetra):
          Etetra[:,i+1]=E_neigh[tuple(p)]
          weights+=weights_all_bands_1tetra(Etetra,Ef)
    
    average_degen(E,weights)
    return weights

    