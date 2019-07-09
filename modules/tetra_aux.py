import numpy as np


class Transformer():
    def __init__(self,reclattice):
        scale_a=np.max(abs(reclattice[:,0]))
        scale_c=np.max(abs(reclattice[:,2]))
        scale=np.array([scale_a,scale_a,scale_c])
        b1=reclattice[0]/scale
        b2=reclattice[1]/scale
        if b1.dot(b2)<0:b2=-b2
        cross=2*np.cross(b1,b2)/np.linalg.norm(b1)*np.sqrt(2.)
#        print "cross=",cross
        b3=(b1+b2+cross)/3.
        b4=(b1+b2-cross)/3.
#        print "b3,b4=",b3,b4
        B=np.array([b3,b4]).dot(np.linalg.inv(reclattice))
        b3=[b3,b4][np.argmin( np.abs(B-np.round(B)).max(axis=1) )]
#        b3=[b3,b4][np.argmin( np.linalg.norm([b3,b4],axis=1)-np.linalg.norm(b1)) ]
        b123=np.array([b1,b2,b3])*scale[None,:]
        B= np.array(b123.dot(np.linalg.inv(reclattice)).round(),dtype=int)
#        print "reclattice,B,b123",reclattice,B,b123
        assert abs(np.linalg.det(B))==1
        self._trans=np.linalg.inv([[0,1,1],[1,0,1],[1,1,0]]).dot(B)

    def __call__(self,a):
        return np.dot(a,self._trans)


### function to construct tetrahedra for a bcc lattice
def construct_tetra_bcc1(reclattice):
    ###  first find three rec lattice vectors which have positive scalar product
    trans=Transformer(reclattice)
    print ("constructing tetrahedra for bcc1 lattice")
    #let's write the point in coordinates where b1=[0,1,1], b2=[1,0,1], b3=[1,1,0]
    points=np.array([np.roll([s1,s2,0],i) for i in range(3) for s1 in +1,-1 for s2 in +1,-1])
    nump=len(points)
    # now express in the  original reciprocal lattice vectors
    points_rec=trans(points)
    tetrahedra=set(tuple(sorted((i,j,k))) for i in range(nump) for j in range(nump) for k in range(nump)  
        if np.linalg.det( [points[i],points[j],points[k]] )==2  and 
           points[i].dot(points[j])>0  and 
           points[k].dot(points[j])>0  and 
           points[i].dot(points[k])>0 
            )
    tetrahedra=[points_rec[tetra,:] for tetra in tetrahedra]
    assert(len(tetrahedra)==8)
    weights=[1./24]*8
    for i in range(3):
      for s1 in +1,-1:
        for s2 in +1,-1:
          for s3 in +1,-1:
            tetrahedra.append(trans(np.roll([[s1,0,0],[s1,s2,0],[s1,0,s3]],i,axis=1)))
            weights.append(1./48)
    return tetrahedra,weights

def get_bcc_shift(reclattice):
    trans=Transformer(reclattice)
    return trans([1,0,0])

def construct_tetra_bcc2(reclattice):
    ###  first find three rec lattice vectors which have positive scalar product
    trans=Transformer(reclattice)
    print ("constructing tetrahedra for bcc2 lattice")
    #let's write the point in coordinates where b1=[0,1,1], b2=[1,0,1], b3=[1,1,0]
    #we construct 8 tetrahedra
    tetrahedra=[ np.diag([s1,s2,s3]) for s1 in +1,-1  for s2 in +1,-1  for s3 in +1,-1 ]
    # now express in the  original reciprocal lattice vectors
    tetrahedra=[trans(t) for t in tetrahedra]
    return tetrahedra,np.ones(8,dtype=float)/48

### function to construct tetrahedra for a generic lattice
print ("constructing tetrahedra for generic lattice")
def construct_tetra_generic():
    points=np.array([     [ 0, 0, 0],
                          [ 0, 0, 1],
                          [ 0, 1, 0],
                          [ 0, 1, 1],
                          [ 1, 0, 0],
                          [ 1, 0, 1],
                          [ 1, 1, 0],
                          [ 1, 1, 1] ])
        
    tetra_1box=np.array([    [1,2,3,6],
                             [2,3,4,6],
                             [1,3,5,6],
                             [3,4,6,8],
                             [3,6,7,8],
                             [3,5,6,7] ])-1
        
    TETRA_NEIGHBOURS=[]

    for ix in 0,-1:
      for iy in 0,-1:
        for iz in 0,-1:
            shift=np.array([ix,iy,iz],dtype=int)
            for tetra in tetra_1box:
                tet=[tuple(points[ip]+shift) for ip in tetra]
                if tuple( (0,0,0) ) in tet:
                   TETRA_NEIGHBOURS.append(np.array([t for t in tet if t!=(0,0,0)]))
    ntetra=len(TETRA_NEIGHBOURS)
    return  TETRA_NEIGHBOURS,np.ones(ntetra,dtype=float)/ntetra


       