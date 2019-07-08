import numpy as np

points=np.array([ [ 0, 0, 0],
                  [ 0, 0, 1],
                  [ 0, 1, 0],
                  [ 0, 1, 1],
                  [ 1, 0, 0],
                  [ 1, 0, 1],
                  [ 1, 1, 0],
                  [ 1, 1, 1] ])

tetra_1box=np.array([ [1,2,3,6],
                      [2,3,4,6],
                      [1,3,5,6],
                      [3,4,6,8],
                      [3,6,7,8],
                      [3,5,6,7] ])-1

tetra_all=[]

for ix in 0,-1:
  for iy in 0,-1:
    for iz in 0,-1:
      shift=np.array([ix,iy,iz],dtype=int)
      for tetra in tetra_1box:
        tet=[tuple(points[ip]+shift) for ip in tetra]
        if tuple( (0,0,0) ) in tet:
           tetra_all.append(np.array([t for t in tet if t!=(0,0,0)]))
           
#print len(tetra_all)
#print tetra_all
points=set(tuple(p) for t in tetra_all for p in t)
print points
print len(points)


       