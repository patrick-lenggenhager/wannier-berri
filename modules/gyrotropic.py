import numpy as np
from collections import Iterable

from berry import fac_ahc, fac_morb, calcImf_K, calcImfgh_K
from spin import calcSpin


def calcAHC(data,Efermi,degen_thresh=None):
    def function(AHC,Efermi,data,degen,E_K_av,delE_K_av,ik):
        imf=calcImf_K(data,degen,ik )
        for e,f,ndeg in zip(E_K_av,imf,degen):
            AHC[Efermi>e]+=f*ndeg[2]
    return __calcSmth_band(data,Efermi,(3,) ,function,   degen_thresh=degen_thresh)*fac_ahc/(data.NKFFT_tot*data.cell_volume)


def calcMorb(data,Efermi,degen_thresh=None):
    def function(Morb,Efermi,data,degen,E_K_av,delE_K_av,ik):
        imf,img,imh=calcImfgh_K(data,degen,ik )
        for e,f,g,h,ndeg in zip(E_K_av,imf,img,imh,degen):
            sel= Efermi>e
            Morb[sel]+=((g+h)-2*f[None,:]*Efermi[sel,None])*ndeg[2]
    return __calcSmth_band(data,Efermi,(3,) ,function,   degen_thresh=degen_thresh)*fac_morb/(data.NKFFT_tot)


def calcGyrotropic(data,Efermi,degen_thresh=None):
    def function(GYRO,Efermi,data,degen,E_K_av,delE_K_av,ik):
        imf,img,imh=calcImfgh_K(data,degen,ik )
        imgh=img-imh
        spin=calcSpin(data,degen)
        for e,de,f,gh,s,ndeg in zip(E_K_av,delE_K_av,imf,imgh,spin,degen):
            sel= Efermi>e
#            print (de.shape,e.shape,GYRO.shape)
            GYRO[sel,0]+=de[:,None]*de[None,:]*ndeg[2]
            GYRO[sel,1]+=de[:,None]*f [None,:]*ndeg[2]
            GYRO[sel,2]+=de[:,None]*gh[None,:]*ndeg[2]
            GYRO[sel,3]+= f[:,None]*gh[None,:]*ndeg[2]
            GYRO[sel,4]+= f[:,None]*s [None,:]*ndeg[2]
            GYRO[sel,5]+=de[:,None]*s [None,:]*ndeg[2]
    return __calcSmth_band(data,Efermi,(6,3,3) ,function,   degen_thresh=degen_thresh)*fac_morb/(data.NKFFT_tot)




##  a general procedure to evaluate smth band-by-band in the Fermi sea
def __calcSmth_band(data,Efermi,shape,function, degen_thresh=None):
    if degen_thresh is None:
        degen_thresh=-1
    degen_bands,E_K_av,delE_K_av=__get_degen_bands(data.E_K,degen_thresh,data.delE_K)

    if not(isinstance(Efermi, Iterable)): 
        Efermi=np.array([Efermi])

    RES=np.zeros( (len(Efermi),)+tuple(shape) )
    for ik in range(data.NKFFT_tot) :
        function(RES,Efermi,data,degen_bands[ik],E_K_av[ik],delE_K_av[ik],ik)
    return RES 


def __get_degen_bands(E_K,degen_thresh,delE_K=None):
    A=[np.hstack( ([0],np.where(E[1:]-E[:1]>degen_thresh)[0]+1, [E.shape[0]]) ) for E in E_K ]
    deg= [[(ib1,ib2,ib2-ib1) for ib1,ib2 in zip(a,a[1:])] for a in A]
    Eav= [ np.array( [E[b1:b2].mean() for b1,b2,nd in deg  ]) for E,deg in zip(E_K,deg)]
    if delE_K is None:
        print ("delE_K is None")
        return deg,Eav
    else:
        print ("delE_K has shape {0}".format(delE_K.shape))
        delEav= [ np.array( [delE[b1:b2].mean(axis=0) for b1,b2,nd in deg  ]) for delE,deg in zip(delE_K,deg) ]
        return deg,Eav,delEav

def __average_degen(X,degen):
    np.array( [X[b1:b2].mean(axis=0) for b1,b2,nd in degen  ]) 