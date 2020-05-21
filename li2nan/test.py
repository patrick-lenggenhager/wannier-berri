#!/usr/bin/env python3

num_proc=4

# load numpy
import numpy as np

# load wannierberri
import os
if 'wannierberri' not in os.listdir() :
       os.symlink("../wannierberri","wannierberri")
import wannierberri as wberri

# define paramters
#omega = np.linspace(0.,1,100)
Efermi = np.linspace(-2,2,100)

# create system
system = wberri.System(tb_file='Li2NaN_strain=0_tb.dat',getAA=True)

# symmetries
SYM = wberri.symmetry
generators = [SYM.C6z, SYM.Mx*SYM.Inversion, SYM.TimeReversal, SYM.Inversion]


wberri.integrate(system,
    NK=100,
    Efermi=Efermi,
    #omega=omega
    smearEf=10,
    quantities=["dos","cumdos"],
    numproc=num_proc,
    adpt_num_iter=10,
    fout_name='Li2NaN_strain=0',
    symmetry_gen=generators,
    restart=False,
)
