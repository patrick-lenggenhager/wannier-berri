#!/usr/bin/env python3


## these linesline if you want to use the git version of the code, instead of the one installed by pip
local_code=True
num_proc=0

import os

if local_code:
   if 'wannierberri' not in os.listdir() :
       os.symlink("../wannierberri","wannierberri")
else:
   if 'wannierberri' in os.listdir() :
       os.remove('wannierberri')

if 'Fe_tb.dat' not in os.listdir():
    os.system('tar -xvf ../data/Fe_tb.dat.tar.gz') 


import wannierberri as wberri

import numpy as np


SYM=wberri.symmetry

omega=np.linspace(0.,1.,101)
system=wberri.System(tb_file='Fe_tb.dat',getAA=True)

generators=[SYM.Inversion,SYM.C4z,SYM.TimeReversal*SYM.C2x]


wberri.integrate(system,
            NK=100,
            omega=omega,
            smearEf=10,
            quantities=["opt_conductivity"],
            numproc=num_proc,
            adpt_num_iter=10,
            fout_name='Fe',
            symmetry_gen=generators,
            restart=False,
            parameters={'mu':12.6283}
            )