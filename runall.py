import os
import sys
sys.path.append(os.path.dirname(__file__)+'/bib/')

import glob
import ctl


""" Run all runscripts in folder "runscripts" """

run_level = 5 # 5:full

path = os.path.dirname(__file__)+'/runscripts/'

files = sorted(glob.glob(path+'*_run*.py'))
ctl.p(files)

for f in files:
    ctl.p(f)
    exec(open(f).read())

ctl.p('finished')
