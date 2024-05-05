import glob
import os
import shutil
import sys

sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import lib

# import SymProFold libraries
path_symprofold = lib.get_main_dir()
sys.path.append(path_symprofold+'lib/')

import ctl
import model_preparation.analyze_prediction_dir


work_dir = os.path.dirname(os.path.realpath(__file__))+'/'
existing_dirs = sorted(glob.glob(work_dir+'*_o*'))

if len(existing_dirs) > 150:
    ctl.e(existing_dirs)
    ctl.e(len(existing_dirs))
    ctl.error('ERROR: more than 150 existing directories with analyzed '+ \
                                                              'predictions')
    sys.exit()

else:            
    for f in existing_dirs:
        if len(f) < 5:
            print('ERROR: path length too short')
            sys.exit()

        try:
            shutil.rmtree(f+'/')
        except IOError:
            pass


prediction_dirs0 = sorted(glob.glob(work_dir+'*_*x*'))
prediction_dirs = []

for d in prediction_dirs0:
    if os.path.isdir(d):
        prediction_dirs.append(d+'/')


prediction_dirs_for_processing = []

for d in prediction_dirs:

    # append directory with unrelaxed files
    prediction_dirs_for_processing.append(d)

    if os.path.isdir(d+d.split('/')[-2]):
        files_ranked = sorted(glob.glob(d+d.split('/')[-2]+ \
                                    '/rank*_*.*.pdb'))

        if len(files_ranked) >= 5:
            # if directory with complete set of relaxed files available,
            # also append directory with relaxed files
            prediction_dirs_for_processing.append(d+d.split('/')[-2]+'/')


for d in prediction_dirs_for_processing:
    model_preparation.analyze_prediction_dir.analyze(d, session)

print('finished')
