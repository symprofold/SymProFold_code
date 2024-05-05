import ctl

import glob


'''
Module providing functions for model preparation.
'''

def get_input_files(path):
    '''
    Get input files.
    '''
    files0 = sorted(glob.glob(path+'rank*.pdb'))
    files = []

    for f in files0:
        fn = f.split('/')[-1]

        if len(fn) >= 13:
            files.append(f)

    if len(files) == 0:
        files = sorted(glob.glob(path+'unrelaxed_rank*.pdb'))

        if len(files) == 0:
            files = sorted(glob.glob(path+'*.pdb'))

    return files
