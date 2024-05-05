import glob

import ctl
import metadata
import symplot.prediction_scenario

'''
Module providing functions to handle nnpTM scores.
'''

def get_nnptm(path_file, root_path):
    '''
    Get nnpTM score of a complex model from path.
    '''
    nnptms = get_nnptms(path_file, root_path)

    if len(nnptms) > 0:
        nnptm = max(nnptms)
    else:
        nnptm = -1

    return nnptm


def get_nnptms(path_file, root_path):
    '''
    Get list of nnpTM scores of next neighbor pairs.
    '''
    filename = symplot.prediction_scenario.get_coord_filename(path_file)

    if not '.pdb' in filename:
        return -1

    filename_part1 = '_'.join(filename.split('_')[:-1])

    path_prediction_set = symplot.prediction_scenario. \
                                                get_path(path_file, root_path)

    files = sorted(glob.glob(path_prediction_set+filename_part1+ \
                                '_ch*_iptm*.txt'))
    nnptm_vals = []

    for f in files:

        # get nnpTM score of a next neighbor pair
        nnptm = metadata.get_chainpair_nnptm(f)
        nnptm_vals.append(nnptm)

    return nnptm_vals
