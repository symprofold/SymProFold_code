import glob

import ctl


'''
Module providing functions for loading and handling prediction scenarios.
Prediction scenario: configuration of a oligomer prediction:
specific species/gene, subchain, multiplicity.
'''

def get_folders(path):
    '''
    Get folders of all prediction scenarios (one folder per scenario) under
    a given path.
    '''
    folders = sorted(glob.glob(path+'*_o'))

    # remove non-regular folder names, e.g. with '_m'
    folders = clean_folders(folders)


    for i,folder in enumerate(folders):
        or_search = sorted(glob.glob(folder+'r'))

        if or_search != []:
            folders[i] = or_search[0]

    return folders


def clean_folders(folders):
    '''
    Remove non-regular folder names, e.g. with '_m'.
    '''
    folders_cleaned = []
    to_exclude = ['_m', '__o']

    for f in folders:
        excluded = False
        
        for e in to_exclude:
            if e in f:
                excluded = True
                break

        if not excluded:
            folders_cleaned.append(f)

    return folders_cleaned


def get_path(file_in_path_predictionscenario, root_path):
    '''
    Get main path of prediction scenario from a coord filename.
    '''
    root_path_depth = len(root_path.split('/'))
    path_prediction_set = '/'.join(file_in_path_predictionscenario.split('/') \
                                                        [:root_path_depth])+'/'

    return path_prediction_set


def get_coord_filename(metadata_path):
    '''
    Get filename of the coord file that belongs to a given metadata path of a
    prediction scenario.
    '''
    metadata_filename = metadata_path.split('/')[-1]

    if '_cla' in metadata_filename:
        part1 = metadata_filename.split('_cla')[0]

    coord_filename = part1+'.pdb'

    return coord_filename
