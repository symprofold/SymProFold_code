import glob

import ctl


'''
Module providing functions to retrieve and handle metadata of coord files.

This includes coding information using terms/abbreviations or extracting
information (e.g. terms/abbreviations in directory/file names).
'''

def get_subchain_abbr(path):
    '''
    Get subchain abbreviation of a complex model from path.
    '''
    path_name = path.split('/')[-1]
    if path_name == '':
        path_name = path.split('/')[-2]
    if path_name == '':
        ctl.error('get_subchain_abbr')

    sc = path_name.split('_')[1].split('x')[0]
    if sc == '':
        sc = 'FL'

    return sc


def get_symm_order(path):
    '''
    Get order of rotational symmetry in a complex model from path.
    '''
    dir_name = path.split('/')[-3]
    mult = -1

    if 'symm_180' in dir_name:
        mult = 2
    elif 'symm_120' in dir_name:
        mult = 3       
    elif 'symm_090' in dir_name:
        mult = 4
    elif 'symm_072' in dir_name:
        mult = 5
    elif 'symm_060' in dir_name:
        mult = 6
    elif 'symm_051' in dir_name:
        mult = 7
    elif 'symm_045' in dir_name:
        mult = 8

    return mult


def get_prediction_scenario(path):
    '''
    Get prediction scenario name of a complex model from path.
    '''
    prediction_scenario = path.split('_o')[0].split('/')[-1]

    return prediction_scenario


def get_rank(path):
    '''
    Get rank of a complex model from path.
    '''
    filename = path.split('/')[-1]

    if not '.pdb' in filename:
        return -1

    rank_str = filename.split('rank')[1][:2]
    rank = int(rank_str)

    return rank


def get_score(path):
    '''
    Get score of a complex model from path.
    '''
    filename = path.split('/')[-1]
    score = -1

    if not '.pdb' in filename:
        return -1

    if not '_cla' in filename:
        score_str = filename.split('.pdb')[-2]
        score_str = score_str.split('_')[-1]
    else:
        score_str = filename.split('_cla')[0]
        score_str = score_str.split('_')[-1]

    try:
        score = float(score_str)
    except ValueError: 
        return -1

    if score%1 == 0: # if score is a full integer number, e.g. model number
        return -1

    return score


def subchainabbr_to_subchains(subchain_abbr):
    '''
    Get subchain range (first domain, last domain) from subchain abbreviation.

    Args:
        subchain_abbr, e.g. "c24"
        
    Returns:
        subchain_range, e.g. (2, 4)
    '''
    subchain_range = ()

    # remove prefix letter (e.g. 'c' in 'c24' etc.)
    if len(subchain_abbr) > 1:
        if subchain_abbr[0] in 'abcdefghijklmnopqrstuvwxyz':
            subchain_abbr = subchain_abbr[1:]

    if len(subchain_abbr) == 0:
        pass
    elif len(subchain_abbr) == 1: # e.g. '2' -> (2,)
        subchain_range = (int(subchain_abbr[0]),)
    elif len(subchain_abbr) == 2: # e.g. '24' -> (2, 4)
        subchain_range = (int(subchain_abbr[0]), int(subchain_abbr[1]))
    elif len(subchain_abbr) == 3: # e.g. '212' -> (2, 12)
        subchain_range = (int(subchain_abbr[0]), int(subchain_abbr[1:3]))
    elif len(subchain_abbr) == 4: # e.g. '1012' -> (10, 12)
        subchain_range = (int(subchain_abbr[0:2]), int(subchain_abbr[2:4]))
    else:
        ctl.e('subchains_from_subchaindef')

    return subchain_range


def get_chainpair_nnptm(path):
    '''
    Get nnpTM score of a chain pair from path.
    '''
    filename = fn.split('/')[-1]
    nnptm = float(filename.split('_iptm')[1].split('.txt')[0])

    return nnptm


def get_clashes(path):
    '''
    Get clashes per 100aa of a complex model from path.
    '''
    filename = path.split('/')[-1]
    cla = filename.split('_cla')[1].split('.pdb')[0]
    cla = float(cla.replace('-', '.'))

    if cla < 0:
        ctl.e(filename)
        ctl.e(cla)
        ctl.error('ERROR: get_clashes')

    return cla


def get_intermol_betasheets(path):
    '''
    Get fraction of beta sheet residues that is involved in
    intermolecular sheets from path.
    '''
    filename = path.split('/')[-1]

    path_dir = '/'.join(path.split('/')[:-1])+'/../sheetintermol/'
    files_fraction = sorted(glob.glob( \
                                path_dir+filename.split('_cla')[0]+'*.txt'))

    if len(files_fraction) != 1:
        ctl.e(path_dir+filename.split('_cla')[0]+'*.txt')
        ctl.e(path_dir)
        ctl.e(files_fraction)
        ctl.error('get_intermol_betasheets: len(files_fraction) != 1')

    fraction_minlen = files_fraction[0].split('_sheetintermol')[1]. \
                          split('_')[0]
    fraction = files_fraction[0].split('_sheetintermol')[1].split('_')[1]
    betasheet_res_n = int(files_fraction[0].split('_sheetintermol')[1]. \
                          split('_')[2].split('.txt')[0])
    fraction_minlen = float(fraction_minlen)
    fraction = float(fraction)

    if fraction < 0:
        ctl.e(filename)
        ctl.e(fraction)
        ctl.error('get_intermol_betasheets: fraction < 0')

    if betasheet_res_n < 0:
        ctl.e(filename)
        ctl.e(betasheet_res_n)
        ctl.error('get_intermol_betasheets: betasheet_res_n < 0')

    return fraction_minlen, fraction, betasheet_res_n


def get_roll_clashes(path):
    '''
    Get clashes per 100aa (rolling range of 200aa) of a complex model from
    path.
    '''
    filename = path.split('/')[-1]
    path_dir = '/'.join(path.split('/')[:-1])+'/'
    files_roll = sorted(glob.glob(path_dir+filename.split('_cla')[0]+'*.txt'))

    if len(files_roll) != 1:
        ctl.error('get_roll_clashes')
        
    filename = files_roll[0].split('/')[-1]
    cla = filename.split('_rollcla')[1].split('_')[0]
    cla = float(cla.replace('-', '.'))

    if cla < 0:
        ctl.e(filename)
        ctl.e(cla)
        ctl.error('get_roll_clashes')

    return cla


def get_rotsymm_ang(file_path):
    '''
    Get rotational symmetry angle of a oligomer model from path.
    '''
    filename = file_path.split('/')[-1]
    folder = '/'.join(file_path.split('/')[:-2])+'/'

    files_rotang = sorted(glob.glob(folder+filename.split('_cla')[0]+ \
                                        '_rotang*.txt'))

    if 'symm_unclassified' in file_path or 'no_symm_found' in file_path:
        return []

    if len(files_rotang) != 1:
        ctl.e(file_path)
        ctl.e(folder+filename.split('_cla')[0]+'_rotang*.txt')
        ctl.e(files_rotang)
        ctl.error('get_rotang')

    if not '_rotang' in files_rotang[0]:
        ctl.e(files_rotang)
        ctl.error('get_rotsymm_ang: rotang metadata not found in path')

    filename = files_rotang[0].split('/')[-1]
    rotang_txt = filename.split('_rotang')[1].split('.txt')[0]

    if '_' in rotang_txt:
        rotangs = [float(rotang_txt.split('_')[0].replace('-', '.')), \
                   float(rotang_txt.split('_')[1].replace('-', '.'))]
    else:
        rotangs = [int(rotang_txt.split('-')[0]), \
                   int(rotang_txt.split('-')[1])]

    if min(rotangs) < 10:
        ctl.e(rotangs)
        ctl.error('get_rotsymm_ang: min(rotangs) < 10')

    return rotangs
