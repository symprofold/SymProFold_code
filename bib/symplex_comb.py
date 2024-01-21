import ctl
import metadata

import glob


'''
Module providing functions to determine combinations of SymPlexes.
'''

def symplex_comb(ax):
    '''
    Get combinations of SymPlexes.
    '''
    comb = []

    for i,f0 in enumerate(ax[0].files):
        for j,f1 in enumerate(ax[1].files):
            files0 = sorted(glob.glob(ax[0].path+'**/'+f0))
            files1 = sorted(glob.glob(ax[1].path+'**/'+f1))

            if len(files0) == 1 and len(files1) == 1:
                comb.append((i, j))
            else:
                ctl.p('not possible:')
                ctl.p(files0)
                ctl.p(files1)

    return comb


def get_symm_order(ax):
    '''
    Get order of rotational symmetry in a complex model from path.
    '''
    path = ax.model_active_path
    fn = path.split('/')[-1]

    files = sorted(glob.glob(ax.path+'**/'+fn))

    if len(files) != 1:
        return -1
    else:
        path = files[0]

    dir_name = path.split('/')[-2]
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


def overlapping_domains(symplex0_folder, symplex1_folder, conf):
    '''
    Get list of overlapping domains.
    '''
    ov_domains = []

    sc0 = metadata.get_subchain_abbr(symplex0_folder)
    sc1 = metadata.get_subchain_abbr(symplex1_folder)

    if sc0 == 'FL':
        sc0 = (1, len(conf.domains))
    else:
        sc0 = metadata.subchainabbr_to_subchains(sc0)

    if sc1 == 'FL':
        sc1 = (1, len(conf.domains))
    else:
        sc1 = metadata.subchainabbr_to_subchains(sc1)

    if len(sc0) == 1:
        sc0 = (sc0[0], sc0[0])

    if len(sc1) == 1:
        sc1 = (sc1[0], sc1[0])

    for d0 in range(sc0[0], sc0[1]+1):
        for d1 in range(sc1[0], sc1[1]+1):
            if d0 == d1:
                if d0 not in ov_domains:
                    ov_domains.append(d0)

    return ov_domains


def seq_coverage(sc0, sc1, conf):
    '''
    Check if SymPlex combination covers full sequence.
    '''

    # check if N terminal end is included in sequence
    if not(sc0 == 'FL' or sc1 == 'FL' or \
           int(sc0[0]) == 1 or int(sc1[0]) == 1):
        return False

    # check if C terminal end is included in sequence
    if not(sc0 == 'FL' or sc1 == 'FL' or \
           int(sc0[-1]) == len(conf.domains) or \
           int(sc1[-1]) == len(conf.domains)):
        return False

    return True


def subchain_order(symplex0_folder, symplex1_folder, conf):
    '''
    Get upstream/downstream order of subchains.

    Return:
        int: SymPlex index with subchain that is located more upstream.
    '''
    sc0 = metadata.get_subchain_abbr(symplex0_folder)
    sc1 = metadata.get_subchain_abbr(symplex1_folder)

    if sc0 == 'FL':
        sc0 = (1, len(conf.domains))
    else:
        sc0 = metadata.subchainabbr_to_subchains(sc0)

    if sc1 == 'FL':
        sc1 = (1, len(conf.domains))
    else:
        sc1 = metadata.subchainabbr_to_subchains(sc1)

    if len(sc0) == 1:
        sc0 = (sc0[0], sc0[0])

    if len(sc1) == 1:
        sc1 = (sc1[0], sc1[0])

    if sc0[0] < sc1[0] or sc0[1] < sc1[1]:
        return 0

    if sc0[0] > sc1[0] or sc0[1] > sc1[1]:
        return 1

    if sc0[0] == sc1[0] or sc0[1] == sc1[1]:
        return -1

    ctl.error('subchain_order')

    return
