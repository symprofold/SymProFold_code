import ctl
import glob


'''
Module providing functions to handle input parameters.
'''

def get_multiplicities(path, microorganism_type):
    '''
    Get given set of multiplicities.
    '''
    def_files = sorted(glob.glob(path+'rotsymmorder_*.txt'))

    if len(def_files) == 1:
        filename = def_files[0].split('/')[-1]
        multiplicities_txt = filename.split('rotsymmorder_')[1]. \
                                     split('.txt')[0]
        multiplicities = [int(m) for m in multiplicities_txt]
    elif microorganism_type == 0: # 0: bacteria and archaea, 1: virus
        multiplicities = [2, 3, 4, 6]
    elif microorganism_type == 1:
        multiplicities = [2, 3, 4, 5, 6]
    else:
        ctl.error('get_multiplicities: set of multiplicities unclear')

    return multiplicities
