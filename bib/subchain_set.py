import glob
import metadata


'''
Module providing functions to handle subchain sets.
'''


def get_subchain_set(path):
    '''
    Get subchain set defined by fasta files given path.
    '''
    subchain_set = []
    files = sorted(glob.glob(path))

    if len(files) == 0:
        return []
    
    for f in files:
        sc = metadata.get_subchain_abbr(f)
        subchain_set.append(sc)

    return subchain_set
