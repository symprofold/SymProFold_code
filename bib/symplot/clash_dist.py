import glob

import ctl


'''
Module providing functions for handling clash distributions.
'''

def load(file):
    '''
    Import clash distribution file.
    '''
    files = sorted(glob.glob(file))

    if len(files) != 1:
        return {}

    cl = {}
    f = filesystem.get_file(file)

    for i,l in enumerate(f):
        l2 = l.split(';')

        if len(l2) >= 2:
            cl[float(l2[0])] = int(l2[1])

    return cl


def get_score(cl_dist):
    '''
    Get clash score from distribution file.
    '''
    cl_score = 9 # if no value found

    if 2.0 in cl_dist:
        cl_score = cl_dist[2.0]

    if 2 in cl_dist:
        cl_score = cl_dist[2]

    return cl_score
