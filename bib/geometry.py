import bibpdb
import ctl

import math
import os


def dist(x1, x2):
    ''' Get distance of 2 input vectors. '''

    if len(x1) != 3 or len(x1) != 3:
        ctl.e('dist_vect: input vector does not have dimension 3')

    d = math.sqrt( (x2[0]-x1[0])**2 + \
                   (x2[1]-x1[1])**2 + \
                   (x2[2]-x1[2])**2 )

    return d


def dist_xy(x1, x2):
    ''' Get xy distance of 2 input vectors. '''

    if len(x1) != 3 or len(x1) != 3:
        ctl.e('dist_vect: input vector does not have dimension 3')

    d = math.sqrt( (x2[0]-x1[0])**2 + \
                   (x2[1]-x1[1])**2 )

    return d


def crossproduct(a, b):
    ''' Calculate crossproduct. '''

    if len(a) != 3 or len(b) != 3:
        ctl.e('crossproduct: input vector does not have dimension 3')

    crossp = [a[2-1]*b[3-1]-a[3-1]*b[2-1], \
              a[3-1]*b[1-1]-a[1-1]*b[3-1], \
              a[1-1]*b[2-1]-a[2-1]*b[1-1]]

    return crossp


def dotproduct(a, b):
    ''' Calculate dot product. '''

    if len(a) != 3 or len(b) != 3:
        ctl.e('dotproduct: input vector does not have dimension 3')

    dotp = a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

    return dotp


def norm(a):
    '''
    Normalize vector.
    '''
    n = 1/dist(a, [0, 0, 0])
    a_norm = [ai*n for ai in a]

    return a_norm


def rotation_around_z(vector, angle):
    ''' Rotate given vector around z axis. '''

    rotated = [0, 0]

    if len(vector) == 3:
        rotated = [0, 0, 0]

    angle = angle/180*math.pi

    rotated[0] = round(vector[0]*math.cos(angle)-vector[1]*math.sin(angle), 5)
    rotated[1] = round(vector[0]*math.sin(angle)+vector[1]*math.cos(angle), 5)

    if len(vector) == 3:
        rotated[2] = vector[2]

    return rotated


def get_center(coord, chainid0='A', res_range=0):
    ''' Get center of CA atoms of given coordinates. '''

    center = [0,0,0]
    sum0 = [0,0,0]
    sumN = 0

    if res_range == 0:
        res_range = [1, 10000]

    for l in coord:
        if l[0:4] == 'ATOM' or l[0:6] == 'HETATM':
           name = bibpdb.get_f(l, 'name')
           chainid = bibpdb.get_f(l, 'chainid')
           resnr = bibpdb.get_f(l, 'resnr')

        if resnr < res_range[0] or resnr > res_range[1]:
            continue

        if (chainid == chainid0 or chainid0 == 'all') and name == 'CA':
            sumN += 1
            sum0[0] += bibpdb.get_f(l, 'x')
            sum0[1] += bibpdb.get_f(l, 'y')
            sum0[2] += bibpdb.get_f(l, 'z')     

    center[0] = sum0[0]/sumN
    center[1] = sum0[1]/sumN
    center[2] = sum0[2]/sumN
    
    return center


def get_extrema(coord, axis, chainid0='A', res_range=0):
    '''
    Get residue id with the lowest and the highest coordinate (CA atom) on
    given axis ("x", "y", "z"). Source is the content of a pdb coordinate file
    as line list.
    '''

    extrema = [[-1, 1000000], [-1, -1]]
        # [[residue id, lowest coordinate], [residue id, highest coordinate]]

    if res_range == 0:
        res_range = [1, 10000]

    for l in coord:
        if l[0:4] == 'ATOM' or l[0:6] == 'HETATM':
           name = bibpdb.get_f(l, 'name')
           chainid = bibpdb.get_f(l, 'chainid')
           resnr = bibpdb.get_f(l, 'resnr')

        if resnr < res_range[0] or resnr > res_range[1]:
            continue

        if (chainid == chainid0 or chainid0 == 'all') and name == 'CA':
            if bibpdb.get_f(l, 'z') < extrema[0][0]:
                extrema[0] = [bibpdb.get_f(l, 'z'), resnr]

            if bibpdb.get_f(l, 'z') > extrema[1][0]:
                extrema[1] = [bibpdb.get_f(l, 'z'), resnr]
    
    return extrema
