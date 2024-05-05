import ctl
import numpy as np


'''
Module providing functions for data format conversions.
'''


def matr_to_ranges(matr):
    '''
    Conversion from matrix (first res is 0) to list of ranges
    (first res is 1).
    Each range is represented by a list with 2 entries:
    [first redidue, last residue].
    '''
    dom_ranges = []
    matr_sums = np.sum(matr, axis=1)

    in_range = False
    for i,r in enumerate(matr_sums):

        if r > 0 and not in_range:
            dom_ranges.append([i+1, -1])
            in_range = True
        elif r > 0 and in_range:
            continue
        elif r == 0 and in_range:
            dom_ranges[-1][1] = i
            in_range = False
        elif r == 0 and not in_range:
            continue

    if len(dom_ranges) > 0:
        if dom_ranges[-1][1] == -1:
            dom_ranges[-1][1] = len(matr_sums)

    return dom_ranges


def ranges_to_resid_list(ranges):
    '''
    Conversion from ranges (start with 1) to list or residue ids
    (start with 1).
    '''
    list_ = []

    for range_ in ranges:
        for r in range(range_[0], range_[1]+1):
            list_.append(r)

    return list_


def resids_to_ranges(resids):
    '''
    Conversion from resids (start with 1) to ranges (start with 1).
    '''
    mask = np.zeros((10000, 1))
        #start with 0

    for r in resids:
        mask[r-1][0] = 1

    ranges = matr_to_ranges(mask)

    return ranges
