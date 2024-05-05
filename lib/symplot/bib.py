import os

'''
Module providing general functions for the creation of SymPlots.
'''


def get_microorganism_type(path):
    '''
    Get type of microorganism (bacteria, virus).
    '''
    microorganism_type = 0

    if 'virus' in path:
        microorganism_type = 1
        # 0: bacteria, 1: virus

    return microorganism_type


def get_plottype(path):
    '''
    Determine type of plot.
    Types of plot: specific species, general plot with many species
    '''
    if os.path.exists(path+'plottype_full.txt'):
        plottype = 1 # general plot with many species
    else:
        plottype = 0 # specific species

    return plottype
