import ctl
import filesystem
import importlib
import os
import sys


'''
Basic configuration of assemblies.
'''

def assembly_config():
    '''
    Load basic configuration of assemblies.
    '''
    root_path_assemblies = filesystem.clean_path( \
                os.path.dirname(os.path.realpath(__file__))+ \
                '/../../conf/dir_predictions/')

    config = [
    ['Desulfurococcus mucosus', 'E8R795', '', \
            root_path_assemblies+'Dmuc/23/', 1, 0, [], [5], [], 1], \
    ]
    # list for each assembly
    #       0:species name, 1:gene id, 2:symm group
    #       3:main folder of predictions
    #       4:include in combined plot [0:no, 1:yes,part1, 2:yes,part2]
    #       5:large plot for full ensemble w/o interface clustering
    #               [0:no, 1:yes]
    #       6:filter params [lowest ranking score]
    #       7:check params [max number of models per prediction scenario]
    #       8:display params [number of clusters in reduced mode,
    #               [ranks of clusters to exclude]]
    #       9:clustering_mode, 0:no clustering, 1:clustering


    # load additional config in root_path_assemblies, if existing
    sys.path.append(root_path_assemblies)

    if os.path.exists(root_path_assemblies+'config_additional.py'):
        config_additional = \
            importlib.import_module('config_additional', \
                                root_path_assemblies+'config_additional.py')
        config = config_additional.assembly_config()

    return config
