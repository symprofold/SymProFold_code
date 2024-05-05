import glob
import os
import sys

import ctl
import filesystem


'''
Module providing general libraries for the program flow.
'''

def get_preassemblies_dir(conf, execution_file):
    '''
    Determine preassemblies directory for a given species by evaluating
    the path information files (preassemblies_path*.txt).
    '''

    # Determine path of 'path information file' (preassemblies_path*.txt)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    path_folders = os.path.dirname(execution_file).split('/')

    for i in range(0, len(path_folders)):
        path_folders_pathfile = path_folders[0:len(path_folders)-i]
        folder_pathfile = '/'.join(path_folders_pathfile)+'/'
        files = sorted(glob.glob(folder_pathfile+'preassemblies_path*.txt'))

        if len(files) > 0:
            for f in files:
                # read relative path information from 'path information file'
                path = filesystem.get_file(f)[0].strip()

                path = filesystem.clean_path(folder_pathfile+path)

                if search_species_dir(path, conf):
                    return path

    ctl.error('get_preassemblies_dir: '+ \
                                'path to preassemblies directory not found')

    return False


def search_species_dir(directory, conf):
    '''
    Search species directory in preassemblies directory.
    '''
    subdirs = sorted(glob.glob(directory+'*/'))

    for subdir in subdirs:
        species_dir = subdir.split('/')[-2]

        if species_dir == conf.species:
            return True

    return False
