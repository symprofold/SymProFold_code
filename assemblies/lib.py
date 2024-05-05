import glob
import os
import sys


'''
Module providing general libraries for the program flow.
'''

def get_main_dir():
    '''
    Determine main directory of the SymProFold installation by evaluating
    the path information file (symprofold_path.txt).
    '''

    # Determine path of 'path information file' (symprofold_path.txt)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    path_folders = os.path.dirname(os.path.realpath(__file__)).split('/')
    path_pathfile = None

    for i in range(0, len(path_folders)):
        path_folders_pathfile = path_folders[0:len(path_folders)-i]
        folder_pathfile = '/'.join(path_folders_pathfile)+'/'
        files = sorted(glob.glob(folder_pathfile+'symprofold_path.txt'))

        if len(files) == 1:
            path_pathfile = files[0]
            break

    if path_pathfile == None:
        print('path to main directory \'SymProFold/\' of SymProFold '+ \
              'installation not found')
        sys.exit()


    # Read relative path information from 'path information file'
    # (symprofold_path.txt)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    with open(path_pathfile, 'r', encoding='utf-8') as file:
        file_txt = file.read()

    rel_path = file_txt.split("\n")[0].strip()


    # Combine folder_pathfile and rel_path to main_path
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    main_path = folder_pathfile+rel_path


    # Clean path string
    # ~~~~~~~~~~~~~~~~~
    path_cleaned = []    
    path_arr = main_path.split('/')
    for p in path_arr:
        if p == '..':
            path_cleaned.pop()
        else:
            path_cleaned.append(p)

    main_path = '/'.join(path_cleaned)

    return main_path
