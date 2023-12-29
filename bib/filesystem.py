import ctl

import time
import os


def clean_path(path):
    ''' Clean path string. '''

    path_cleaned = []    
    path_arr = path.split('/')
    for p in path_arr:
        if p == '..':
            path_cleaned.pop()
        else:
            path_cleaned.append(p)

    path_cl = '/'.join(path_cleaned)

    return path_cl


def create_folder(paths):
    '''
    Create folder and subfolders for each entry in a list with path names.
    Input path names can contain filenames.
    '''

    for path in paths:
        path_cl = clean_path(path)

        # get directory path without filename
        path_cl = path_cl.split('/')
        path_cl = '/'.join(path_cl[0:-1])

        if not os.path.exists(path_cl):
            os.makedirs(path_cl)

    return


def get_file(file, attempt=0):
    ''' Open file and create list with lines. '''

    attempt += 1
    ra = 5
 
    try:
        f = open(file, 'r')
        ra = f.readable()
        if ra == True:
            lin = f.readlines()
        else:
            ra = 5
        f.close()

    except IOError:
        ctl.d("ERROR: file:"+file+", attempt:"+str(attempt))
        ctl.d("retry")

        time.sleep(0.1*attempt)
        return get_file(file)

    if ra != True:
        ctl.d(ra) 

    return lin
