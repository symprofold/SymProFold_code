import ctl

import os

from bs4 import BeautifulSoup


def get_dssp_from_log(file_path):
    '''
    Extract secondary structure from dssp text output of ChimeraX.
    '''
    with open(file_path, 'r', encoding='utf-8') as f:
        txt = f.read()
    
    txt_object = BeautifulSoup(txt, 'html.parser')
    txt_wohtml = txt_object.get_text()

    return txt_wohtml


def get_dssp(model_id, root_path, sess, additional_paramaters=''):
    '''
    Get secondary structure using dssp in ChimeraX.
    '''
    dssp = sess.run('log clear')
    dssp = sess.run('dssp #'+str(model_id)+' report true '+ \
                                                additional_paramaters)
    dssp = sess.run('log save "'+root_path+'tmplog.htm"')
    dssp_txt = get_dssp_from_log(root_path+'tmplog.htm')
    os.unlink(root_path+'tmplog.htm')

    dssp = dssp_txt.split("\n")

    return dssp
