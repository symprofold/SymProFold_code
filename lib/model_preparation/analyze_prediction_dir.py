import glob
import os
import sys

import ctl
import bib
import bibfasta
import chimerax_api
import model_preparation.bib
import model_preparation.prepare_file

from structure.modelreg import ModelReg


def analyze(prediction_dir, session):
    '''
    Preprocessing of predicted SymPlex (Symmetric Protein Complexes)
    candidate models
    - alignment of whole model to xy plane
    - counterclockwise renaming of chain ids
    - determination of the symmetry of SymPlex candidate
    '''
    sess = chimerax_api.ChimeraxSession(session, 1)
    model_reg = ModelReg()
    sess.set_model_reg(model_reg)

    sess.run('close session')
    sess.run('set bgColor white')

    path_export_postfix_running = '_running'
                # postfix for temporary name of export folder during runtime

    path_export = prediction_dir[0:-1]+'_o'+path_export_postfix_running+'/'
    files = model_preparation.bib.get_input_files(prediction_dir)

    # get gene id from prediction_dir
    gene_id = prediction_dir.split('/')[-2].split('_')[0]

    if prediction_dir.split('/')[-2] != prediction_dir.split('/')[-3]:
        fasta_path = '/'.join(prediction_dir.split('/')[:-2])+'/'
    else: 
        fasta_path = '/'.join(prediction_dir.split('/')[:-3])+'/'
        prefix = prediction_dir.split('/')[-2]
        path_export = fasta_path+prefix+'_or'+path_export_postfix_running+'/'

    try:
        os.mkdir(path_export)
    except IOError:
        pass


    #fasta_file = sorted(glob.glob(fasta_path+'*.fa*'))
    fasta_file = sorted(glob.glob(fasta_path+gene_id+'.fa*'))

    if len(fasta_file) != 1:
        f = open(prediction_dir+'ERROR_no_unique_fasta.txt', 'w')
        f.write('') 
        f.close()
        ctl.error('no unique fasta')

    seq = bibfasta.get_seq_from_fasta(fasta_file[0])


    # iterate through all model input files.
    for f0 in files:
        model_preparation.prepare_file.prepare_file(f0, seq, \
                    path_export, sess)


    # rename export folder after completion
    os.rename(path_export, \
              path_export.replace(path_export_postfix_running, ''))

    return
