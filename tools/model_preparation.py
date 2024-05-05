import os
import sys
sys.path.append(os.path.dirname(__file__)+'/../lib/')

import ctl
import bib
import bibfasta
import chimerax_api
import filesystem
import model_preparation.bib
import model_preparation.prepare_file

import glob
import pathlib

from structure.modelreg import ModelReg





# Preprocessing of predicted SymPlex (Symmetric Protein Complexes)
# candidate models
# - alignment of whole model to xy plane
# - counterclockwise renaming of chain ids
# - determination of the symmetry of SymPlex candidate
# ================================================================

sess = chimerax_api.ChimeraxSession(session, 0)
model_reg = ModelReg()
sess.set_model_reg(model_reg)

current_model_id = 0
meta = {} 

# use path of starting script
path = str(pathlib.Path().absolute())+'/'

if len(path) < 20:
    # manual path setting
    path = filesystem.clean_path(os.path.dirname(os.path.realpath(__file__))+ \
                '/../conf/dir_predictions/'+'Dmuc/23/E8R795_68x4/')


path_export_postfix_running = '_running'
                # postfix for temporary name of export folder during runtime

path_export = path[0:-1]+'_o'+path_export_postfix_running+'/'
files = model_preparation.bib.get_input_files(path)

# get gene id from path
gene_id = path.split('/')[-2].split('_')[0]

if path.split('/')[-2] != path.split('/')[-3]:
    fasta_path = '/'.join(path.split('/')[:-2])+'/'
else: 
    fasta_path = '/'.join(path.split('/')[:-3])+'/'
    prefix = path.split('/')[-2]
    path_export = fasta_path+prefix+'_or'+path_export_postfix_running+'/'

try:
    os.mkdir(path_export)
except IOError:
    pass


fasta_file = sorted(glob.glob(fasta_path+gene_id+'.fa*'))

if len(fasta_file) != 1:
    f = open(path+'ERROR_no_unique_fasta.txt', 'w')
    f.write('') 
    f.close()
    ctl.error('no unique fasta')

seq = bibfasta.get_seq_from_fasta(fasta_file[0])


# iterate through all model input files.

for f0 in files:
    model_preparation.prepare_file.prepare_file(f0, seq, \
                path_export, sess)


# rename export folder after completion
os.rename(path_export, path_export.replace(path_export_postfix_running, ''))
