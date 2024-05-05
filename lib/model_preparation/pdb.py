import bib
import ctl


'''
Module providing functions to handle pdb files during model preparation.
'''

def align_to_fasta(input_file, output_file, seq, sess):
    ''' Align residue ids of pdb file to given sequence. '''

    sess.run('close session')
    sess.run('set bgColor white')

    current_model_id = 0
    meta = {}
    current_model_id, meta = bib.open_model(sess, \
                                            input_file, current_model_id, \
                                            meta, 1)
    sess.run('split #'+str(current_model_id))
    
    seq_model = sess.get_seq(str(current_model_id)+'.1')
    offset = seq.find(seq_model)

    if offset == -1:
        ctl.error('align_pdb_to_fasta: offset invalid')
    
    sess.run('close session')
    sess.run('set bgColor white')
    current_model_id = 0
    meta = {}
    current_model_id, meta = bib.open_model(sess, input_file, \
                                            current_model_id, meta, 1)

    if (offset != 0):
        sess.run('renumber #'+str(current_model_id)+' start '+str(offset+1))

    sess.run('save "'+output_file+'" #'+str(current_model_id))

    sess.run('close session')
    sess.run('set bgColor white')

    return
