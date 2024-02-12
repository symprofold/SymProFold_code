import ctl
import bib
import bibpdb
import geometry
import interface_matrix
import interface_matrix_signed
import model_preparation.betasheet
import model_preparation.clashes
import model_preparation.chain_numbering
import model_preparation.pdb
import model_preparation.rotsymm_ax

from model_preparation.rotsymm_ax import RotSymmAxes

import os
import params
import shutil


'''
Module providing functions for model preparation of a single coord file.
'''

def prepare_file(f0, seq, path_export, sess, verbous=False, iteration=0):
    '''
    Prepare a single coord file.
    '''
    f = f0
    filename = f.split('/')[-1]
    ctl.d(filename)


    # get type of microorganism (bacteria, virus)
    microorganism_type = 0 # 0: bacteria and archaea, 1: virus

    if 'virus' in f:
        microorganism_type = 1


    model_preparation.pdb.align_to_fasta(f, path_export+filename, seq, sess)
    f = path_export+filename

    sess.run('close session')
    sess.run('set bgColor white')

    current_model_id = 0
    meta = {} 
    bib.place_plane(sess)

    symmaxis = True # unique axis of rotational symmetry found


    # preprocessing of pdb input file
    model_preparation.chain_numbering.realign_chain_ids(f0, f, f)
        
    current_model_id, meta = bib.open_model(sess, f, current_model_id, \
                                            meta, 1, \
                                            termini_with_signalsequence=False)

    # continue if whole model disordered (large rmsds)
    if meta[current_model_id][1][0] == -1 or \
       meta[current_model_id][1][1] == 10000:

        return


    # determination of number of chains in 2 ways and check for agreement
    resids = sess.resids(current_model_id)
    res_n = sess.get_res_n((current_model_id,))
    mol_count = res_n/(resids[-1]-resids[0]+1)

    if mol_count%1 > 0.001:
        ctl.e(mol_count)
        ctl.error('unclear number of chains')
    else:
        mol_count = int(mol_count)

    mol_count2 = bib.get_multimer_n(f)

    if mol_count != mol_count2:
        ctl.e(mol_count)
        ctl.e(mol_count2)
        ctl.error('unclear number of chains')

    # continue if not at least 2 chains
    if mol_count <= 1:

        return


    # get clashes per 100aa in a rolling range of 200aa along the whole
    # chain.
    roll_clash_postfix, zero = model_preparation.clashes. \
                                    clashes_per_100_rolling(
                                        current_model_id, 200, mol_count, sess)

    # get clashes per 100aa
    clash_postfix, zero = model_preparation.clashes. \
                                    clashes_per_100(current_model_id, sess)


    # check if parts of the protein chain run through the sphere of a
    # wrong monomer
    #
    # This check evaluates the fraction of residues involved in medium-length
    # intermolecular beta sheet strands.
    sheetintermol_infix, zero, zero = \
            model_preparation.betasheet.intermolecular_betasheet_fraction( \
                                current_model_id, path_export, sess, verbous)


    # alignment of whole model to xy plane
    # ------------------------------------
    bib.model_to_plane(meta, current_model_id, sess)
    sess.run('save "'+f+'" #'+str(current_model_id))

    coord = bibpdb.open_pdb(f, False)
    cen = geometry.get_center(coord, 'all', meta[current_model_id][1])

    sess.run('move x '+str((-1)*cen[0])+' models #'+str(current_model_id))
    sess.run('move y '+str((-1)*cen[1])+' models #'+str(current_model_id))
    sess.run('save "'+f+'" #'+str(current_model_id))

    if cen[2] > 0:
        sess.run('turn x 180 models #'+str(current_model_id))
        sess.run('save "'+f+'" #'+str(current_model_id)+' ')

        current_model_id, meta = bib.open_model( \
                sess, f, current_model_id, meta, 1, \
                termini_with_signalsequence=False)
         
        coord = bibpdb.open_pdb(f, False)
        cen = geometry.get_center(coord, 'all', meta[current_model_id][1])
        sess.run('move x '+str((-1)*cen[0])+' models #'+str(current_model_id))
        sess.run('move y '+str((-1)*cen[1])+' models #'+str(current_model_id))
        sess.run('save "'+f+'" #'+str(current_model_id)+' ')  

    sess.run('marker #110 position '+ \
             str(cen[0])+','+str(cen[1])+','+str(cen[2])+ \
             ' color yellow radius 1')


    # counterclockwise renaming of chain ids
    # --------------------------------------
    if mol_count >= 3:
        coord = bibpdb.open_pdb(f, False)     
        complexcenter = geometry.get_center(coord, 'all', \
                                            meta[current_model_id][1])
        order = model_preparation.chain_numbering. \
                    determine_subchain_order(coord, meta[current_model_id], \
                                             mol_count)

        # reverse order of chain ids if orientation is not counterclockwise
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if order == -1:
            # preparation: flip order when clockwise
            multimer_n = mol_count

            order_flipped = []
            for i in range(multimer_n-1, -1, -1):
                order_flipped.append(i)

            bibpdb.reorder_chainids(path_export+filename, order_flipped)


        coord = bibpdb.open_pdb(f, False)
        complexcenter = geometry.get_center(coord, 'all', \
                                            meta[current_model_id][1])

        order = model_preparation.chain_numbering. \
                    determine_subchain_order(coord, meta[current_model_id], \
                                             mol_count)

        if order == -1:
            ctl.d('order')
            ctl.d(order)
            ctl.d('still: order == -1, no symm axis detectable')
            symmaxis = False

        ctl.d('symmaxis')
        ctl.d(symmaxis)     

        if symmaxis:
          ctl.d(order)
          bibpdb.reorder_chainids(f, order)
     
    ctl.d(mol_count)
    ctl.d(symmaxis)


    try:
        os.mkdir(path_export+'clashes/')
    except IOError:
        pass

    f_cl = open(path_export+'clashes/'+ \
                filename[:-4]+clash_postfix+'.pdb', 'w')
    f_cl.write('')
    f_cl.close()

    f_rollcl = open(path_export+'clashes/'+ \
                filename[:-4]+roll_clash_postfix+'.txt', 'w')
    f_rollcl.write('')
    f_rollcl.close()


    # determination of the symmetry of SymPlex candidate
    # ==================================================

    # determination  of rotation angles between all monomers in the
    # SymPlex candidate
    # ------------------------------------------------------------------
    # Rotation angles are determined by superposition of the interfaces.

    current_model_id, meta = \
            bib.open_model(sess, f, current_model_id, meta, \
                           termini_with_signalsequence=False)
    current_model_id, meta = \
            bib.open_model(sess, f, current_model_id, meta, \
                           termini_with_signalsequence=False)

    # determination of interface residues
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    interface_residues = sess.get_interface_residues( \
                            (current_model_id, 1), (current_model_id, 2))

    interface_residues = sorted([r[1] for r in interface_residues])

    if len(interface_residues) < 1:
        interface_res_range = meta[current_model_id][1]
    else:
        interface_res_range = [interface_residues[0], \
                               interface_residues[-1]]

        if interface_res_range[1]-interface_res_range[0] < 5:
            interface_res_range = meta[current_model_id][1]


    rotang, rot_axis, rot_axis_point, rmsd, mate_distances = \
            model_preparation.geometry.get_nn_rotation(current_model_id-1, \
                        current_model_id, interface_res_range, mol_count, \
                        sess)

    rot_axes_point_average, uniqueness = \
            model_preparation.rotsymm_ax.check_uniqueness( \
                                    rot_axis_point, rot_axis, f, \
                                    True if mol_count >= 3 else False)

    if mol_count >= 3:
        if uniqueness == False:
            ctl.d('no unique rotation symmetry axis')
            f_info = open(f+'_no_unique_rot_axis.txt', 'w')
            f_info.write('')
            f_info.close()
            symmaxis = False

    ctl.d(rotang)
    ctl.d(rot_axis)
    ctl.d(rot_axis_point)


    # determination of rotation angles
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if not symmaxis:
        # no unique axis of rotational symmetry found

        ctl.d('no symm axis')
        ctl.d(symmaxis)  

        rotang_mat, rot_axis_mat, rot_axis_point_mat = \
                model_preparation.geometry.get_pairwise_rotation( \
                        current_model_id-1, \
                        current_model_id, interface_res_range, mol_count, 
                        sess)

        if verbous:
            ctl.p(rotang)
            ctl.p(rot_axis)
            ctl.p(rot_axis_point)
            ctl.p('rotang_mat')
            ctl.p(rotang_mat)
            ctl.p('rot_axis_mat')
            ctl.p(rot_axis_mat)
            ctl.p('rot_axis_point_mat')
            ctl.p(rot_axis_point_mat)

        try:
            os.mkdir(path_export+'no_symm_found/')
        except IOError:
            pass
        try:
            os.mkdir(path_export+'no_symm_found/clashes/')
        except IOError:
            pass       
        try:
            os.mkdir(path_export+'no_symm_found/sheetintermol/')
        except IOError:
            pass

        shutil.copyfile(f, path_export+'no_symm_found/'+filename)
        shutil.copyfile(f, path_export+'no_symm_found/clashes/'+ \
                        filename[:-4]+clash_postfix+'.pdb')

        f_rollcl = open(path_export+'no_symm_found/clashes/'+ \
                        filename[:-4]+roll_clash_postfix+'.txt', 'w')
        f_rollcl.write('')
        f_rollcl.close()

        f_cl = open(path_export+'no_symm_found/sheetintermol/'+ \
                    filename[:-4]+sheetintermol_infix+'.txt', \
                    'w')
        f_cl.write('')
        f_cl.close()


        if mol_count > 2 and microorganism_type == 1 and iteration == 0:

            # determine interface residues to avoid processing of oligomers
            # w/o contact
            interface_res_mat = interface_matrix_signed. \
                    get_pairwise_interface_residues( \
                            current_model_id, mol_count, sess)

            # determination of rotational symmetry axes
            rot_axes = RotSymmAxes(rot_axis_point_mat, rotang_mat, \
                                   interface_res_mat)
            rot_axes.path_export = path_export
            rot_axes.path = f
            rot_axes.filename = filename
            rot_axes.coinciding_axes()

            # export all SymPlexes in coord file f to separate coord files
            files_written = rot_axes.write_symplexes(current_model_id, sess)

            for file in files_written:
                prepare_file(file, seq, path_export, sess, False, iteration+1)

    else:
        # unique axis of rotational symmetry found

        # export of rmsd between two monomers in the SymPlex candidate
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        try:
            os.mkdir(path_export+'interfaces/')
        except IOError:
            pass  

        f_info = open(path_export+'interfaces/'+filename[:-4]+'_rmsd.txt', 'w')
        f_info.write(str(rmsd))
        f_info.close()

        current_model_id, meta = bib.open_model( \
                                        sess, f, \
                                        current_model_id, meta, \
                                        termini_with_signalsequence=False)
        sess.run('move x '+str((-1)*rot_axes_point_average[0])+ \
                 ' models  #'+str(current_model_id))
        sess.run('move y '+str((-1)*rot_axes_point_average[1])+ \
                 ' models  #'+str(current_model_id))  
        sess.run('save "'+f+'" #'+str(current_model_id)+' ')  
        bibpdb.clean_pdb(f)

        if mol_count >= 2:
            multiplicities = params.get_multiplicities(path_export+'../', \
                                                microorganism_type)

            symm_found = []

            for multiplicity in multiplicities:
                symm_ang = 360/multiplicity
                symm_ang_infix = ('000'+str(round(symm_ang)))[-3:]

                rotang_cleaned = model_preparation.geometry. \
                                clean_rotang(rotang, multiplicity, mol_count)

                d0 = 15
                d1 = 15
                if 5 in multiplicities:
                    if multiplicity == 5:
                        d0 = 6
                        d1 = 9
                    elif multiplicity == 6:
                        d1 = 6
                    elif multiplicity == 4:
                        d0 = 9

                if 7 in multiplicities:
                    if multiplicity == 7:
                        d1 = 4
                    elif multiplicity == 6:
                        d0 = 4

                # rule out that e.g. 4mer results in 2fold axis
                if mol_count <= multiplicity:
                    if abs(max(rotang_cleaned)-symm_ang) < d0 and \
                       abs(min(rotang_cleaned)-symm_ang) < d1:
                        try:
                            os.mkdir(path_export+'symm_'+symm_ang_infix+'/')
                        except IOError:
                            pass
                        try:
                            os.mkdir(path_export+'symm_'+symm_ang_infix+ \
                                     '/clashes/')
                        except IOError:
                            pass
                        try:
                            os.mkdir(path_export+'symm_'+symm_ang_infix+ \
                                     '/sheetintermol/')
                        except IOError:
                            pass

                        shutil.copyfile(f, path_export+'symm_'+ \
                                        symm_ang_infix+'/'+filename)
                        f_cl = open(path_export+'symm_'+ \
                                    symm_ang_infix+'/clashes/'+ \
                                    filename[:-4]+clash_postfix+'.pdb', 'w')
                        f_cl.write('')
                        f_cl.close()
                        
                        rotang_postfix = '_rotang'+ \
                                     str(round(min(rotang_cleaned))). \
                                         replace('.', '-')+'_'+ \
                                     str(round(max(rotang_cleaned))). \
                                         replace('.', '-')
                        f_cl = open(path_export+'symm_'+ \
                                    symm_ang_infix+'/'+ \
                                    filename[:-4]+rotang_postfix+'.txt', 'w')
                        f_cl.write('')
                        f_cl.close()

                        f_cl = open(path_export+'symm_'+ \
                                    symm_ang_infix+'/clashes/'+ \
                                    filename[:-4]+roll_clash_postfix+'.txt', \
                                    'w')
                        f_cl.write('')
                        f_cl.close()

                        f_cl = open(path_export+'symm_'+ \
                                    symm_ang_infix+'/sheetintermol/'+ \
                                    filename[:-4]+sheetintermol_infix+'.txt', \
                                    'w')
                        f_cl.write('')
                        f_cl.close()

                        try:
                            os.mkdir(path_export+'symm_'+symm_ang_infix+ \
                                     '/interfaces/')
                        except IOError:
                            pass  


                        # calculate and export non-zero interface
                        # matrix elements for the first two monomers in the
                        # SymPlex candidate
                        # -------------------------------------------------
                        current_model_id, meta = bib.open_model( \
                                sess, f, current_model_id, meta, \
                                termini_with_signalsequence=False)

                        interface_residues = sess.get_interface_residues( \
                                (current_model_id, 1), (current_model_id, 2))

                        interface_mat = interface_matrix.create( \
                                interface_residues, mate_distances, sess)

                        fn = path_export+ \
                                'symm_'+symm_ang_infix+'/interfaces/'+ \
                                filename[:-4]+'_interface.txt'
                        interface_matrix.export_to_txt(interface_mat, fn)


                        # calculate and export non-zero interface
                        # matrix (signed interface matrix) elements for
                        # monomer pairs in the SymPlex candidate
                        # ---------------------------------------------
                        rmsds = bibpdb.get_rmsds(f)
                        interface_mat_signed = {}

                        for monomer0 in range(1, mol_count):
                            interface_residues = sess.get_interface_residues( \
                                (current_model_id, monomer0), \
                                (current_model_id, monomer0+1))

                            interface_mat_signed_ = interface_matrix_signed. \
                                        create(interface_residues, rmsds, sess)

                            interface_mat_signed.update(interface_mat_signed_)

                        if mol_count >= 3:
                            interface_residues = \
                                sess.get_interface_residues( \
                                (current_model_id, 1), \
                                (current_model_id, mol_count))

                            interface_mat_signed_ = interface_matrix_signed. \
                                        create(interface_residues, rmsds, sess)

                            interface_mat_signed.update(interface_mat_signed_)


                        fn = path_export+ \
                                'symm_'+symm_ang_infix+'/interfaces/'+ \
                                filename[:-4]+'_interface_signed.txt'

                        interface_matrix_signed.export_to_txt( \
                                interface_mat_signed, fn)


                        symm_found.append(multiplicity)
                        break


            if symm_found == []:
                try:
                    os.mkdir(path_export+'symm_unclassified/')
                except IOError:
                    pass
                try:
                    os.mkdir(path_export+'symm_unclassified/clashes/')
                except IOError:
                    pass            
                try:
                    os.mkdir(path_export+'symm_unclassified/sheetintermol/')
                except IOError:
                    pass

                shutil.copyfile(f, path_export+'symm_unclassified/'+ \
                        filename[:-4]+'_'+ \
                        str(round(min(rotang)))+'-'+str(round(max(rotang)))+ \
                        '.pdb')
                shutil.copyfile(f, path_export+'symm_unclassified/clashes/'+ \
                        filename[:-4]+'_'+ \
                        str(round(min(rotang)))+'-'+str(round(max(rotang)))+ \
                        clash_postfix+'.pdb')

                f_rollcl = open(path_export+'symm_unclassified/clashes/'+ \
                        filename[:-4]+'_'+ \
                        str(round(min(rotang)))+'-'+str(round(max(rotang)))+ \
                        roll_clash_postfix+'.txt', 'w')
                f_rollcl.write('')
                f_rollcl.close()

                f_cl = open(path_export+'symm_unclassified/sheetintermol/'+ \
                            filename[:-4]+sheetintermol_infix+'.txt', \
                            'w')
                f_cl.write('')
                f_cl.close()
