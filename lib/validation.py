import bib
import ctl
import export
import filesystem

import math
import os


def validation(model_id, model_id_restrict, axes, export_path, conf, sess):
    '''
    Validation of created models and export of validation data.
    '''
    model_id = sess.model_reg.convert_model_id(model_id)
    model_id_restrict = sess.model_reg.convert_model_id(model_id_restrict)

    model_id_str = sess.model_reg.convert_model_id_to_str(model_id)
    model_id_restrict_str = sess.model_reg. \
                                convert_model_id_to_str(model_id_restrict)

    model_res_n = sess.get_model_res_n(model_id)

    if export_path[0:len(conf.export_path)] != conf.export_path:
        ctl.e(export_path)
        ctl.e(conf.export_path)
        ctl.error('validation: path error')

    path_part2 = export_path[len(conf.export_path):]

    path_validation_file = conf.get_struct_coll_meta_path()+ \
                           conf.species+'_'+conf.species_name+'/'+ \
                           path_part2+'_'+conf.meta_file_prefix

    if os.path.exists(conf.path_ax_predictions+conf.species+'/'+ \
                                        conf.symplex_path+'setting_meta.txt'):
        filesystem.create_folder([path_validation_file])


    # determination of clashes
    # ~~~~~~~~~~~~~~~~~~~~~~~~
    clashes = sess.run('clashes #'+str(model_id_str)+' & ~@H* '+ \
                       'restrict #'+str(model_id_restrict_str)+' & ~@H* '+ \
                       'intraModel true interModel true interSubmodel true')
            # exclude clashes involving hydrogen atoms (performance)
            #
            # notice: a ChimeraX model including the clashes is created

    clash_cutoff = model_res_n*2*2

    clashes_residues_n_allatoms, clashes_residues_n_woh, \
            clashes_atoms_n_allatoms, clashes_atoms_n_woh, \
            clashes_backboneatoms_n_woh, clash_cutoff_exceeded = \
                                        analyze_clashes(clashes, clash_cutoff)
            # clashes involving hydrogen atoms are excluded anyhow in this
            # function (performance)

    if clashes_residues_n_allatoms != clashes_residues_n_woh:
        ctl.e('clashes_residues_n_allatoms')
        ctl.e(clashes_residues_n_allatoms)

        ctl.e('clashes_residues_n_woh')
        ctl.e(clashes_residues_n_woh)

        ctl.e('clashes_atoms_n_allatoms')
        ctl.e(clashes_atoms_n_allatoms)

        ctl.e('clashes_atoms_n_woh')
        ctl.e(clashes_atoms_n_woh)

        ctl.error('validation: '+ \
                  'clash number not consistent with excluded hydrogens')


    if clashes_atoms_n_allatoms != clashes_atoms_n_woh and \
       clash_cutoff_exceeded == False:
        ctl.e('clashes_residues_n_allatoms')
        ctl.e(clashes_residues_n_allatoms)

        ctl.e('clashes_residues_n_woh')
        ctl.e(clashes_residues_n_woh)

        ctl.e('clashes_atoms_n_allatoms')
        ctl.e(clashes_atoms_n_allatoms)

        ctl.e('clashes_atoms_n_woh')
        ctl.e(clashes_atoms_n_woh)

        ctl.error('validation: '+ \
                  'clash number not consistent with excluded hydrogens')


    if clashes_atoms_n_allatoms%2 == 1 and clash_cutoff_exceeded == False:
        ctl.e(clashes_atoms_n_allatoms)
        ctl.error('validation: odd number of clashing atoms')

    if clashes_atoms_n_woh%2 == 1 and clash_cutoff_exceeded == False:
        ctl.e(clashes_atoms_n_woh)
        ctl.error('validation: odd number of clashing atoms')

    if clashes_backboneatoms_n_woh%2 == 1 and clash_cutoff_exceeded == False:
        ctl.e(clashes_backboneatoms_n_woh)
        ctl.error('validation: odd number of clashing atoms')

    clashes_between_atoms_n_allatoms = int(clashes_atoms_n_allatoms/2)
    clashes_between_atoms_n_woh = int(clashes_atoms_n_woh/2)
    clashes_between_backboneatoms_n_woh = int(clashes_backboneatoms_n_woh/2)


    clashes = clashes_between_atoms_n_woh
    clashes_per_100residues = 100*clashes/model_res_n
    clashes_per_residue = clashes/model_res_n

    backbone_clashes_per_residue = \
                            clashes_between_backboneatoms_n_woh/model_res_n


    # determination of layer bending
    # ------------------------------
    ax1_tilt_sum_delta_sq = 0
            # ax1_tilt: axial tilt between rotational symmetry axis and
            # global z axis
    ax1_tilt_sum_n = 0
    ax1_tilt_delta_str = ''
    ax1_tilt_rmsd = 0
    ax1_tilt_score = 0

    sum_delta_sq = 0
    sum_n = 0
    delta_str = ''
    rmsd = 0

    if len(axes) > 1:

        # calculate rmsd of axial tilt of ax1
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for ax1_mainmodel_id in range(2, 2+axes[0].fold):
            rotsymm_axis = axes[1].get_representation(ax1_mainmodel_id). \
                                                            rotsymm_ax.axis

            ax1_tilt = axes[1].get_representation(ax1_mainmodel_id). \
                                                       rotsymm_ax.get_tilt()
            ax1_tilt_sum_delta_sq += ax1_tilt**2
            ax1_tilt_sum_n += 1

            ax1_tilt_deg = math.degrees(ax1_tilt)
            ax1_tilt_delta_str += str(round(ax1_tilt_deg, 1))+'° '

        ax1_tilt_rmsd = math.sqrt( ax1_tilt_sum_delta_sq/ax1_tilt_sum_n )
        ax1_tilt_score = ax1_tilt_rmsd/(math.pi/2) # score (normalized)

        if max(conf.flatten_modes) >= 0:

            # calculate rmsd of deflection in z-direction from the xy plane
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            for ax0_mainmodel_id in range(2+axes[0].fold, \
                                          2+axes[0].fold+axes[0].fold):
                z = axes[0].get_representation(ax0_mainmodel_id).trans_vect[2]
                sum_delta_sq += z**2
                sum_n += 1
                delta_str += str(round(z, 1))+'Å '

            rmsd = math.sqrt( sum_delta_sq/sum_n )


    quality_score = clashes_per_residue+ax1_tilt_score
    overall = clashes_per_100residues+rmsd

    if os.path.exists(conf.path_ax_predictions+conf.species+'/'+ \
                                        conf.symplex_path+'setting_meta.txt'):
        f = open(path_validation_file+ \
                'qual'+str(round(quality_score, 3))+ \
                '_cl'+str(round(clashes_per_residue, 3))+ \
                '_be'+str(round(ax1_tilt_score, 3))+ \
                 '.txt', 'w', encoding='utf-8')
        f.write('clashes: '+str(clashes)+"\r\n")  
        f.write('intermolecular clashes per 100 residues: '+ \
                str(round(clashes_per_100residues, 2))+"\r\n")
        f.write('number of residues: '+str(model_res_n)+"\r\n")  

        if len(axes) > 1:
            f.write("\r\n")  
            f.write('score_clash (intermolecular clashes per residue): '+ \
                            str(round(clashes_per_residue, 4))+"\r\n")
            f.write('score_bend (layer bending, axis tilt): '+ \
                            str(round(ax1_tilt_score, 3))+"\r\n")
            f.write('individual values layer bending (axis tilt): '+ \
                            ax1_tilt_delta_str+"\r\n")
            f.write('quality score (score_clash+score_bend): '+ \
                            str(round(quality_score, 3))+"\r\n")

            if max(conf.flatten_modes) >= 0:
                f.write("\r\n")  
                f.write('score_bendz (layer bending, z deflection): '+ \
                                str(round(rmsd, 1))+"\r\n")
                f.write('individual values layer bending (z deflection): '+ \
                                delta_str+"\r\n")
                f.write('quality score (score_clash100+score_bendz): '+ \
                                str(round(overall, 1))+"\r\n")

        f.close()

    return [clashes_per_100residues, \
            [quality_score, clashes_per_residue, ax1_tilt_score, \
            backbone_clashes_per_residue]]


def analyze_clashes(clashes, clash_cutoff=None, without_h=True):
    '''
    Determine number of clashes for different clash types.
    '''
    clashes_residues_n_allatoms = 0
            # number of clashing residues,
            # considering all types of atoms (including hydrogens)

    clashes_residues_n_woh = 0
            # number of clashing residues,
            # considering all types of atoms except hydrogens

    clashes_atoms_n_allatoms = 0
            # number of clashing atoms,
            # considering all types of atoms (including hydrogens)

    clashes_atoms_n_woh = 0
            # number of clashing atoms,
            # considering all types of atoms except hydrogens

    clashes_backboneatoms_n_woh = 0
            # number of clashing backbone atoms,
            # considering all types of atoms except hydrogens

    clash_cutoff_exceeded = False
            # if clash cutoff exceeded, continue with backbone clashes only

    for clash_res in clashes:
        perres_clashes_atoms_n_allatoms = 0
                # number of clashing atoms per residue (clash_res),
                # considering all types of atoms (including hydrogens)

        perres_clashes_atoms_n_woh = 0
                # number of clashing atoms per residue (clash_res),
                # considering all types of atoms except hydrogens

        perres_clashes_backboneatoms_n_woh = 0
                # number of clashing backbone atoms per residue (clash_res),
                # considering all types of atoms except hydrogens

        # check if clash cutoff is exceeded (performance)
        if clash_cutoff != None and clashes_atoms_n_woh > clash_cutoff:
            clash_cutoff_exceeded = True

        clash_res_ = str(clash_res).strip().split(' ')
        clash_res_name = clash_res_[-3].strip().lower()
                # name of clashing residue, e.g. 'ala'
        clash_atom = clash_res_[-1].strip()
                # clashing atom, e.g. 'C'

        if without_h == True and clash_atom[0] == 'H':
            ctl.error('analyze_clashes: hydrogen (H atom) in list of clashes')

        if clash_cutoff_exceeded == False or \
            clash_cutoff_exceeded == True and clash_atom in ('N', 'CA', 'C'):
                    # in the case that clash cutoff is exceeded,
                    # continue with backbone clashes only (performance)

            for clash_partner in clashes[clash_res]:
                clash_partner_ = str(clash_partner).strip().split(' ')
                clash_partner_atom = clash_partner_[-1].strip()

                # skip clashes between cysteines to avoid counting clashes
                # caused by unrecognized disulfide bonds
                if clash_res_name == 'cys':
                    if clash_partner_[-3].strip().lower() == 'cys':
                        continue

                if without_h == True and clash_partner_atom[0] == 'H':
                    ctl.error('analyze_clashes: '+ \
                              'hydrogen (H atom) in list of clashes')

                perres_clashes_atoms_n_allatoms += 1

                if clash_atom[0] != 'H' and clash_partner_atom[0] != 'H':
                    perres_clashes_atoms_n_woh += 1

                if clash_atom in ('N', 'CA', 'C') and \
                                clash_partner_[-1].strip() in ('N', 'CA', 'C'):
                    perres_clashes_backboneatoms_n_woh += 1


            clashes_atoms_n_allatoms += perres_clashes_atoms_n_allatoms
            clashes_atoms_n_woh += perres_clashes_atoms_n_woh
            clashes_backboneatoms_n_woh += perres_clashes_backboneatoms_n_woh

            clashes_residues_n_allatoms += \
                                        min(1, perres_clashes_atoms_n_allatoms)
            clashes_residues_n_woh += min(1, perres_clashes_atoms_n_woh)


    if clash_cutoff != None and clashes_atoms_n_woh > clash_cutoff:
        clashes_atoms_n_woh = clash_cutoff

    return clashes_residues_n_allatoms, clashes_residues_n_woh, \
            clashes_atoms_n_allatoms, clashes_atoms_n_woh, \
            clashes_backboneatoms_n_woh, clash_cutoff_exceeded
