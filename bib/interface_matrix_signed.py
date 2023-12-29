import ctl
import filesystem
import geometry
import molmodel

import glob
import os


'''
Module providing functions for handling signed interface matrices.
'''

def create(interface_res, rmsds, sess):
    '''
    Calculate non-zero interface matrix (signed interface matrix) elements from
    given residues (list interface_res).
    Distances above a cutoff value are not considered.
    '''
    interface_residues = molmodel.interface_residues_per_model(interface_res)

    matr = {}
    cutoff = 10

    for model_id0 in interface_residues:
        for model_id1 in interface_residues:
            if model_id1[1] <= model_id0[1]:
                continue

            for r0 in interface_residues[model_id0]:
                for r1 in interface_residues[model_id1]:
                    r0_coords =  sess.get_xyz(model_id0, r0)
                    r1_coords =  sess.get_xyz(model_id1, r1)
                    d = geometry.dist(r0_coords, r1_coords)
                    key = ((model_id0[1], model_id1[1]), r0, r1)

                    if d <= cutoff:
                        matr[key] = [round(d, 2), rmsds[r0], rmsds[r1]]

    matr_sort = dict(sorted(matr.items(), key=lambda item: item[0]))

    return matr_sort


def get_pairwise_interface_residues(current_model_id, chain_n, sess):
    '''
    Get interface residues between the monomers in relation to each other.
    '''
    mat_interface_residues = \
                [[None for j in range(0, chain_n)] for i in range(0, chain_n)]

    for i in range(0, chain_n):
        for j in range(0, chain_n):
            if j >= i:
                continue

            interface_residues = sess.get_interface_residues( \
                    (current_model_id, i+1), (current_model_id, 1+j))
            mat_interface_residues[i][j] = interface_residues

    return mat_interface_residues


def get_corr_coefficient(distogram0, distogram1):
    '''
    Calculate correlation coefficients between 2 given interface distograms.
    '''

    # The order of the axes of both interface distograms can be different
    # because the order of the underlying interface monomers is not defined.
    # Therefore, correlations for both possible overlay orientations are
    # calculated.


    # get correlation of overlay orientation 1
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nonzero_n = 0
    corr_sum = 0


    # get correlation of overlay in orientation 2 (reverse)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nonzero_n1 = 0
    corr_sum1 = 0


    for d0 in distogram0:
        for d1 in distogram1:

            # get correlation of overlay orientation 1
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if d0 == d1:
                rel = min(distogram0[d0][0], distogram1[d1][0])/ \
                      max(distogram0[d0][0], distogram1[d1][0])

                if rel < 0:
                    ctl.error('get_corr_coefficient: interface_correlation')

                corr_sum += rel
                nonzero_n += 1


            # get correlation of overlay in orientation 2 (reverse)
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if d0[0] == d1[1] and d0[1] == d1[0]:
                rel = min(distogram0[d0][0], distogram1[d1][0])/ \
                      max(distogram0[d0][0], distogram1[d1][0])

                if rel < 0:
                    ctl.error('get_corr_coefficient: interface_correlation')

                corr_sum1 += rel
                nonzero_n1 += 1


    # get correlation of overlay orientation 1
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if nonzero_n == 0:
        corr_coeff = 0
    else:
        corr_coeff = corr_sum/nonzero_n


    # get correlation of overlay in orientation 2 (reverse)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if nonzero_n1 == 0:
        corr_coeff1 = 0
    else:
        corr_coeff1 = corr_sum1/nonzero_n1


    if corr_coeff1 > corr_coeff:
        corr_coeff = corr_coeff1
        nonzero_n = nonzero_n1

    if nonzero_n == 0:
        return 0, 0


    max_corr_coeff = 50 # 50 not reached in typical predictions

    if corr_coeff > max_corr_coeff:
        ctl.e('corr_coeff')
        ctl.e(corr_coeff)
        ctl.error('interface_correlation: corr_fact > max_corr_coeff')

    return corr_coeff, nonzero_n


def corr_coefficients(interface_distograms):
    '''
    Calculate correlation coefficients between all combinations of given
    interface distograms.
    '''
    corr_coefficients = {}

    for i,d0 in enumerate(interface_distograms):
        for j,d1 in enumerate(interface_distograms):
            if j <= i:
                continue

            corr_coefficient, nonzero_n = get_corr_coefficient( \
                    interface_distograms[d0], interface_distograms[d1])

            corr_coefficients[(d0, d1)] = (corr_coefficient, nonzero_n)

    return corr_coefficients


def export_to_txt(matrix_elements, path_file):
    '''
    Export signed interface matrix to textfile.
    Each line contains a nonzero matrix element. Empty (zero) matrix elements
    are not exported.
    '''
    f = open(path_file, 'w')

    for key in matrix_elements:
        f.write(str(key[0][0])+";"+str(key[0][1])+";"+ \
                str(key[1])+";"+str(key[2])+";"+ \
                str(matrix_elements[key][0])+";"+ \
                str(matrix_elements[key][1])+";"+ \
                str(matrix_elements[key][2])+";"+"\r\n")

    f.close()

    return


def load(path, excl_res=[]):
    '''
    Load interface distogram.
    '''
    fn = path.split('/')[-1]
    fn_part1 = '_'.join(fn.split('_')[0:-1])
    path_dir = os.path.dirname(path)+'/'
    path_dir2 = filesystem.clean_path(path_dir+'../interfaces/')

    file = path_dir2+fn_part1+'_interface_signed.txt'
    files = sorted(glob.glob(file))

    if len(files) != 1:
        return {}

    f = filesystem.get_file(file)

    matr = {}

    for i,l in enumerate(f):
        l2 = l.split(';')
        if len(l2) >= 2:
            if int(l2[2]) not in excl_res and int(l2[3]) not in excl_res:
                key = (int(l2[2]), int(l2[3]))
                matr[key] = [round(float(l2[4]), 2), \
                             round(float(l2[5]), 2), \
                             round(float(l2[6]), 2)]

    return matr
