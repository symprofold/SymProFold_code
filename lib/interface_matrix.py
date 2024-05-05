import ctl
import filesystem
import geometry

import glob
import os


'''
Module providing functions for handling interface matrices.
'''

def create(interface_res, mate_distances, sess):
    '''
    Calculate non-zero interface matrix elements from given residues
    (list interface_res).
    Distances above a cutoff value are not considered.
    '''
    matr = {}
    cutoff = 10

    for r0 in interface_res:
        for r1 in interface_res:
            if r1[1] <= r0[1]:
                continue

            r0_coords =  sess.get_xyz(r0[0], r0[1])
            r1_coords =  sess.get_xyz(r1[0], r1[1])
            d = geometry.dist(r0_coords, r1_coords)
            key = (r0[1], r1[1])

            if d <= cutoff:
                matr[key] = [round(d, 2), \
                             mate_distances[r0[1]], mate_distances[r1[1]]]

    matr_sort = dict(sorted(matr.items(), key=lambda item: item[0]))

    return matr_sort


def get_corr_coefficient(distogram0, distogram1):
    '''
    Calculate correlation coefficients between 2 given interface distograms.
    '''
    nonzero_n = 0
    corr_sum = 0
    
    for d0 in distogram0:
        for d1 in distogram1:
            if d0 == d1:
                rel = min(distogram0[d0][0], distogram1[d1][0])/ \
                      max(distogram0[d0][0], distogram1[d1][0])

                if rel < 0:
                    ctl.error('get_corr_coefficient: interface_correlation')

                corr_sum += rel
                nonzero_n += 1

    if nonzero_n == 0:
        return 0, 0

    corr_coeff = corr_sum/(nonzero_n**(1/2))

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
    Export interface matrix to textfile.
    Each line contains a nonzero matrix element. Empty (zero) matrix elements
    are not exported.
    '''
    f = open(path_file, 'w')

    for key in matrix_elements:
        f.write(str(key[0])+";"+str(key[1])+";"+ \
                str(matrix_elements[key][0])+";"+ \
                str(round(matrix_elements[key][1], 2))+";"+ \
                str(round(matrix_elements[key][2], 2))+";"+"\r\n")

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

    file = path_dir2+fn_part1+'_interface.txt'
    files = sorted(glob.glob(file))

    if len(files) != 1:
        return {}

    f = filesystem.get_file(file)

    matr = {}

    for i,l in enumerate(f):
        l2 = l.split(';')
        if len(l2) >= 2:
            if int(l2[0]) not in excl_res and int(l2[1]) not in excl_res:
                key = (int(l2[0]), int(l2[1]))
                matr[key] = [round(float(l2[2]), 2), \
                             round(float(l2[3]), 2), \
                             round(float(l2[4]), 2)]

    return matr
