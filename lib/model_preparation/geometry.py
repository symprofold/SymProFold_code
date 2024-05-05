import chimerax_api
import ctl
import geometry
import math
import molmodel
import statistics


'''
Module providing gemetric functions for model preparation.
'''

def get_distances(points, a):
    ''' Get distance of point a to each point in list "points". '''

    dist = []
        
    for c in points:
        l = geometry.dist(c, a)
        if l == 0:
            l = 1000000

        dist.append(l)

    return dist


def get_pairwise_rotation(model_id, model_id2, res_range, mult, sess):
    '''
    Get the rotation of the monomers in relation to each other.

    list of rotation axes, axis points, rotation angles.
    '''
    res_start = res_range[0]
    res_end = res_range[1]


    # calculate matrices with rotation axes, axis points, rotation angles
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rotang_mat = [[None for j in range(0, mult)] for i in range(0, mult)]
    rot_axis_mat = [[None for j in range(0, mult)] for i in range(0, mult)]
    rot_axis_point_mat = \
                    [[None for j in range(0, mult)] for i in range(0, mult)]

    # model_id: model to use for movements caused by match command
    # model_id2: reference model
    for i in range(1, mult+1):
        for j in range(1, mult+1):
            if j >= i:
                continue

            # reset monomer of model_id to the original position
            sess.run('match #'+str(model_id)+'.'+str(i)+':'+ \
                     str(res_start)+'-'+str(res_end)+ \
                     ' to #'+str(model_id2)+'.'+str(i))

            # measure rotation axes, axis points, rotation angles
            sess.run('match #'+str(model_id)+'.'+str(i)+':'+ \
                     str(res_start)+'-'+str(res_end)+ \
                     ' to #'+str(model_id2)+'.'+str(j))
            re = sess.run('measure rotation #'+ \
                          str(model_id)+'.'+str(i)+':'+ \
                          str(res_start)+'-'+str(res_end)+ \
                          ' toModel #'+ \
                          str(model_id2)+'.'+str(j)+' showAxis false')

            rotangle = re.rotation_angle()/math.pi*180
            rotang_mat[i-1][j-1] = rotangle
            rot_axis_mat[i-1][j-1] = chimerax_api.get_rot_axis( \
                                                re.description())
            rot_axis_point_mat[i-1][j-1] = chimerax_api.get_rot_axis_point( \
                                                re.description())

    return rotang_mat, rot_axis_mat, rot_axis_point_mat


def get_nn_rotation(model_id, model_id2, res_range, mult, sess):
    '''
    Get the rotation of next neighbor monomers in relation to each other.
    '''
    res_start = res_range[0]
    res_end = res_range[1]

    rotang = []
    rot_axis = []
    rot_axis_point = []

    for i in range(mult, 1, -1):
        sess.run('match #'+str(model_id)+'.'+str(i)+':'+ \
                 str(res_start)+'-'+str(res_end)+ \
                 ' to #'+str(model_id)+'.'+str(i-1))
        re = sess.run('measure rotation #'+ \
                      str(model_id)+'.'+str(i)+':'+ \
                      str(res_start)+'-'+str(res_end)+ \
                      ' toModel #'+ \
                      str(model_id)+'.'+str(i-1)+' showAxis false')
        rotangle = re.rotation_angle()/math.pi*180

        rotang.append(rotangle)
        rot_axis.append(chimerax_api.get_rot_axis(re.description()))
        rot_axis_point.append(chimerax_api. \
                              get_rot_axis_point(re.description()))

    sess.run('match #'+str(model_id)+'.'+str(1)+':'+ \
             str(res_start)+'-'+str(res_end)+ \
             ' to #'+str(model_id2)+'.'+str(mult))
    re = sess.run('measure rotation #'+ \
             str(model_id)+'.'+str(1)+':'+ \
             str(res_start)+'-'+str(res_end)+' toModel #'+ \
             str(model_id2)+'.'+str(mult)+' showAxis false')

    rotangle = re.rotation_angle()/math.pi*180
    rotang.append(rotangle)
    rot_axis.append(chimerax_api.get_rot_axis(re.description()))
    rot_axis_point.append(chimerax_api.get_rot_axis_point(re.description()))

    rmsd, mate_distances = molmodel.get_sibling_rmsd( \
                        (model_id, 1), (model_id2, mult), [res_range], sess)

    # rotang: list is in inverse order!
    return rotang, rot_axis, rot_axis_point, rmsd, mate_distances


def clean_rotang(rotang, multiplicity, mult):
    '''
    Clean rotation angles.

    Remove rotation angle(s) with the maximal difference to median.
    '''
    rotang_cleaned = []

    if multiplicity == mult:
        return rotang

    median = statistics.median(rotang)
    diffs_to_median = [abs(r-median) for r in rotang]
    max_diff = diffs_to_median.index(max(diffs_to_median))

    for i,r in enumerate(rotang):
        if i != max_diff:
            rotang_cleaned.append(r)

    return rotang_cleaned
