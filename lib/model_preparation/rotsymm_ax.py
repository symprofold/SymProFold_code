import bib
import bibpdb
import ctl
import geometry
import metadata

import glob


'''
Module providing a class and functions to handle and analyze
rotational symmetry axes.
'''

class RotSymmAxes():
    '''
    This class describes the set of rotational symmetry axes in a
    oligomer model.
    '''

        
    def __init__(self, mat_rot_axis_point, mat_rotang, interface_res_mat):
        '''
        Initialization of the RotSymmAxes class.
        '''
        self.mat_rot_axis_point = mat_rot_axis_point
        self.mat_rotang = mat_rotang
        self.interface_res_mat = interface_res_mat
        self.filename = None
        self.path = None
        self.path_export = None

        self.axes = [] # [[axis points], [chain ids]]
                # chain ids starting with 0

        self.axis_maxdist = 5
                # cutoff value: max distance of axes to be considered as
                # coinciding

        return


    def get_infix_file_num_next(self):
        '''
        Get infix for file numbering.
        '''
        files = sorted(glob.glob(self.path_export+'*.pdb'))
        rank_max = 0

        for f in files:
            ctl.p(f)
            rank = metadata.get_rank(f)
            rank_max = max(rank_max, rank)

        return rank_max


    def write_symplexes(self, current_model_id, sess):
        '''
        Write (export) all SymPlexes to separate coord files.
        '''
        files_written = []

        meta = {} 
        current_model_id, meta = bib.open_model(sess, self.path, \
                                                current_model_id, meta)
        number = 0

        for ax in self.axes:
            model_ids = []

            for chain_id in ax[1]:
                model_ids.append((current_model_id, chain_id+1))

            number = max(self.get_infix_file_num_next()+1, 20)
            filename_parts = self.filename.split('rank')
            filename_ = filename_parts[0]+'rank'+('00'+str(number))[-2:]+ \
                        filename_parts[1][2:]

            sess.save_model_ids(model_ids, self.path_export+filename_, 'pdb')
            bibpdb.clean_pdb(self.path_export+filename_)
            files_written.append(self.path_export+filename_)

        return files_written


    def coinciding_axes(self):
        '''
        Identify coinciding rotation axes.
        '''
        for i in range(0, len(self.mat_rot_axis_point)):
            for j in range(0, len(self.mat_rot_axis_point)):
                point = self.mat_rot_axis_point[i][j]
                rotang = self.mat_rotang[i][j]

                if point != None:
                    if self.interface_res_mat[max(i, j)][min(i, j)] == None:
                        ctl.e(i)
                        ctl.e(j)
                        ctl.error('RotSymmAxes: coinciding_axes: '+ \
                            'interface_res_mat[max(i, j)][min(i, j)] == None')
                    elif self.interface_res_mat[max(i, j)][min(i, j)] != []:
                        self.add_axis(point, [i, j], rotang)

        return


    def add_axis(self, point, chain_ids, rotang):
        '''
        Add axis to list of already assigned axes.

        chain ids starting with 0
        '''
        assigned = False

        for ax in self.axes:
            for ax_point in ax[0]:
                d = geometry.dist(ax_point, point)

                if d <= self.axis_maxdist:
                    if assigned:
                        ctl.error('RotSymmAxes: add_axis: '+ \
                                  'double asignment of axis possible')

                    ax[0].append(point)
                    ax[1].append(chain_ids[0])
                    ax[1].append(chain_ids[1])
                    ax[1] = sorted(ax[1])
                    ax[2].append(rotang)
                    assigned = True
                    break

        if not assigned:
            self.axes.append([[point], sorted([chain_ids[0], chain_ids[1]]), \
                              [rotang]])

        return


def check_uniqueness(rot_axis_points, rot_axes, file, verbous=False):
    '''
    Check if all rotation axes are oriented in z direction and have similar
    rotation axis point.
    '''

    # cutoff value: max deviation of max_rot_axis_point in xy direction from
    # (0, 0, 0)
    rot_axis_point_max_dev_xy = 5
    
    rot_axes_point_average = [0, 0, 0]
    rot_axes_point_average_n = 0
    rot_axes_points_z = []

    uniqueness = True

    for rot_axis in rot_axes:
        rot_axes_points_z.append(rot_axes_point_average[2])

        # check all rotation axes are oriented in z direction
        if abs(rot_axis[2]) < 0.95:
            ctl.e(rot_axis)
            f = open(file+'_ERROR_rot_axis_not_perpendicular.txt', 'w')
            f.write('') 
            f.close()  
            ctl.d('check_rot_axis_uniqueness: rot_axis[2]) < 0.95')
            uniqueness = False

    for i,rot_axis_point in enumerate(rot_axis_points):
        if abs(rot_axis_point[0]) > 0.05 or abs(rot_axis_point[1]) > 0.05:
            ctl.d('check_rot_axis: deviation >0.05')

            rot_axes_point_average[0] += rot_axis_point[0]
            rot_axes_point_average[1] += rot_axis_point[1]
            rot_axes_point_average[2] += rot_axis_point[2]
            rot_axes_point_average_n += 1

    for i,rot_axis_point in enumerate(rot_axis_points):

        if i == 0:
            continue

        rot_ax_deviation = geometry.dist( \
                            (rot_axis_point[0], rot_axis_point[1], 0), \
                            (rot_axis_points[0][0], rot_axis_points[0][1], 0))

        if rot_ax_deviation > rot_axis_point_max_dev_xy:
            if verbous:
                f_info = open(file+'_rotax_dev_'+ \
                              str(round(rot_ax_deviation, 3))+ \
                              '_exceeds_cutoff.txt', 'w')
                f_info.write('') 
                f_info.close()  
            ctl.d('rotaxes deviation exceeds cutoff' )
            ctl.d(rot_ax_deviation)
            uniqueness = False

    if rot_axes_point_average_n != 0:
        rot_axes_point_average[0] = \
                rot_axes_point_average[0]/rot_axes_point_average_n
        rot_axes_point_average[1] = \
                rot_axes_point_average[1]/rot_axes_point_average_n
        rot_axes_point_average[2] = \
                rot_axes_point_average[2]/rot_axes_point_average_n

    return rot_axes_point_average, uniqueness
