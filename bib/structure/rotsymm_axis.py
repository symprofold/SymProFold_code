import ctl
import geometry
import math


class RotSymmAxis():
    '''
    This class describes the rotational symmetry axis of a SymPlex.
    '''


    def __init__(self, molecule):
        '''
        Initialization of the RotSymmAxis class.
        '''
        self.molecule = molecule # molecule complex
        self.frame = self.molecule.frame

        self.axis = []
                # rotational symmetry axis when get_rotsymm_axis() is called
        self.tilt = None
                # axial tilt in degrees
        self.orientation_initial = None
                # sign that was needed to orient rotsymm_axis of a SymPlex with
                # 2-fold symmetry in a way that the z component is >0

        self.get_rotsymm_axis()


        # check if rotational symmetry axis is perpendicular to xy plane
        if abs(self.axis[2]) < 0.99:
            ctl.e(self.molecule.model_id)
            ctl.e(self.axis)
            raise Exception('RotSymmAxis: init: '+ \
                      'rotational symmetry axis not perpendicular to xy plane')

        # check if z component of rotational symmetry axis is >0.
        if self.axis[2] <= 0:
            ctl.e(self.molecule.model_id)
            ctl.e(self.axis)
            raise Exception('RotSymmAxis: init: '+ \
                      'z component of rotational symmetry axis is <=0')

        return


    def get_tilt(self):
        '''
        Get axial tilt relative to the z axis (in rad).
        '''
        tilt = math.acos(self.axis[2])
        self.tilt = tilt

        return tilt


    def get_rotsymm_axis(self, mode_readout_angle_change=False):
        '''
        Get rotational symmetry axis.

        Args:
            mode_readout_angle_change: mode to readout the current rotsymm_axis
                    for calculation of the angle change compared to last call
                    of get_rotsymm_axis()
        '''
        if self.molecule.multimer_n <= 1:
                    ctl.error('RotSymmAxis: get_rotsymm_axis: multimer_n <= 1')

        elif self.molecule.multimer_n <= 2:

            # standard mode
            if not mode_readout_angle_change:
                rp0, quality0 = self.frame.get_reference_points_2fold()
                rp1, quality1, type_ = self.frame. \
                                            get_reference_points_2fold_zdist()

                if quality0 == 1:
                    rp = rp0
                    basisvect = [None, None, None]
                    basisvect[0] = geometry.norm(rp[1]-rp[0])
                    basisvect[1] = geometry.norm(rp[3]-rp[2])
                    basisvect[2] = geometry.norm(geometry.\
                                    crossproduct(basisvect[0], basisvect[1]))

                    rotaxis = basisvect[2]
                    self.frame.reference_points_type = 1

                elif quality1 == 1:
                    rp = rp1
                    vect = [None, None]
                    vect[0] = rp[1]-rp[0]
                    vect[1] = rp[3]-rp[2]

                    rotaxis = [(rp[2][i]+vect[1][i]*0.5) - \
                               (rp[0][i]+vect[0][i]*0.5) for i in range(0,3)]
                    rotaxis = geometry.norm(rotaxis)
                    self.frame.reference_points_type = type_

                else:
                    raise Exception('RotSymmAxis: get_rotsymm_axis: '+ \
                            'determination of rotational symmetry axis '+ \
                            'not possible')

                if rotaxis[2] < 0:
                    self.orientation_initial = -1
                    rotaxis = [self.orientation_initial*x for x in rotaxis]
                else:
                    self.orientation_initial = 1

            # mode_readout_angle_change
            else:
                if self.frame.reference_points_type == 1:
                    rp, zero = self.frame.get_reference_points_2fold()

                    basisvect = [None, None, None]
                    basisvect[0] = geometry.norm(rp[1]-rp[0])
                    basisvect[1] = geometry.norm(rp[3]-rp[2])
                    basisvect[2] = geometry.norm(geometry.\
                                    crossproduct(basisvect[0], basisvect[1]))
                    rotaxis = basisvect[2]

                elif self.frame.reference_points_type >= 2:
                    rp, zero, zero = self.frame. \
                                    get_reference_points_2fold_zdist( \
                                        self.frame.reference_points_type)

                    vect = [None, None]
                    vect[0] = rp[1]-rp[0]
                    vect[1] = rp[3]-rp[2]

                    rotaxis = [(rp[2][i]+vect[1][i]*0.5) - \
                               (rp[0][i]+vect[0][i]*0.5) for i in range(0,3)]
                    rotaxis = geometry.norm(rotaxis)

                else:
                    ctl.error('RotSymmAxis: get_rotsymm_axis: '+ \
                              'mode_readout_angle_change, '+ \
                              'determination of rotational symmetry axis '+ \
                              'not possible')

                rotaxis = [self.orientation_initial*x for x in rotaxis]

        elif self.molecule.multimer_n >= 3:
            rp = self.frame.get_reference_points()

            if self.frame.reference_points_type != 0:
                ctl.error('RotSymmAxis: get_rotsymm_axis: '+ \
                          'self.frame.reference_points_type != 0')

            basisvect = [None, None, None]
            basisvect[0] = geometry.norm(rp[1]-rp[0])
            basisvect[1] = geometry.norm(rp[2]-rp[0])
            basisvect[2] = geometry.norm(geometry. \
                                    crossproduct(basisvect[0], basisvect[1]))
            rotaxis = basisvect[2]

        self.axis = rotaxis

        return rotaxis
