import ctl
import geometry

import math


class ReferenceFrame():
    '''
    This class describes a (local) reference frame for a molecule complex, e.g.
    a SymPlex.
    '''


    def __init__(self, molecule):
        '''
        Initialization of the ReferenceFrame class.
        '''
        self.reference_res = [None, None, None]
                # 0:center residue (pos 0.5)
                # 1:residue at rel. pos 0.25)
                # 2:residue at rel. pos 0.25)
        self.reference_points = [None, None, None]
        self.reference_points_type = None
        self.basis_vect = [None, None, None]

        self.molecule = molecule # molecule complex
        self.chimerax_session = self.molecule.chimerax_session

        self.get_reference_res()

        return


    def get_reference_res(self):
        '''
        Get reference residues as anchor points for the reference frame.
        '''
        center_res = self.molecule.get_center_res()
        center_res_25 = self.molecule.get_center_res(0.25)
        center_res_75 = self.molecule.get_center_res(0.75)

        self.reference_res = [center_res, center_res_25, center_res_75]

        return self.reference_res


    def get_reference_points(self):
        '''
        Get reference points as anchor points for the reference frame of a
        n-fold SymPlex with n>=3.
        '''
        ref_points = []

        if self.molecule.multimer_n < 3:
                ctl.error('ReferenceFrame: get_reference_points: ax1.fold < 3')

        elif self.molecule.multimer_n >= 3:
            center_res = self.reference_res[0]

            submodel_ids = [1, \
                            1+round(self.molecule.multimer_n*(1/3)), \
                            1+round(self.molecule.multimer_n*(2/3))]

            c0 = self.chimerax_session.get_xyz(\
                    (self.molecule.id[0], submodel_ids[0]), center_res)
            c1 = self.chimerax_session.get_xyz(\
                    (self.molecule.id[0], submodel_ids[1]), center_res)
            c2 = self.chimerax_session.get_xyz(\
                    (self.molecule.id[0], submodel_ids[2]), center_res)

            self.reference_points = [c0, c1, c2]
            self.reference_points_type = 0

        return self.reference_points


    def get_reference_points_2fold(self):
        '''
        Get reference points as anchor points for the reference frame of a
        2-fold SymPlex.
        '''
        ref_points = []
        quality = 0

        if self.molecule.multimer_n != 2:
                    ctl.error('ReferenceFrame: get_reference_points_2fold: '+ \
                              'multimer_n != 2')

        elif self.molecule.multimer_n == 2:
            center_res = self.reference_res
            submodel_ids = [1, 2]

            c0 = self.chimerax_session.get_xyz( \
                    (self.molecule.id[0], submodel_ids[0]), center_res[1])
            c1 = self.chimerax_session.get_xyz( \
                    (self.molecule.id[0], submodel_ids[1]), center_res[1])
            c2 = self.chimerax_session.get_xyz( \
                    (self.molecule.id[0], submodel_ids[0]), center_res[2])
            c3 = self.chimerax_session.get_xyz( \
                    (self.molecule.id[0], submodel_ids[1]), center_res[2])

            c0c1 = geometry.norm([c1[0]-c0[0], c1[1]-c0[1], c1[2]-c0[2]])
            c2c3 = geometry.norm([c3[0]-c2[0], c3[1]-c2[1], c3[2]-c2[2]])
            
            dp = geometry.dotproduct(c0c1, c2c3)
            dp_min = dp
            ref_points = [c0, c1, c2, c3]

            if abs(dp_min) > 0.5:
                c0 = self.chimerax_session.get_xyz( \
                        (self.molecule.id[0], submodel_ids[0]), \
                        center_res[0])
                c1 = self.chimerax_session.get_xyz( \
                        (self.molecule.id[0], submodel_ids[1]), \
                        center_res[0])
                c2 = self.chimerax_session.get_xyz( \
                        (self.molecule.id[0], submodel_ids[0]), \
                        center_res[2])
                c3 = self.chimerax_session.get_xyz( \
                        (self.molecule.id[0], submodel_ids[1]), \
                        center_res[2])

                c0c1 = geometry.norm([c1[0]-c0[0], c1[1]-c0[1], c1[2]-c0[2]])
                c2c3 = geometry.norm([c3[0]-c2[0], c3[1]-c2[1], c3[2]-c2[2]])
                
                dp = abs(geometry.dotproduct(c0c1, c2c3))
                if dp < dp_min:
                    dp_min = dp
                    ref_points = [c0, c1, c2, c3]

            if abs(dp_min) > 0.5:
                c0 = self.chimerax_session.get_xyz( \
                        (self.molecule.id[0], submodel_ids[0]), \
                        center_res[0])
                c1 = self.chimerax_session.get_xyz( \
                        (self.molecule.id[0], submodel_ids[1]), \
                        center_res[0])
                c2 = self.chimerax_session.get_xyz( \
                        (self.molecule.id[0], submodel_ids[0]), \
                        center_res[1])
                c3 = self.chimerax_session.get_xyz( \
                        (self.molecule.id[0], submodel_ids[1]), \
                        center_res[1])

                c0c1 = geometry.norm([c1[0]-c0[0], c1[1]-c0[1], c1[2]-c0[2]])
                c2c3 = geometry.norm([c3[0]-c2[0], c3[1]-c2[1], c3[2]-c2[2]])
                
                dp = abs(geometry.dotproduct(c0c1, c2c3))
                if dp < dp_min:
                    dp_min = dp
                    ref_points = [c0, c1, c2, c3]

            if abs(dp_min) > 0.8:
                quality = 0
            else:
                quality = 1

        return ref_points, quality


    def get_reference_points_2fold_zdist(self, type_given=None):
        '''
        Get reference points as anchor points for the reference frame of a
        2-fold SymPlex.
        The reference points must have a minimum distance on the z axis.
        '''
        ref_points = []
        submodel_ids = [1, 2]
        ref_res = self.reference_res

        quality = 0
        type_ = 2

        if self.molecule.multimer_n != 2:
                    ctl.error('ReferenceFrame: get_reference_points_2fold: '+ \
                              'multimer_n != 2')

        elif self.molecule.multimer_n == 2:

            # reference points of type 2
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~
            c0 = self.chimerax_session.get_xyz( \
                    (self.molecule.id[0], submodel_ids[0]), ref_res[1])
            c1 = self.chimerax_session.get_xyz( \
                    (self.molecule.id[0], submodel_ids[1]), ref_res[1])
            c2 = self.chimerax_session.get_xyz( \
                    (self.molecule.id[0], submodel_ids[0]), ref_res[2])
            c3 = self.chimerax_session.get_xyz( \
                    (self.molecule.id[0], submodel_ids[1]), ref_res[2])

            dz_max = abs(c0[2]-c2[2])
            ref_points = [c0, c1, c2, c3]

            # reference points of type 3
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~
            c0 = self.chimerax_session.get_xyz( \
                    (self.molecule.id[0], submodel_ids[0]), ref_res[0])
            c1 = self.chimerax_session.get_xyz( \
                    (self.molecule.id[0], submodel_ids[1]), ref_res[0])
            c2 = self.chimerax_session.get_xyz( \
                    (self.molecule.id[0], submodel_ids[0]), ref_res[2])
            c3 = self.chimerax_session.get_xyz( \
                    (self.molecule.id[0], submodel_ids[1]), ref_res[2])

            if (abs(c0[2]-c2[2]) > dz_max and type_given == None) or \
               type_given == 3:
                type_ = 3
                dz_max = abs(c0[2]-c2[2])
                ref_points = [c0, c1, c2, c3]

            # reference points of type 4
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~
            c0 = self.chimerax_session.get_xyz( \
                    (self.molecule.id[0], submodel_ids[0]), ref_res[0])
            c1 = self.chimerax_session.get_xyz( \
                    (self.molecule.id[0], submodel_ids[1]), ref_res[0])
            c2 = self.chimerax_session.get_xyz( \
                    (self.molecule.id[0], submodel_ids[0]), ref_res[1])
            c3 = self.chimerax_session.get_xyz( \
                    (self.molecule.id[0], submodel_ids[1]), ref_res[1])

            if (abs(c0[2]-c2[2]) > dz_max and type_given == None) or \
               type_given == 4:
                type_ = 4
                ref_points = [c0, c1, c2, c3]

            if abs(ref_points[0][2]-ref_points[2][2]) >= 8:
                quality = 1

        if type_given != None:
            if type_given != type_:
                ctl.e(type_given)
                ctl.e(type_)
                ctl.error('ReferenceFrame: '+ \
                          'get_reference_points_2fold_zdist: '+ \
                          'type_given != type_')

        return ref_points, quality, type_
