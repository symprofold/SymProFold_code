import ctl
import geometry

import math
from structure.primitive_unit_cell import PrimitiveUnitCell


class Layer():
    '''
    This class describes a layer. A layer is build from SymPlexes (axes).
    '''


    def __init__(self, axes=[]):
        ''' Initialization of the Axis class. '''

        self.symmgroup = ''
        self.symmgroups_compatible = []
        self.representations = {}
        self.axes = []
        self.primitive_unit_cell = PrimitiveUnitCell(self)
        self.lattice_constant = 0
        self.refpoint_constant = 0
        self.lattice_constant_raw = ''
        self.refpoint_constant_raw = ''
        self.flipped = False
            # indication if the layer has been flipped
            # This is important for e.g. reference points, which then also have
            # to be determined for the flipped case.

        self.add_axes(axes)

        return
        

    def set_symmgroup(self, symmgroup):
        ''' Set symmetry group. '''
        self.symmgroup = symmgroup
        return


    def add_representation(self, representation):
        ''' Add representation. '''
        self.representations[representation.id] = representation
        return


    def get_representation(self, model_id):
        ''' Get representation. '''

        model_reg = self.axes[0].model_reg
        model_id = model_reg.convert_model_id(model_id)

        representation = self.representations[model_id]

        return representation


    def get_representations(self):
        ''' Get representations. '''

        model_reg = self.axes[0].model_reg

        representations = {}
        for r in self.representations:
            representations[r] = self.representations[r]

        return representations


    def add_axis(self, axis):
        ''' Add axis. '''
        self.axes.append(axis)
        return


    def add_axes(self, axes):
        ''' Add list of axes. '''

        for axis in axes:
            self.axes.append(axis)

        return


    def ax_models(self, ax):
        '''
        Get all representants of a given Axis object ax that are assigned to
        this layer.
        '''
        ax_model_ids = sorted(list( \
                set(self.get_representations().keys()) & \
                set(ax.get_representations().keys()) ))

        ax_models_ = []
        for model_id in ax_model_ids:
            ax_models_.append(self.get_representation(model_id))

        return ax_models_


    def calc_refpoint_constant(self, lc_offset=0):
        ''' Calculate reference point constant. '''

        model_n = 0
        distsum = 0
        distances = ''

        for i,model in enumerate(self.ax_models(self.axes[0])):
            if i == 0:
                continue

            model_n += 1   
            model_xyz = model.trans_vect

            d = geometry.dist(model_xyz, [0, 0, 0])
            distances = distances+str(round(d/10, 3))+'nm, '
            distsum += d


        if model_n == 0:
            self.refpoint_constant = 0
            self.refpoint_constant_raw = ''
        else:
            self.refpoint_constant = distsum/model_n + lc_offset
            self.refpoint_constant_raw = distances[:-2]

        return self.refpoint_constant


    def calc_lattice_constant(self, lc_offset=0):
        ''' Calculate lattice constant from reference point constant. '''

        model_n = 0
        distsum = 0
        distances = ''

        ax0_models = self.ax_models(self.axes[0])

        for i,model in enumerate(ax0_models):
            if i == 0:
                model_neighbour_xyz = ax0_models[-1].trans_vect
                continue

            model_n += 1   
            model_xyz = model.trans_vect

            # symmetry group p6 with 3fold axis0 and 2fold axis1
            if self.axes[0].fold == 3 and self.axes[1].fold == 2:
                d = geometry.dist(model_xyz, model_neighbour_xyz)
            else:
                d = geometry.dist(model_xyz, [0, 0, 0])

            distances = distances+str(round(d/10, 3))+'nm, '
            distsum += d

            model_neighbour_xyz = model_xyz


        if model_n == 0:
            self.lattice_constant = 0
            self.lattice_constant_raw = ''
        else:
            self.lattice_constant = distsum/model_n + lc_offset
            self.lattice_constant_raw = distances[:-2]

        return self.lattice_constant


    def seqlen(self):
        '''
        Determine sequence length of combined monomer, without N terminus and
        C terminus.
        '''
        seq_range = [-1, -1]

        ax0_repr = self.axes[0].get_representations()[(1,)]

        if len(self.axes) > 1:
            ax1_repr = self.axes[1].get_representations()[(2,)]
            seq_range[0] = min(ax0_repr.termini[0], ax1_repr.termini[0])
            seq_range[1] = max(ax0_repr.termini[1], ax1_repr.termini[1])
        else:
            seq_range[0] = ax0_repr.termini[0]
            seq_range[1] = ax0_repr.termini[1]

        seq_len = seq_range[1]-seq_range[0]+1

        return seq_len


    def export_param(self, path, lc_offset=0):
        ''' Export parameter to biochem_properties file. '''

        ctl.d('export_param')
        ctl.d(self.get_representations())
        self.calc_lattice_constant(lc_offset)

        # area of unit cell
        area = 0
        alpha = 90
        area_str = 'not determined'
        if self.symmgroup == 'p3':
            area = math.sqrt(3)/2*(self.lattice_constant**2)
        if self.symmgroup == 'p4':
            area = self.lattice_constant**2
        if self.symmgroup == 'p6':
            area = math.sqrt(3)/2*(self.lattice_constant**2)

        if area != 0:
            area_str = str(int(round(area/100, 0)))+'nm^2'

        seqlen = self.seqlen()

        f = open(path+'_biochem_properties.txt', 'w')

        # symmetry group
        symmgroup_str = 'not determined'
        if self.symmgroup != '':
            symmgroup_str = self.symmgroup
        f.write('symmetry group:   '+symmgroup_str+"\r\n")

        # lattice constant
        if self.lattice_constant != 0:
            f.write('a = '+str(round(self.lattice_constant/10,1))+'nm'+ \
                    "   (calculated using ab initio models)\r\n")

            relation = min(self.lattice_constant, self.refpoint_constant) / \
                       max(self.lattice_constant, self.refpoint_constant)
            if relation > 0.99:
                f.write("\r\n"+'lattice constant calculated as the average '+ \
                        "of the following values:\r\n")
                f.write(self.lattice_constant_raw+"\r\n")

        # area of unit cell
        f.write("\r\n"+'area of unit cell:   '+area_str+"\r\n")  

        # mol per unit cell
        mol_per_unit_cell = 'not determined'
        if symmgroup_str != 'not determined':
            mol_per_unit_cell = self.axes[0].fold

        f.write(f'molecules / unit cell:   {mol_per_unit_cell}'+"\r\n")
        f.write('sequence length (monomer, without n-term and c-term):   '+ \
                str(seqlen)+" residues\r\n")  

        # density
        if area != 0:
            density = (self.axes[0].fold*seqlen)/area
            f.write(f"residues per nm^2:   {round(density*100, 1)}\r\n")

        f.close()

        return


    def get_possible_symmgroups(self, axes):
        '''
        Determine possible symmetry groups using fold information of given
        axes.

        axes: list of the axis' respective folds
        '''
        input_folds = [(0), (0), 0, 0, 0, (0), 0]
        sg_str = ['', 'p1', 'p2', 'p3', 'p4', '', 'p6']

        compatible = []
        determined = []

        for a in axes:
            input_folds[a] += 1

        # check for p1
        if input_folds[2] == 0 and \
           input_folds[3] == 0 and \
           input_folds[4] == 0 and \
           input_folds[6] == 0:

           compatible.append(1)

        # check for p2
        if input_folds[2] in (1, 2, 3) and \
           input_folds[3] == 0 and \
           input_folds[4] == 0 and \
           input_folds[6] == 0:

           compatible.append(2)

        if input_folds[2] == 3 and \
           input_folds[3] == 0 and \
           input_folds[4] == 0 and \
           input_folds[6] == 0:

           determined.append(2)

        # check for p3
        if input_folds[2] == 0 and \
           input_folds[3] in (1, 2, 3) and \
           input_folds[4] == 0 and \
           input_folds[6] == 0:

           compatible.append(3)

        if input_folds[2] == 0 and \
           input_folds[3] in (2, 3) and \
           input_folds[4] == 0 and \
           input_folds[6] == 0:

           determined.append(3)

        # check for p4
        if input_folds[2] in (0, 1) and \
           input_folds[3] == 0 and \
           input_folds[4] in (0, 1, 2) and \
           input_folds[6] == 0 and \
          (input_folds[2]+input_folds[4]) in (1, 2, 3):

            compatible.append(4)

        if input_folds[2] in (0, 1) and \
           input_folds[3] == 0 and \
           input_folds[4] in (1, 2) and \
           input_folds[6] == 0 and \
          (input_folds[2]+input_folds[4]) in (2, 3):

            determined.append(4)

        # check for p6
        if input_folds[2] in (0, 1) and \
           input_folds[3] in (0, 1) and \
           input_folds[4] == 0 and \
           input_folds[6] in (0, 1) and \
          (input_folds[2]+input_folds[3]+input_folds[6]) in (1, 2, 3):

           compatible.append(6)

        if input_folds[2] in (0, 1) and \
           input_folds[3] in (0, 1) and \
           input_folds[4] == 0 and \
           input_folds[6] in (0, 1) and \
          (input_folds[2]+input_folds[3]+input_folds[6]) in (2, 3):

            determined.append(6)


        if len(determined) > 1:
            ctl.error('get_possible_symmgroups: '+ \
                      +'more than 1 determined symmgroup')


        if len(determined) == 1:
            self.symmgroup = sg_str[determined[0]]

        self.symmgroups_compatible = compatible

        ret = [compatible, determined]

        return ret


    def export_possible_symmgroups(self, axes, conf):
        ''' Export possible symmetry groups to text file. '''

        sg_str = ['', 'p1', 'p2', 'p3', 'p4', '', 'p6']
        sg = self.get_possible_symmgroups(axes)

        if len(sg[1]) == 1:
            f = open(conf.export_path+'spacegroup_predicted_'+ \
                     sg_str[sg[1][0]]+'.txt', 'w')
            f.write('')
            f.close()

        if len(sg[1]) > 1:
            f = open(conf.export_path+ \
                     'spacegroup_determined_ERROR_more_than_1'+'.txt', 'w')
            f.write('')
            f.close()

        if len(sg[0]) > 0 and len(sg[1]) != 1: 
            sglist = ''
            for s in sg[0]:
                sglist = sglist+str(sg_str[s])+'_'
            
            f = open(conf.export_path+'spacegroups_possible_'+sglist[:-1]+ \
                     '.txt', 'w')
            f.write('')
            f.close()

        if len(sg[1]) >= 1 and (sg[0] != sg[1]):
            ctl.d(sg)
            f = open(conf.export_path+'spacegroup_ERROR_CODE01.txt', 'w')
            f.write(str(sg[0])+"\r\n")
            f.write(str(sg[1])+"\r\n")
            f.close()

        return     


    def orient_points_primitive_unit_cell(self):
        '''
        Get orientation points of primitive unit cell.
        '''

        lc = [self.refpoint_constant]
        p = []

        if self.symmgroup == 'p1' or self.symmgroup == 'p2':
            lc = [self.refpoint_constant, self.refpoint_constant]
            al = 90
            da = math.cos(math.radians(al))*lc[1]
            db = math.sin(math.radians(al))*lc[0]
            p.append( [        0,   0, 0] )
            p.append( [    lc[0],   0, 0] )

            p.append( [     0+da,  db, 0] )
            p.append( [ lc[0]+da,  db, 0] )

        if self.symmgroup == 'p3' or self.symmgroup == 'p6':
            d = math.sqrt(3)/2

            p.append( [        0,          0, 0] )
            p.append( [  d*lc[0],   -lc[0]/2, 0] )
            p.append( [  d*lc[0],    lc[0]/2, 0] )
            p.append( [2*d*lc[0],          0, 0] )

        if self.symmgroup == 'p4':
            p.append( [     0,      0, 0] )
            p.append( [ lc[0],      0, 0] )
            p.append( [ lc[0],  lc[0], 0] )
            p.append( [ lc[0],  lc[0], 0] )


        if self.symmgroup == '':

            symmgroup = ''
            for sc in self.symmgroups_compatible:
                if sc == 2:
                    symmgroup = 'p2'

            if symmgroup == 'p2':
                lc = [self.refpoint_constant]
                p.append( [-lc[0],      0, 0] )
                p.append( [     0,      0, 0] )
                p.append( [ lc[0],      0, 0] )

        return p


    def orient_points(self):
        '''
        Get orientation points for axis 0 aligned to primitive unit cell.
        '''
        lc = [self.refpoint_constant]
        p = []

        if self.symmgroup == 'p1' or self.symmgroup == 'p2':
            lc = [self.refpoint_constant, self.refpoint_constant]
            al = 90
            da = math.cos(math.radians(al))*lc[1]
            db = math.sin(math.radians(al))*lc[0]
            p.append( [-lc[0]-da, -db, 0] )
            p.append( [     0-da, -db, 0] )
            p.append( [ lc[0]-da, -db, 0] )

            # p.append( [-lc[0],      0, 0] )
            p.append( [     0,      0, 0] )
            # p.append( [ lc[0],      0, 0] )

            # p.append( [-lc[0]+da,  db, 0] )
            p.append( [     0+da,  db, 0] )
            # p.append( [ lc[0]+da,  db, 0] )


        if self.symmgroup == 'p3':
            d = math.sqrt(3)/2
            
            # p.append( [       0,   -lc[0]  , 0] )

            p.append( [-d*lc[0],   -lc[0]/2, 0] )
            p.append( [ d*lc[0],   -lc[0]/2, 0] )

            p.append( [       0,          0, 0] )

            # p.append( [-d*lc[0],    lc[0]/2, 0] )
            # p.append( [ d*lc[0],    lc[0]/2, 0] )

            p.append( [       0,    lc[0]  , 0] )


        if self.symmgroup == 'p4':
            p.append( [     0, -lc[0], 0] )

            p.append( [-lc[0],      0, 0] )
            p.append( [     0,      0, 0] )
            p.append( [ lc[0],      0, 0] )

            p.append( [     0,  lc[0], 0] )


        if self.symmgroup == 'p6':

            if self.axes[0].fold == 6:

                d = math.sqrt(3)/2
                
                p.append( [       0,   -lc[0]  , 0] )

                p.append( [-d*lc[0],   -lc[0]/2, 0] )
                p.append( [ d*lc[0],   -lc[0]/2, 0] )

                p.append( [       0,          0, 0] )

                p.append( [-d*lc[0],    lc[0]/2, 0] )
                p.append( [ d*lc[0],    lc[0]/2, 0] )

                p.append( [       0,    lc[0]  , 0] )

            elif self.axes[0].fold == 3:
                d = math.sqrt(3)/2
                
                # p.append( [       0,   -lc[0]  , 0] )

                p.append( [-d*lc[0],   -lc[0]/2, 0] )
                p.append( [ d*lc[0],   -lc[0]/2, 0] )

                p.append( [       0,          0, 0] )

                # p.append( [-d*lc[0],    lc[0]/2, 0] )
                # p.append( [ d*lc[0],    lc[0]/2, 0] )

                p.append( [       0,    lc[0]  , 0] )


        if self.symmgroup == '':

            symmgroup = ''
            for sc in self.symmgroups_compatible:
                if sc == 2:
                    symmgroup = 'p2'

            if symmgroup == 'p2':
                lc = [self.refpoint_constant]
                p.append( [-lc[0],      0, 0] )
                p.append( [     0,      0, 0] )
                p.append( [ lc[0],      0, 0] )


        if self.flipped == True:
            for p_ in p:
                p_[1] *= (-1)

        return p


    def orient_points_2fold(self):
        '''
        Get orientation points for axis 1 (2-fold axis) aligned to primitive
        unit cell.
        '''
        lc = [self.refpoint_constant]
        p = [ ]

        if self.symmgroup == 'p4':
            p.append( [       0, -lc[0]/2, 0] )

            p.append( [-lc[0]/2,        0, 0] )
            p.append( [ lc[0]/2,        0, 0] )

            p.append( [       0,  lc[0]/2, 0] )


        if self.symmgroup == 'p6':
            d = math.sqrt(3)/2
            
            p.append( [         0,   -lc[0]  /2, 0] )

            p.append( [-d*lc[0]/2,   -lc[0]/2/2, 0] )
            p.append( [ d*lc[0]/2,   -lc[0]/2/2, 0] )

            p.append( [-d*lc[0]/2,    lc[0]/2/2, 0] )
            p.append( [ d*lc[0]/2,    lc[0]/2/2, 0] )

            p.append( [         0,    lc[0]  /2, 0] )


        if self.symmgroup == '':

            symmgroup = ''
            for sc in self.symmgroups_compatible:
                if sc == 2:
                    symmgroup = 'p2'

            if symmgroup == 'p2':
                lc = [self.refpoint_constant]
                p.append( [-lc[0]/2,      0, 0] )
                p.append( [ lc[0]/2,      0, 0] )

        return p


    def orient_points_3fold(self):
        '''
        Get orientation points for axis 1 (3-fold axis) aligned to primitive
        unit cell.
        '''
        lc = [self.refpoint_constant]
        p = [ ]

        if self.symmgroup == 'p3':
            d = math.sqrt(3)/2
            r = 1/math.sqrt(3)

            p.append( [-r*lc[0]/2, -lc[0]/2, 0] )
            p.append( [ r*lc[0]/2, -lc[0]/2, 0] )

            p.append( [-r*lc[0]  ,        0, 0] )
            p.append( [ r*lc[0]  ,        0, 0] )

            p.append( [-r*lc[0]/2,  lc[0]/2, 0] )
            p.append( [ r*lc[0]/2,  lc[0]/2, 0] )

        elif self.symmgroup == 'p6':
            d = math.sqrt(3)/2
            r = 1/math.sqrt(3)

            p.append( [-r*lc[0]/2, -lc[0]/2, 0] )
            p.append( [ r*lc[0]/2, -lc[0]/2, 0] )

            p.append( [-r*lc[0]  ,        0, 0] )
            p.append( [ r*lc[0]  ,        0, 0] )

            p.append( [-r*lc[0]/2,  lc[0]/2, 0] )
            p.append( [ r*lc[0]/2,  lc[0]/2, 0] )

        return p


    def orient_points_4fold(self):
        '''
        Get orientation points for axis 1 (4-fold axis) aligned to primitive
        unit cell.
        '''
        lc = [self.refpoint_constant]
        p = [ ]

        if self.symmgroup == 'p4':
            p.append( [-lc[0]/2, -lc[0]/2, 0] )
            p.append( [ lc[0]/2, -lc[0]/2, 0] )
            p.append( [-lc[0]/2,  lc[0]/2, 0] )
            p.append( [ lc[0]/2,  lc[0]/2, 0] )

        return p


    def get_ref_point(self, coords):
        ''' Get nearest reference point (for axis 0) to given coordinates. '''

        ref_points = self.orient_points()
        cand_min = [0,0,0]
        dist_min = 10000

        for r in ref_points:
            d = geometry.dist(coords, r)
            if d < dist_min:
                cand_min = r
                dist_min = d

        return cand_min, dist_min


    def get_ref_point_ax1(self, coords, axis):
        ''' Get nearest reference point (for axis 1) to given coordinates. '''

        if axis.fold == 2:
            rp, d = self.get_ref_point_2fold(coords)
            return rp, d

        if axis.fold == 3:
            rp, d = self.get_ref_point_3fold(coords)
            return rp, d

        if axis.fold == 4:
            rp, d = self.get_ref_point_4fold(coords)
            return rp, d

        ctl.e(axis.fold)
        ctl.error('Layer: get_ref_point_ax1: no ref point available.')

        return


    def get_ref_point_2fold(self, co):
        '''
        Get nearest reference point (for axis 1 (2), 2-fold axis) to given
        coordinates.
        '''
        ref_points = self.orient_points_2fold()
        cand_min = [0, 0, 0]
        dist_min = 10000

        for r in ref_points:
            d = geometry.dist(co, r)
            if d < dist_min:
                cand_min = r
                dist_min = d

        return cand_min, dist_min


    def get_ref_point_3fold(self, co):
        '''
        Get nearest reference point (for axis 1 (2), 3-fold axis) to given
        coordinates.
        '''
        ref_points = self.orient_points_3fold()
        cand_min = [0, 0, 0]
        dist_min = 10000

        for r in ref_points:
            d = geometry.dist(co, r)
            if d < dist_min:
                cand_min = r
                dist_min = d

        return cand_min, dist_min


    def get_ref_point_4fold(self, co):
        '''
        Get nearest reference point (for axis 1 (2), 4-fold axis) to given
        coordinates.
        '''
        ref_points = self.orient_points_4fold()
        cand_min = [0, 0, 0]
        dist_min = 10000

        for r in ref_points:
            d = geometry.dist(co, r)
            if d < dist_min:
                cand_min = r
                dist_min = d

        return cand_min, dist_min
