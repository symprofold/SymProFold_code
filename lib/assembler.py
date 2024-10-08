import structure.complex_
from structure.layer import Layer
import structure.monomer

import bib
import align_layer
import chimerax_api
import complete_outer_chains
import ctl
import export
import filesystem
import geometry
import math
import os
import structure.primitive_unit_cell
import proteinmerge
import snapin
import validation

import glob
import json


class Assembler():
    '''
    This class builds the layer model from the parts.
    '''


    def __init__(self, axes, conf, export_file_infix='', \
                 contact_submodel_orientation=1, \
                 primitive_unit_cell_molecules=[], \
                 lc_offset=False):
        ''' Initialization of the Assembler class. '''

        self.axes = axes

        layer = Layer(self.axes)
        layer_flat = Layer(self.axes)

        if primitive_unit_cell_molecules != []:
            for primitive_unit_cell_molecule in primitive_unit_cell_molecules:
                layer_flat.primitive_unit_cell.add_molecule( \
                                                primitive_unit_cell_molecule)

        self.layers = [layer, layer_flat]

        self.conf = conf
        self.alignment_pivot_pos = 1
            # 0: pivot residue at the upstream side of the alignment section
            #    (e.g. domain)
            # 1: pivot residue at the downstream side of the alignment section
            #    (e.g. domain)

        self.export_file_infix = export_file_infix

        self.contact_submodel_orientation = contact_submodel_orientation
                # relative position of second ax0 representation/representative
                # regarding ax1 representation/representative

        self.model_reg = self.axes[0].model_reg

        self.lc_offset = lc_offset
        self.lc_var = [0, 30, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, \
                       22, 24, 26, 28]
        # self.lc_var = [0, 2, 4]

        if lc_offset == False:
            self.lc_var = [0]

        # test if gene contains only one domain
        if len(conf.domains) == 1 and len(self.axes) > 1:
            conf.set_export_path_filters([
                '../raw/'+conf.export_file_prefix+'_showfullsuperposition', \
                '../'+conf.export_file_prefix, \
                '../raw/'+conf.export_file_prefix+'_show'+ \
                str(self.axes[0].fold)+'foldax'])

        return


    def run(self, debug_mode=False):
        '''
        Run whole assembly process.

        Return:
            int: 1:assembly possible
                 2:assembly not possible because tilt >45°
                 3:assembly not possible because gap between both SymPlexes
                    too large
        '''
        folds = [ax.fold for ax in self.axes]
        self.layers[1].export_possible_symmgroups(folds, self.conf)

        try:
            self.build_layer()
        except Exception as e:
            if str(e) == 'build_layer: tilt of ax1 not within range':
                ctl.p(str(e))
                return False, 2, str(e)

            elif str(e) == 'RotSymmAxis: init: '+ \
                            'rotational symmetry axis not perpendicular to '+ \
                            'xy plane':
                ctl.p(str(e))
                return False, 4, str(e)

            elif str(e) == 'RotSymmAxis: get_rotsymm_axis: '+ \
                            'determination of rotational symmetry axis '+ \
                            'not possible':
                ctl.p(str(e))
                return False, 6, str(e)

            elif str(e) == 'RotSymmAxis: init: '+ \
                            'z component of rotational symmetry axis is <=0':
                ctl.p(str(e))
                return False, 7, str(e)

            elif str(e) == 'process_layer_lc: no symmetry group determined':
                ctl.p(str(e))
                return False, 9, str(e)

            elif str(e) == 'Assembler: alignment_pivot_res: '+ \
                            'no coincident residues':
                ctl.p(str(e))
                return False, 12, str(e)

            else:
                ctl.error('run: Exception')

        self.layers[1].export_param(self.conf.export_path+ \
                                    self.conf.export_file_prefix+ \
                                    self.export_file_infix, \
                                    lc_offset=self.lc_offset)

        try:
            clashes = self.process_layer()
        except Exception as e:
            if str(e) == 'combine_chains: dist_check failed':
                ctl.p(str(e))
                return False, 3, str(e)

            elif str(e) == 'snapin_layer: snapin distance too large':
                ctl.p(str(e))
                return False, 5, str(e)

            elif str(e) == 'snapin_layer: snapin rotation too large':
                ctl.p(str(e))
                return False, 8, str(e)

            elif str(e) == 'snapin_layer: snapin rotation step too large':
                ctl.p(str(e))
                return False, 10, str(e)

            elif str(e) == 'process_layer_specific: '+ \
                            'ax surface completely within termini':
                ctl.p(str(e))
                return False, 11, str(e)

            elif str(e) == 'process_layer_specific: '+ \
                            'monomers not complete':
                ctl.p(str(e))
                return False, 13, str(e)

            else:
                ctl.error('run: Exception')

        ctl.d(clashes)

        if not debug_mode:
            self.conf.delete_rawfiles()

        return True, 1, 'run: finished'


    def get_coincident_residues(self, id1, id2, res_ov, session, max_dist=1):
        '''
        Get corresponding residues in residue overlap that have a distance <=1.
        '''
        res_d = []
        coin_res = []

        for r in res_ov:
            id1_r = session.get_xyz(id1, r)
            id2_r = session.get_xyz(id2, r)

            d = geometry.dist(id1_r, id2_r)
            res_d.append(d)

            if d <= max_dist:
                coin_res.append(r)

        return coin_res, res_d


    def ax_a_model_add_ax_b(self, ax_a, ax_b, ax_a_model_id, ax_b_model_id, \
                            mode='default', set_connection=True):
        '''
        Match a model of ax_b to a model of ax_a by superposition and create a
        connection.
        '''
        ax_a_model_id = self.model_reg.convert_model_id(ax_a_model_id)
        ax_b_model_id = self.model_reg.convert_model_id(ax_b_model_id)

        ax_a_id = ax_a_model_id[0]
        ax_b_id = ax_b_model_id[0]

        ax_a.chimerax_session.match(ax_b_model_id, \
                                    (ax_b.domains[0][0], ax_b.domains[0][1]), \
                                    ax_a_model_id, \
                                    ax_b_id)

        if set_connection == True:
            ax_a.get_representation(ax_a_id).set_connection(ax_a_model_id, \
                                                         ax_b_model_id, mode)

        return


    def ax1_contact_submodel_0(self, binding_site, ax1fold):
        '''
        Determine axis 1 submodel_id when connected to axis 0.

        binding_site: binding_site of axis0, starting with 1
        '''
        submodel_id = (binding_site-1)%ax1fold+1

        return submodel_id


    def ax1_contact_submodel_1(self, binding_site, orientation, ax0fold):
        '''
        Determine axis 1 submodel_id when further axis 0 is attached.
        '''
        submodel_id = (ax0fold+ \
                       self.ax1_contact_submodel_0(binding_site, ax0fold)+ \
                       orientation-1)%ax0fold+1

        return submodel_id


    def build_layer(self):
        '''
        Place axis representations/representatives so that they form the layer.

        4 main steps:
        - step 1: place central ax0
        - step 2: place ax1 models around central ax
        - step 3: place ax0 representations/representatives around ax1
        - step 4: place ax0 models flattened around ax1
        '''
        if self.layers[1].symmgroup == '' and \
           self.layers[1].symmgroups_compatible == []:
            raise Exception('process_layer_lc: no symmetry group determined')

        sess = self.axes[0].chimerax_session
        sess.init()

        ax0 = self.axes[0]

        ax0_pairing_offset = 0
        if ax0.fold%2 == 0:
            ax0_pairing_offset = int(ax0.fold/2)
        elif ax0.fold == 3:
            ax0_pairing_offset = 2


        # step 1: place central ax0
        # -------------------------     
        ax0_current_model = ax0.open_model()

        for l in self.layers:
            l.add_representation(ax0_current_model)

        # set name to submodels of ax0 representation/representative
        for i in range(1, 1+ax0.fold):
            sess.rename_model((ax0_current_model.id[0], i), \
                    'ax'+ax0_current_model.idstr+'mol'+str(i))


        # continue, if more than 1 axis object is provided as input
        if len(self.axes) > 1:
            ax1 = self.axes[1]

            layer = self.layers[0]
            layer_flat = self.layers[1]

            # no influence on general assembly, but better orientation of
            # ax1 representations/representatives in first steps
            if ax1.fold > 2:
                orientation_offset = self.contact_submodel_orientation*(-1)

            if ax1.fold == 2:
                orientation_offset = 0


            # step 2: place ax1 models around central ax0
            # -------------------------------------------
            for ax0_submodel_id in range(1, 1+ax0.fold):
                        # 1 because first ax1 element has model id .1

                ax1_current_model = ax1.open_model()
                layer.add_representation(ax1_current_model)
                layer_flat.add_representation(ax1_current_model)

                renamed = False # avoid 2x renaming of chain ids

                if ax1.model_active_orient == -1:
                    intermediate_ids = []
                    ax1_submodel_id_fixed = 1
                    
                    for ii in range(1, ax1.fold):

                        # inverse chain ids for flipped symplex
                        # ax1_submodel_id is fixed
                        sess.change_model_id( \
                                (ax1_current_model.id[0], \
                                    (ax1_submodel_id_fixed+ii-1)%ax1.fold+1), \
                                (ax1_current_model.id[0], \
                                    10+(ax1.fold+ax1_submodel_id_fixed-ii-1)% \
                                    ax1.fold+1) )
                        intermediate_ids.append( \
                                    10+(ax1.fold+ax1_submodel_id_fixed-ii-1)% \
                                    ax1.fold+1)

                    for ii in intermediate_ids:
                        sess.change_model_id( \
                                (ax1_current_model.id[0], ii), \
                                (ax1_current_model.id[0], (ii-10)) )

                    renamed = True # avoid 2x renaming of chain ids


                # determine rotational symmetry axis of ax1 before matching to
                # ax0
                ax1_current_model.rotsymm_ax.get_rotsymm_axis()


                # add ax1 rep to corresponding ax0 rep
                # use correct submodels of each rep
                # setting connections is done by ax_a_model_add_ax_b()
                ax1_submodel_id = self.ax1_contact_submodel_0( \
                                    ax0_submodel_id, ax1.fold)

                ax_0_model_id = (1, ax0_submodel_id)
                ax_1_model_id =  (ax1_current_model.id[0], ax1_submodel_id)

                self.ax_a_model_add_ax_b(ax0, ax1, \
                                    ax_0_model_id, ax_1_model_id)


                # get pivot residue for alignment of 2 SymPlexes
                self.alignment_pivot_res(ax_0_model_id, ax_1_model_id, sess)
                pivot_res = self.axes[1].alignment_pivot_res[ \
                                                    self.alignment_pivot_pos]


                model_id_tmp = sess.open_model(ax1.model_active_path)
                sess.split_model(model_id_tmp)

                model_tmp_center = structure.complex_.get_center(
                        model_id_tmp, \
                        ax1.fold, \
                        ax1_current_model.termini, \
                        sess)
                diff = ax1_current_model.get_center()-model_tmp_center

                current_model_alignment_fixppoint_xyz = \
                    sess.get_xyz(ax_1_model_id, pivot_res)

                model_tmp_alignment_fixppoint_xyz = \
                    sess.get_xyz((model_id_tmp[0], ax1_submodel_id), pivot_res)

                diff_fixppoint = current_model_alignment_fixppoint_xyz- \
                                 model_tmp_alignment_fixppoint_xyz

                sess.move_model(model_id_tmp, \
                                [diff[0], diff[1], diff_fixppoint[2]])
                        # diff_fixppoint[2]: align z regarding pivot_res


                # determine rotational symmetry axis of ax1 SymPlex
                rotsymm_axis = ax1_current_model.rotsymm_ax. \
                            get_rotsymm_axis(mode_readout_angle_change=True)


                # flip model_tmp from bottom to top if ax1 SymPlex got flipped
                # during matching to ax0 SymPlex
                if rotsymm_axis[2] < 0 or ax1.model_active_orient == -1:
                            # the z orientation of the rotational symmetry axis
                            # is used as indicator for a bottom/top flip
                    ax1_current_model.flipped = True
                    ax1.model_active_orient = -1

                    # flip from bottom to top
                    sess.turn_model(model_id_tmp, 0, 180, model_tmp_center)

                    model_tmp_alignment_fixppoint_xyz_ = \
                        sess.get_xyz((model_id_tmp[0], ax1_submodel_id), \
                                     pivot_res)

                    diff_fixppoint_ = current_model_alignment_fixppoint_xyz- \
                                     model_tmp_alignment_fixppoint_xyz_

                    # correct y shift inverted by 180 degree turn around x,
                    # correct z shift after 180 degree turn around x
                    sess.move_model(model_id_tmp, \
                                [0, 2*diff[1], diff_fixppoint_[2]])


                # determine rotation angle after flip that may have occurred
                rotang, dist = sess.get_rot_angle_z( \
                            (ax1_current_model.id[0], 1), (model_id_tmp[0], 1))

                ctl.d('rotang_after')
                ctl.d(rotang)
                ctl.d(ax1_current_model.idstr)
                ctl.d(model_id_tmp)
                ax1_current_model.set_rot_angle(rotang)

                if ax1_current_model.flipped == True:
                    ax1_current_model.set_rot_angle(rotang*(-1))

                    if renamed == False:
                        intermediate_ids = []

                        for ii in range(1, ax1.fold):

                            # inverse chain ids for flipped symplex
                            # ax1_submodel_id is fixed
                            sess.change_model_id( \
                                    (ax1_current_model.id[0], \
                                        (ax1_submodel_id+ii-1)%ax1.fold+1), \
                                    (ax1_current_model.id[0], \
                                        10+(ax1.fold+ax1_submodel_id-ii-1)% \
                                        ax1.fold+1) )
                            intermediate_ids.append( \
                                10+(ax1.fold+ax1_submodel_id-ii-1)%ax1.fold+1)

                        for ii in intermediate_ids:
                            sess.change_model_id( \
                                    (ax1_current_model.id[0], ii), \
                                    (ax1_current_model.id[0], (ii-10)) )


                # determine rotational symmetry axis of ax1
                rotsymm_axis = ax1_current_model.rotsymm_ax.get_rotsymm_axis()


                # check orientation of ax1 rotational symmetry axis
                if rotsymm_axis[2] < 0:
                    ctl.e(rotsymm_axis)
                    ctl.error('Assembler: build_layer: rotsymm_axis[2] < 0')


                # check if axial tilt of ax1 (rotational symmetry axis) is
                # within range
                ax1_tilt = ax1_current_model.rotsymm_ax.get_tilt()
                ax1_tilt_deg = math.degrees(ax1_tilt)

                if abs(ax1_tilt_deg) > 45:
                    ctl.p(rotsymm_axis)
                    ctl.p(ax1_tilt_deg)
                    ctl.e('Assembler: build_layer: '+ \
                              'tilt of ax1 not within range')
                    raise Exception('build_layer: '+ \
                              'tilt of ax1 not within range')


                # bring ax1 representation/representative in flattened position
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if max(self.conf.flatten_modes) >= 0:
                    model_tmp_center = structure.complex_.get_center(
                        model_id_tmp, \
                        ax1.fold, \
                        ax1_current_model.termini, \
                        sess)
                    
                    model_tmp_rotang = ax1_current_model.rot_angle

                    sess.turn_model(model_id_tmp, 2, model_tmp_rotang, \
                                    model_tmp_center)

                    model_tmp_domain_center = structure.monomer.get_center( \
                            (ax1_current_model.id[0], \
                                    self.ax1_contact_submodel_0( \
                                    ax0_submodel_id, ax1.fold)), \
                            (ax1.domains[0][0], ax1.domains[0][1]), \
                            ax1.chimerax_session)

                    ax0_submodel_domain_center = structure.monomer. \
                        get_center( \
                            (1, ax0_submodel_id), \
                            (ax0.domains[0][0], ax0.domains[0][1]), \
                            ax1.chimerax_session)

                    ax1.chimerax_session.match((ax1_current_model.id[0], 1), \
                                    (ax1.domains[0][0], ax1.domains[0][1]), \
                                    (model_id_tmp[0], 1), \
                                    ax1_current_model.id)

                sess.close_id(model_id_tmp)


            # step 3: place ax0 representations/representatives around ax1
            # ------------------------------------------------------------
            if max(self.conf.flatten_modes) >= 0:
                ax1_models = layer.ax_models(ax1)

                for ax1_model_i,ax1_model in enumerate(ax1_models):
                    ax0_current_model = ax0.open_model()
                    layer.add_representation(ax0_current_model)

                    source_submodel_id = self.ax1_contact_submodel_1( \
                                        ax1_model_i+1, \
                                        self.contact_submodel_orientation, \
                                        ax1.fold)
                    dest_submodel_id = (((ax1_model_i+1)+ \
                        ax0_pairing_offset+orientation_offset-1)%ax0.fold+1)

                    # add ax0 rep to corresponding ax1 rep
                    # use correct submodels of each rep
                    # setting connections is done by ax_a_model_add_ax_b()
                    self.ax_a_model_add_ax_b(
                            ax1, ax0, \
                            (ax1_model.id[0], source_submodel_id), \
                            (ax0_current_model.id[0], dest_submodel_id)
                            )

                    # set name to submodels of ax0 representation/
                    # representative
                    for i in range(1,1+ax0.fold):
                        sess.rename_model((ax0_current_model.id[0], i), \
                                'ax'+str(2+ax1_model_i)+'mol'+str(i))

                    # determine rotation of ax0 representation/representative
                    transl_vec = sess.measure_transl_vec( \
                                        (ax0_current_model.id[0], 1), (1, 1))
                    ax0_current_model.set_trans_vect(transl_vec)


            # step 4: place ax0 models flattened around ax1
            # ---------------------------------------------
            if max(self.conf.flatten_modes) >= 0:
                ax0_models = layer.ax_models(ax0, order=1)
                ax1_models = layer.ax_models(ax1)

                for ax0_model_i,ax0_model in enumerate(ax0_models):
                    ax0_current_model = ax0.open_model()
                    ax0_current_model.set_trans_vect(ax0_model.trans_vect)

                    layer_flat.add_representation(ax0_current_model)


                    # perform rotation if necessary
                    if ax0.fold == 3 and ax1.fold == 2:
                        if ax0.get_representation(1).flipped == True:
                            sign = -1
                        else:
                            sign = 1

                        sess.turn_model(ax0_current_model.id, \
                                        2, -60*sign, \
                                        ax0_current_model.get_center())


                    # set name to submodels of ax0 representation/
                    # representative
                    for i in range(1, 1+ax0.fold):
                        sess.rename_model( \
                            (ax0_current_model.id[0], i), \
                            'ax'+str(2+ax0.fold+ax1_model_i)+'mol'+str(i))

                    ax0_current_model.move_model( \
                        [ax0_model.trans_vect[0], ax0_model.trans_vect[1], 0])

                    # set connections, use correct submodels for each
                    # representation/representative
                    source_submodel_id = ((ax0_model_i+1)+ \
                            ax0_pairing_offset+orientation_offset-1)%ax0.fold+1

                    dest_submodel_id = self.ax1_contact_submodel_1( \
                            ax0_model_i+1, \
                            self.contact_submodel_orientation, \
                            ax1.fold)

                    ax0_current_model.set_connection(
                        (ax0_current_model.id[0], source_submodel_id), \
                        (ax1_models[ax0_model_i].id[0], dest_submodel_id), \
                        'flattened')

            if os.path.exists(self.conf.path_ax_predictions+ \
                            self.conf.species+'/'+self.conf.symplex_path+ \
                            'setting_meta.txt'):
                meta_path = self.conf.get_struct_coll_meta_path()
                export_file_prefix = self.conf.export_file_prefix+ \
                                     self.export_file_infix
                self.export_meta(ax0, meta_path+export_file_prefix+'.txt')

        sess.format_session()
        sess.save_models('all', self.conf.layer_raw_path)

        return


    def process_layer(self, combine_mode=2):
        '''
        Process the layer for different variations of lattice constants.
        '''
        clashes = []

        # variation of lattice constant
        for i,lc_offset in enumerate(self.lc_var):
            self.lc_offset = lc_offset
            if i == 0:
                preserve_connections = False
            else:
                preserve_connections = True

            clashes_ = self.process_layer_lc( \
                            combine_mode, \
                            preserve_connections=preserve_connections)
            clashes.append(clashes_)

            if i > 1 and abs(clashes_-clashes[1]) <= \
                        abs(clashes[0]-clashes[1])*0.1353:
                ctl.d('variation finished')
                break

        return clashes


    def process_layer_lc(self, combine_mode=2, preserve_connections=False):
        '''
        Process the layer for a given lattice constant,
        eg. create cross connections.

        create special setups:
        - snapin_layer (if needed)
        - complete_outer_chains (if needed)

        create cross connections

        export different configurations
        '''
        validation_result_tile = 1000 # return variable
        
        ax0 = self.axes[0]

        if len(self.axes) > 1:
            ax1 = self.axes[1]

        ax0.chimerax_session.init()
        
        self.layers[1].calc_refpoint_constant(self.lc_offset)

        # align layer_flat to ref points and move all other models relative to
        # it
        if len(self.axes) > 1 and max(self.conf.flatten_modes) >= 1:
            ax0.chimerax_session.open_session(self.conf.layer_raw_path)

            align_layer.align_layer(self.axes, \
                                    self.layers[1], \
                                    self.conf, \
                                    ax0.chimerax_session)
        aligned_file = sorted(glob.glob(self.conf.layer_aligned_raw_path))
        

        # snapin ax0s to ref points
        snapin_file = sorted(glob.glob(self.conf.layer_snapin_raw_path))
        if max(self.conf.flatten_modes) >= 2:
            ax0.chimerax_session.open_session(self.conf.layer_aligned_raw_path)
            snapin.snapin_layer(self.axes, self.layers[1], self.conf, \
                                preserve_connections)
            snapin_file = sorted(glob.glob(self.conf.layer_snapin_raw_path))


        meta_path = self.conf.get_struct_coll_meta_path()
        export_file_prefix = self.conf.export_file_prefix+ \
                             self.export_file_infix

        # insert lc_offset infix only when offset different to 0
        if self.lc_offset != 0:
            lc_offset_infix = '_lc'+('00'+str(self.lc_offset))[-2:]
        else:     
            lc_offset_infix = ''

        self.layers[1].export_param(self.conf.export_path+ \
                                    self.conf.export_file_prefix+ \
                                    self.export_file_infix+ \
                                    lc_offset_infix, \
                                    lc_offset=self.lc_offset)

        if os.path.exists(self.conf.path_ax_predictions+ \
            self.conf.species+'/'+self.conf.symplex_path+'setting_meta.txt'):
            self.export_meta(ax0, \
                         meta_path+export_file_prefix+lc_offset_infix+'.txt')


        model_reg_resetpoint = ax0.model_reg

        # flatten_modes 0,1,2,3,4:
        # 0: pure superposition, 
        # 1: flattened, 
        # 2: snapin/tile,
        # 3: completed chains,
        # 4: primitive unit cell,
        # 5: assembly of 3x3 primitive unit cells

        do_snapshot_w_separated_chains = 1

        for deletetermini in self.conf.delete_termini_modes:
            for flatten_mode in self.conf.flatten_modes:
                for filter_for_export in self.conf.filters_for_export:
                    for snapshot in self.conf.snapshot_modes:
                        if snapshot == 1 and \
                                do_snapshot_w_separated_chains == 0:
                            continue
                        if snapshot == 1 and \
                                do_snapshot_w_separated_chains == 1:
                            do_snapshot_w_separated_chains = 0

                        flatten = 0
                        if flatten_mode >= 1:
                            # flatten 0: superposition, 1: others
                            flatten = 1

                        # process flatten_mode 0,1 only once because here no
                        # variation of lattice constant
                        if flatten_mode < 2 and self.lc_offset > 0:
                            continue

                        # process layer for specific setup
                        validation_result_tile_ = \
                                self.process_layer_specific(deletetermini, \
                                        flatten_mode, flatten, \
                                        filter_for_export, \
                                        snapshot, \
                                        do_snapshot_w_separated_chains, \
                                        aligned_file, snapin_file, \
                                        combine_mode, model_reg_resetpoint)

                        if validation_result_tile_ != -1:
                            validation_result_tile = validation_result_tile_

        return validation_result_tile


    def process_layer_specific(self, \
                               deletetermini, flatten_mode, flatten, \
                               filter_for_export, \
                               snapshot, do_snapshot_w_separated_chains, \
                               aligned_file, snapin_file, \
                               combine_mode, model_reg_resetpoint):
        '''
        Process layer for specific setup.
        '''
        merge_likely_related_fragments = True
            # search for protein fragments that are likely part of the
            # same molecule and merge them

        validation_result_tile = -1

        if flatten > len(self.layers):
            ctl.error('Assembler: process_layer_specific: '+ \
                        'flatten > len(layers)')
        layer = self.layers[flatten]


        # for snapin/tile, completed chains, primitive unit cell
        if flatten_mode >= 2 and self.lc_offset != 0:
            lc_offset_infix = 'lc'+('00'+str(self.lc_offset))[-2:]+'_'
        else:     
            lc_offset_infix = ''

        ax0 = self.axes[0]
        if len(self.axes) > 1:
            ax1 = self.axes[1]

        ax0.chimerax_session.init()

        if len(aligned_file) == 1:
            ax0.chimerax_session.open_session( \
                    self.conf.layer_aligned_raw_path)
        else:
            ax0.chimerax_session.open_session( \
                    self.conf.layer_raw_path)

        if flatten_mode >= 2:
            if len(snapin_file) != 1:
                ctl.error('Assembler: process_layer: len(snapin_file) != 1')

            ax0.chimerax_session.open_session( \
                    self.conf.layer_snapin_raw_path)


        # reset connected flag ax0
        for r in ax0.get_representations():
              ax0.get_representation(r).connected = []

        ax0.model_reg = model_reg_resetpoint


        # reset connected flag ax1, but only if ax1 exists
        if len(self.axes) > 1:
            for r in ax1.get_representations():
                ax1.get_representation(r).connected = []

            ax1.model_reg = model_reg_resetpoint


        if deletetermini == 1:

            # raise exception if ax0 surface completely within termini range
            for i in ax0.get_representations():
                ax0_termini = ax0.get_representation(i).termini

                if ax0.surface != None:
                    if ax0.surface[0][1] <= ax0_termini[0]:

                        raise Exception('process_layer_specific: '+ \
                                    'ax surface completely within termini')

                    if ax0.surface[0][0] >= ax0_termini[1]:

                        raise Exception('process_layer_specific: '+ \
                                    'ax surface completely within termini')

            ax0.delete_termini()

            if len(self.axes) > 1:

                # raise exception if ax1 surface completely within
                # termini range
                for i in ax1.get_representations():
                    ax1_termini = ax1.get_representation(i).termini

                    if ax1.surface != None:
                        if ax1.surface[0][1] <= ax1_termini[0]:

                            raise Exception('process_layer_specific: '+ \
                                        'ax surface completely within termini')

                        if ax1.surface[0][0] >= ax1_termini[1]:

                            raise Exception('process_layer_specific: '+ \
                                        'ax surface completely within termini')

                ax1.delete_termini()


        # combine_mode 1
        if snapshot == 0 and combine_mode == 1:
            for r in ax0.get_representations():
                ax0.get_representation(r).combine_chains_all()

        # combine_mode 2
        if snapshot == 0 and combine_mode == 2:

            # determine residue ranges that will be included
            # in the overall model
            ax0.determine_surface()
            if len(self.axes) > 1:
                ax1.determine_surface()

            # modes: 'default', 'flattened', 'snapin'
            modes_connect = ['default']

            if flatten_mode == 0:
                modes_connect = ['default']
            elif flatten_mode == 1:
                modes_connect = ['flattened', 'default']
            elif flatten_mode >= 2:
                modes_connect = ['snapin', 'flattened', 'default']

            if flatten_mode == 0:
                ax0_models = self.layers[0].ax_models(ax0)
                        # ax0 models assigned to layer
            elif flatten_mode >= 1:
                ax0_models = self.layers[1].ax_models(ax0)
                        # ax0 models assigned to layer_flat


            for mode_connect in modes_connect:
                for r in ax0_models:
                    ctl.d('mode_connect:')
                    ctl.d(mode_connect)
                    r.combine_chains_all(mode_connect)


            # search for protein fragments that can be merged and merge them
            if merge_likely_related_fragments:
                self.merge_protein_fragments(flatten, modes_connect[0])


        if flatten_mode <= 2:
            # ensure that passive residues are not deleted in complete chains

            if snapshot == 0:
                # delete in unmerged models the residue ranges that will not be
                # included in the overall model
                ax0.delete_passive_of_unmerged()
                if len(self.axes) > 1:
                    ax1.delete_passive_of_unmerged()

            if snapshot == 1: # snapshot mode
                ax0.delete_passive()
                if len(self.axes) > 1:
                    ax1.delete_passive()

        # complete chains
        if flatten_mode == 3: 
            complete_outer_chains.complete_outer_chains( \
                    self, layer, self.contact_submodel_orientation)

            # search for protein fragments that can be merged and merge them
            if merge_likely_related_fragments:
                self.merge_protein_fragments(flatten, modes_connect[0])


        # model id for combine
        combination_model_id = ax0.chimerax_session.last_id()+1

        # close combination_model_id in case it is empty group
        ax0.chimerax_session.close_id(combination_model_id)


        export_path = \
            [[self.conf.export_path_pure_superposition, \
              self.conf.export_path, \
              self.conf.export_path_tile, \
              self.conf.export_path_complete_chains, \
              self.conf.export_path_primitive_unit_cell, \
              self.conf.export_path_assembly] \
              [flatten_mode], \
            self.conf.export_path_snapshot][snapshot]+ \
            self.conf.export_path_filters[filter_for_export]+ \
            self.export_file_infix+ \
            self.conf.fn_infix_termini[deletetermini]+'_'+ \
            ['', \
            'flat_', \
            'tile_', \
            'complete_chains_', \
            'primitive_unit_cell_', \
            'assembly_3x3_'][flatten_mode]+ \
            lc_offset_infix+'v'+self.conf.version+'.cif'
        export_path = filesystem.clean_path(export_path)

        # superposition
        if flatten_mode < 4:
            if len(self.axes) > 1:

                if filter_for_export == 0:
                    if flatten == 0:
                        layer_models = self.layers[0].get_representations()

                    elif flatten == 1:
                        layer_models = self.layers[1].get_representations()

                    ax0.chimerax_session.combine_models( \
                            layer_models.values(), combination_model_id)
    
                # axis 1
                if filter_for_export == 1:
                    layer_ax1_models = self.layers[0].ax_models(ax1)
                    ax0.chimerax_session.combine_models( \
                            layer_ax1_models, combination_model_id)

                # axis 0
                if filter_for_export == 2:
                    layer_ax0_models = self.layers[0].ax_models(ax0)
                    ax0.chimerax_session.combine_models( \
                            layer_ax0_models, combination_model_id)


            if len(self.axes) == 1:
                  # only axis 0
                  ax0.chimerax_session. \
                      run('combine #1 close false modelId #'+ \
                          str(combination_model_id))

            # validate export combined model
            validation_results = validation.validation( \
                                    (combination_model_id,), \
                                    (combination_model_id,), \
                                    self.axes, \
                                    export_path, \
                                    self.conf, \
                                    ax0.chimerax_session)

            # export
            export.compatibility_cif_export(export_path, \
                    combination_model_id, \
                    ax0.chimerax_session, \
                    self.conf.cif_postprocess)

            # 2: snapin/tile
            if flatten_mode == 2:
                validation_result_tile = validation_results[0]


        # flatten_mode 4: primitive unit cell
        elif flatten_mode == 4:

            # create and export primitive unit cell if possible
            self.layers[1].primitive_unit_cell.definition(self.lc_offset)
            self.layers[1].primitive_unit_cell.suggest_molecules()
            self.layers[1].primitive_unit_cell. \
                                        combine_to_model(combination_model_id)

            self.layers[1].primitive_unit_cell. \
                    export(combination_model_id, export_path, self.conf)            


        # flatten_mode 5: assembly 3x3
        elif flatten_mode == 5:

            # create and export primitive unit cell if possible
            self.layers[1].primitive_unit_cell.definition(self.lc_offset)
            monomers_complete = \
                        self.layers[1].primitive_unit_cell.suggest_molecules()

            # raise exception if monomer(s) of primitive unit cell not complete
            if monomers_complete == False:
                raise Exception('process_layer_specific: '+ \
                                        'monomers not complete')

            model_complete = self.layers[1].primitive_unit_cell. \
                                        combine_to_model(combination_model_id)

            if model_complete != False:
                mates_n = 3
                a = self.layers[1].primitive_unit_cell.a
                b = self.layers[1].primitive_unit_cell.b

                x = a
                y = 0
                ax0.chimerax_session. \
                    run('sym #'+str(combination_model_id)+ \
                        ' shift,'+str(mates_n)+','+str(x)+','+str(y)+',0')

                current_model_id = ax0.chimerax_session.last_id()

                x = a*math.cos(math.radians( \
                                    self.layers[1].primitive_unit_cell.gamma))
                y = b*math.sin(math.radians( \
                                    self.layers[1].primitive_unit_cell.gamma))
                ax0.chimerax_session. \
                    run('sym #'+str(current_model_id)+ \
                        ' shift,'+str(mates_n)+','+str(x)+','+str(y)+',0')

                current_model_id = ax0.chimerax_session.last_id()

                # model id of the center primitive unit cell (puc) in the
                # 3x3 assembly
                central_puc_id = (current_model_id, 2, 2)

                # validate model
                validation_scores = validation.validation(central_puc_id, \
                                                    (current_model_id,), \
                                                    self.axes, \
                                                    export_path, \
                                                    self.conf, \
                                                    ax0.chimerax_session)

                self.layers[1].primitive_unit_cell.validation_scores = \
                                                    validation_scores[1]

                export.compatibility_cif_export(export_path, \
                        current_model_id, \
                        ax0.chimerax_session, \
                        self.conf.cif_postprocess)

        return validation_result_tile


    def merge_protein_fragments(self, flatten, mode='default'):
        '''
        Global search for protein fragments that can be merged and merge
        them.

        modes: 'default', 'flattened', 'snapin'
        '''
        all_joins = []

        if len(self.axes) < 2:
            return all_joins

        ax0 = self.axes[0]
        ax1 = self.axes[1]
        ax0_models = ax0.get_representations()
        ax1_models = ax1.get_representations()

        layer_models = self.layers[0].get_representations()
        layer_flat_models = self.layers[1].get_representations()

        ax0_model_ids = [ax0_models[m].id[0] for m in ax0_models]
        ax1_model_ids = [ax1_models[m].id[0] for m in ax1_models]
        layer_model_ids = [layer_models[m].id[0] for m in layer_models]
        layer_flat_model_ids = [layer_flat_models[m].id[0] \
                                            for m in layer_flat_models]

        # model ids not assigned to a layer (e.g. models used for
        # completion of outer chains)
        ax_not_layer_ids = sorted(list( \
                (set(ax0_model_ids) | set(ax1_model_ids)) - \
                (set(layer_model_ids) | set(layer_flat_model_ids)) ))


        if flatten == 0:

            # axis 0
            layer_ax0_models = self.layers[0].ax_models(ax0)
            ids1 = [m.id[0] for m in layer_ax0_models]

            # others
            layer_ax1_models = self.layers[0].ax_models(ax1)
            ids2 = [m.id[0] for m in layer_ax1_models]+ax_not_layer_ids


        ids1_part2 = []
        if flatten == 1:

            # axis 0
            layer_flat_ax0_models = self.layers[1].ax_models(ax0)
            ids1 = [m.id[0] for m in layer_flat_ax0_models]

            # others
            layer_flat_ax1_models = self.layers[1].ax_models(ax1)
            ids2 = [m.id[0] for m in layer_flat_ax1_models]+ax_not_layer_ids


            # completed chains (completion of ax1)
            layer_flat_ax1_models = self.layers[1].ax_models(ax1)
            ids1_part2 = [m.id[0] for m in layer_flat_ax1_models]
            ids2_part2 = ax_not_layer_ids
                    # model ids not assigned to a layer (e.g. models used for
                    # completion of outer chains)


        for id1 in ids1:
            if not ax0.model_reg.model_exists(id1):
                continue

            to_join = proteinmerge.complex_merge(id1, ids2, \
                                                 ax0.chimerax_session)
            all_joins += to_join

            for tj in to_join:
                if mode == 'snapin':
                    ctl.d(mode)
                    ctl.d(flatten)
                    ctl.d('Assembler: merge_protein_fragments: snapin 1')
                    ctl.d(to_join)

                c1 = ax0.model_reg.get_model(int(tj[0].split('.')[0])). \
                             get_connection(tj[0], mode)
                c2 = ax0.model_reg.get_model(int(tj[1].split('.')[0])). \
                             get_connection(tj[1], mode)                             
                  
                if c1 == False and c2 != False:
                    ctl.d(mode)
                    ctl.d(tj)
                    ctl.d(c1)
                    ctl.d(c2)
                    ctl.d('Assembler: merge_protein_fragments: '+ \
                              'unequal get_connections, 1')

                if c2 == False and c1 != False:
                    # can be error or just large search radius/dist cutoff
                    # that leads to connections that do not fit
                    ctl.d(mode)
                    ctl.d(tj)
                    ctl.d(c1)
                    ctl.d(c2)
                    ctl.d('Assembler: merge_protein_fragments: '+ \
                          'unequal get_connections, 2')
                    continue

                if c1 != False and c2 != False:
                    ctl.d(mode)
                    ctl.d(tj)
                    ctl.d(c1)
                    ctl.d(c2)
                    ctl.d('merge_protein_fragments: '+ \
                          'both already existing, 5')
                    continue

                combined = bib.combine_chains(tj[0], tj[1], \
                                              ax0.chimerax_session)


        for id1_part2 in ids1_part2:
            if not ax0.model_reg.model_exists(id1_part2):
                continue

            to_join = proteinmerge.complex_merge(id1_part2, ids2_part2, \
                                                 ax0.chimerax_session)
            all_joins += to_join

            for tj in to_join:
                if mode == 'snapin':
                    ctl.d(mode)
                    ctl.d(flatten)
                    ctl.d('Assembler: merge_protein_fragments: snapin 2')
                    ctl.d(to_join)

                c1 = ax0.model_reg.get_model(int(tj[0].split('.')[0])). \
                            get_connection(tj[0], mode)
                c2 = ax0.model_reg.get_model(int(tj[1].split('.')[0])). \
                            get_connection(tj[1], mode)

                if c1 == False and c2 != False:
                    ctl.e(tj)
                    ctl.e(c1)
                    ctl.e(c2)
                    ctl.error('Assembler: merge_protein_fragments: '+ \
                              'unequal get_connections, 3')

                if c2 == False and c1 != False:
                    ctl.e(tj)
                    ctl.e(c1)
                    ctl.e(c2)
                    ctl.error('Assembler: merge_protein_fragments: '+ \
                              'unequal get_connections, 4')

                if c1 != False and c2 != False:
                    ctl.d(tj)
                    ctl.d(c1)
                    ctl.d(c2)
                    ctl.d('Assembler: merge_protein_fragments: '+ \
                          'both already existing, 6')
                    continue

                combined = bib.combine_chains(tj[0], tj[1], \
                                              ax0.chimerax_session)

        return all_joins


    def alignment_pivot_res(self, ax0_model_id, ax1_model_id, sess):
        '''
        Determine list of 2 pivot points for alignment of 2 SymPlexes.

        Return:
            pivot_res: [pivot residue at the upstream side of the
                        alignment section (e.g. domain),
                        pivot residue at the downstream side of the
                        alignment section (e.g. domain)]
        '''
        res_ov = sess.get_residue_overlap(ax0_model_id, ax1_model_id)

        coin_res, res_d = self.get_coincident_residues( \
                ax0_model_id, ax1_model_id, res_ov, sess)

        if len(coin_res) == 0:
            coin_res, res_d = self.get_coincident_residues( \
                    ax0_model_id, ax1_model_id, res_ov, sess, 2)

        if len(coin_res) == 0:
            raise Exception('Assembler: alignment_pivot_res: '+ \
                    'no coincident residues')

        if coin_res[0] < coin_res[-1]:
            pivot_res = [coin_res[0], coin_res[-1]]
        else:
            pivot_res = [coin_res[-1], coin_res[0]]

        self.axes[0].alignment_pivot_res = pivot_res
        self.axes[1].alignment_pivot_res = pivot_res

        return pivot_res


    def export_meta(self, ax, file):
        ''' Export meta data. '''

        f = open(file, 'w')

        json_export = []
        for ax_rep in ax.get_representations():
            representation = ax.get_representation(ax_rep)

            json_export.append([
                            representation.termini, \
                            representation.multimer_n, \
                            representation.domains, \
                            representation.surface, \
                            representation.trans_vect \
                            ])

        if self.axes[0].alignment_pivot_res != []:
            json_export.append([self.axes[0].alignment_pivot_res, \
                            self.alignment_pivot_pos, \
                            self.axes[0].alignment_pivot_res \
                            [self.alignment_pivot_pos]])

        json.dump(json_export, f)
        f.close()

        return
