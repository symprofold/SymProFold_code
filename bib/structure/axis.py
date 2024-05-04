import bib
import bibpdb
import ctl
import filesystem
import geometry
import molmodel
from structure.model import Model
from structure.monomer import Monomer
from structure.referenceframe import ReferenceFrame
from structure.rotsymm_axis import RotSymmAxis

import glob
import math
import os
import shutil


class Axis():
    '''
    This class describes in abstract form one of the rotational symmetry
    complexes (SymPlexes), out of which the layer is assembled.
    The concept includes the rotational symmetry axis and the partial molecular
    models which are subject to this symmetry axis.
    E.g. the abstract form of a sixfold axis.
    '''

        
    def __init__(self, model_reg):
        ''' Initialization of the Axis class. '''

        self.model_reg = model_reg
        self.fold = 0
        self.model_active_path = ''
        self.model_active_orient = 1
        self.export_path = ''
        self.meta = {}
        self.representations = {}
        self.representations_order = [{}, {}]
        self.conf = 0

        self.domains = [ [1,10000] ]
                # 10000 to cover all common protein lengths
        self.surface = None
        self.surface_resids = []
        self.passive_resids = []

        self.chimerax_session = None
        self.preferred_bs = 0
        self.alignment_pivot_res = []

        # orientation of chainids (counterclockwise, clockwise)
        self.chainid_orient = 1 # 1: counterclockwise

        return


    def set_session(self, session, conf):
        ''' Set Chimerax session. '''
        self.chimerax_session = session # class ChimeraxSession
        self.conf = conf
        return


    def set_fold(self, fold):
        ''' Set multiplicity of the rotation symmetry axis. '''
        self.fold = fold
        return


    def set_folder(self, folder, orient):
        '''
        Set folder of the input models from which the overall model is
        composed and index the input models.
        '''
        self.pathRaw = self.conf.path_ax_predictions+ \
                       self.conf.species+'/'+folder

        if orient == 1:
            self.path = self.pathRaw[:-1]+'_o/'
        elif orient == 2:
            self.path = self.pathRaw[:-1]+'_or/'
        else:
            self.path = self.pathRaw

        # search for models in given folder according to priority order
        files = sorted(glob.glob(self.path+'rank*_*.*.pdb'))
        if len(files) == 0:
            files = sorted(glob.glob(
                self.path+'unrelaxed_rank*_*.*.pdb'))
        if len(files) == 0:
            files = sorted(glob.glob(
                self.path+'relaxed_model_*_multimer_*_pred_*.pdb'))
        if len(files) == 0:
            files = sorted(glob.glob(
                self.path+'unrelaxed_model_*_multimer_*_pred_*.pdb'))
        if len(files) == 0:
            ctl.e(self.path)
            ctl.error('set_folder: no models found in folder.')


        filenames = []
        for f in files:
            filenames.append(f.split('/')[-1])

        self.files = filenames

        # set first model active
        if self.model_active_path == '':
            self.set_model_active(0)

        return


    def add_representation(self, representation):
        ''' Add representation. '''
        self.representations[representation.id] = representation
        return


    def get_representation(self, model_id):
        ''' Get representation. '''

        model_id = self.model_reg.convert_model_id(model_id)
        representation = self.representations[model_id]

        return representation


    def get_representations(self):
        ''' Get representations. '''

        representations = {}
        for r in self.representations:
            representations[self.model_reg.convert_model_id(r)] = \
                                                self.representations[r]

        return representations


    def set_path(self, path, orient):
        '''
        Set path of the input models directly without indexing the input
        models.
        '''        
        self.pathRaw = path
        self.path = path[:-1]+'_o/' if orient == 1 else path
        return  


    def set_files(self, files):
        ''' Set index of the input models directly. '''
        self.files = files
        return


    def set_domains(self, domains):
        ''' Set domain ranges. '''
        self.domains = domains
        return


    def set_surface(self, surface):
        ''' Set residue range that will be included in the overall model. '''

        self.surface = surface
        _surface_resids = []

        for i in range(0,2000): # 2000: max value for resids
            for s in surface:
                if s[0] <= i <= s[1]:
                    _surface_resids.append(i)
    
        self.surface_resids = list(set(_surface_resids))

        # set/update new values in all representatants
        for ax_rep in self.get_representations():
            self.get_representation(ax_rep).set_surface(self.surface)

        return  


    def determine_surface(self):
        '''
        Determine residue range that will be included in the overall model.
        '''
        if self.surface == []:
            self.set_surface([[1, 2000]])

        return  


    def set_model_active(self, model_active, orient=False):
        ''' Flag the active model. '''
        self.model_active = model_active
        self.model_active_path = self.path+self.files[self.model_active]

        # orientation of chainids (counterclockwise, clockwise)
        if orient == False: # counterclockwise
            self.chainid_orient = 1
        if orient == True: # clockwise
            self.chainid_orient = -1

        self.model_active_orient = 1
        self.fold = bib.get_multimer_n(self.model_active_path)

        return


    def set_preferred_bs(self, preferred_bs):
        ''' Set the preferred binding site. '''
        self.preferred_bs = preferred_bs
        return


    def delete_termini(self):
        '''
        Delete termini (N-termini and C-termini) of all representations.
        '''
        representations = self.get_representations()

        for i in representations:
            representations[i].delete_termini()

        return


    def has_representatative(self, rep_id):
        ''' Check if a apecific representations exists. '''

        rep_id = self.model_reg.convert_model_id(rep_id)
        try:
            r = self.get_representations()[rep_id]

            return True

        except KeyError:

            return False


    def delete_passive(self):
        '''
        Delete residues ranges that will not be included in the overall model.
        '''
        if self.surface != None:
            for i, sur in enumerate(self.surface):
                for r in self.get_representations():
                    if i == 0:
                        self.chimerax_session.delete_res_range( \
                            r, [-1, self.surface[0][0]-1])
                        self.chimerax_session.delete_res_range( \
                            r, [self.surface[-1][1]+1, 10000])

                    if i < len(self.surface)-1:
                        self.chimerax_session.delete_res_range( \
                            r, [self.surface[i][1]+1, self.surface[i+1][0]-1])

        return


    def delete_passive_of_unmerged(self):
        '''
        Delete in unmerged models the residue ranges that will not be
        included in the overall model.
        '''
        for r in self.get_representations():
            subids = self.chimerax_session.get_submodel_ids(r)
         
            for submodel_id in subids:
                subid = self.model_reg.convert_model_id_to_str(submodel_id)
                
                if self.model_reg.get_model(submodel_id). \
                                                modelling_completeness != 2:
                    if self.surface != None:
                        for i, sur in enumerate(self.surface):
                            if i == 0:
                                self.chimerax_session.delete_res_range( \
                                    submodel_id, \
                                    [-1, self.surface[0][0]-1])
                                self.chimerax_session.delete_res_range( \
                                    submodel_id, \
                                    [self.surface[-1][1]+1, 10000])

                            if i < len(self.surface)-1:
                                self.chimerax_session.delete_res_range( \
                                    submodel_id, \
                                    [self.surface[i][1]+1, \
                                     self.surface[i+1][0]-1])

        return


    def export_meta(self, file):
        ''' Not used yet because conformation label not implemented yet. '''

        file = self.conf.get_struct_coll_meta_path()+ \
               conf.export_file_prefix+'.txt'
        f = open(file, 'w')
       
        for ax_rep in self.get_representations():
            f.write(''+str([self.get_representation(ax_rep).termini, \
                            self.get_representation(ax_rep).multimer_n, \
                            self.get_representation(ax_rep).domains, \
                            self.get_representation(ax_rep).surface, \
                            self.get_representation(ax_rep).trans_vect])+ \
                            "\r\n")
        f.close()
        
        return


    def open_model(self, register_writeprotection=True, split=True):
        '''
        Open axis model and create new representation of axis (class Axis_rep).
        '''

        # check if file exists
        if not os.path.exists(self.model_active_path):
            ctl.e(self.model_active_path)
            ctl.error('Axis: open_model: file does not exist.')

        opened_model_id = self.chimerax_session.open_model( \
                                                    self.model_active_path)

        self.chimerax_session.reset_center_of_rotation()    


        if self.conf.export_ax_predictions == True:
            filesystem.create_folder([self.conf.export_path_ax_predictions])

            # path to original predicted model file
            if '_or' in self.model_active_path:
                file_to_export = self.model_active_path. \
                        replace('_or', ''). \
                        replace('_o', ''). \
                        replace('symm_180/', ''). \
                        replace('symm_120/', ''). \
                        replace('symm_090/', ''). \
                        replace('symm_060/', ''). \
                        replace('symm_045/', ''). \
                        replace('_x6_superposed', '_x3')

                folder = file_to_export.split('/')[-2]
                file_to_export = file_to_export.replace( \
                                    folder+'/', folder+'/'+folder+'/')

            else:
                file_to_export = self.model_active_path. \
                        replace('_o', ''). \
                        replace('symm_180/', ''). \
                        replace('symm_120/', ''). \
                        replace('symm_090/', ''). \
                        replace('symm_060/', ''). \
                        replace('symm_045/', ''). \
                        replace('_x6_superposed', '_x3')

            dest_filename = file_to_export.split('/')[-2]+'_'+ \
                                file_to_export.split('/')[-1]


            # check if file already exists to export this model only once
            files = sorted(glob.glob(
                    self.conf.export_path_ax_predictions+dest_filename))

            if len(files) == 0:
                shutil.copy(file_to_export, \
                            self.conf.export_path_ax_predictions+ \
                                dest_filename)
            else:
                ctl.d('file already exists')
                ctl.d(files)

              
        if self.conf.performance == False:
            self.chimerax_session.hide_atoms([opened_model_id])
            self.chimerax_session.show_cartoons(opened_model_id)
        elif self.conf.performance == True:
            self.chimerax_session.hide_atoms([opened_model_id])
            self.chimerax_session.hide_cartoons(opened_model_id)

        if split:
            self.chimerax_session.split_model(opened_model_id)


        rmsds = bibpdb.get_rmsds(self.model_active_path)
        termini = molmodel.get_termini(rmsds)
        resids = self.chimerax_session.resids((opened_model_id[0],1))

        multimer_n = bib.get_multimer_n(self.model_active_path)


        # create Axis_rep object and set it as representative of
        # respective axis
        ax_rep = Axis_rep(self.model_reg)
        ax_rep.set_id(opened_model_id)
        self.add_representation(ax_rep)

        # set abstract Axis object
        ax_rep.axis = self

        # register model in model_reg and assign Monomer objects to
        # Axis_rep object
        self.model_reg.add_model(opened_model_id, ax_rep, \
                                        'axis', register_writeprotection)

        for i in range(1, multimer_n+1):
            current_model_monomer = Monomer(self.model_reg, (opened_model_id[0], i))
            self.model_reg.add_model((opened_model_id[0], i), \
                current_model_monomer, 'monomer', register_writeprotection)
            ax_rep.monomers[(opened_model_id[0], i)] = current_model_monomer


        # set parameter of Axis_rep object
        ax_rep.set_session(self.chimerax_session)
        ax_rep.set_multimer_n(multimer_n)
        ax_rep.set_fold(self.fold)
        ax_rep.resids_initial = resids
        ax_rep.set_termini(termini)
        ax_rep.set_domains(self.domains)
        ctl.d('self.surface')
        ctl.d(self.surface)
        ax_rep.set_surface(self.surface)
        ctl.d(ax_rep.surface)
        ctl.d(ax_rep.surface_resids)
        ax_rep.set_domain_binding_status()
        ax_rep.conf = self.conf
        ax_rep.init_referenceframe()
        ax_rep.init_rotsymm_axis()


        if len(self.get_representations()) == 1:
            self.representations_order[0] = ax_rep
        else:
            self.representations_order[1][opened_model_id] = ax_rep

        return ax_rep


class Axis_rep(Model, Axis):
    '''
    This class describes one single rotational symmetry complex (SymPlex),
    which is one building block of the layer.
    The concept includes the rotational symmetry axis and the partial molecular
    models which are subject to this symmetry axis.    
    E.g. the concrete model of a sixfold axis.
    '''

    def __init__(self, model_reg):
        '''
        Initialization of the Axis_rep class.
        '''
        self.axis = 0 # underlying abstract Axis object
        self.model_reg = model_reg

        self.id = None
        self.resids_initial = []
        self.termini = [-1, -1]
        self.multimer_n = 0
        self.rot_matrix = []
        self.rot_angle = None
        self.flipped = False
        self.trans_vect = [0, 0]
        self.sequence = ''
        self.monomers = {} # ids of monomers

        self.connections = { 'default': {}, \
                             'flattened': {}, \
                             'snapin': {}, \
                             'completechains': {} }
        self.connected = []
        self.domains = []
        self.conf = 0

        self.frame = None # (local) reference frame
        self.rotsymm_ax = None # rotational symmetry axis

        return


    def init_referenceframe(self):
        '''
        Initialize (local) reference frame for the SymPlex.
        '''
        self.frame = ReferenceFrame(self)

        return


    def init_rotsymm_axis(self):
        '''
        Initialize rotational symmetry axis.
        '''
        self.rotsymm_ax = RotSymmAxis(self)

        return


    def get_center_res(self, pos=0.5):
        '''
        Get residue between both termini.
        '''
        if self.termini != [-1, -1]:
            center_res = self.termini[0]+ \
                    round( (self.termini[1]-self.termini[0])*pos )

            return center_res
        else:

            return -1


    def determine_sequence(self):
        ''' Determine and set sequence. '''
        self.sequence = self.chimerax_session.get_seq((self.id[0],1))
        return


    def set_termini(self, termini):
        ''' Set N-terminus and C-terminus. '''
        self.termini = termini

        return


    def set_domains(self, domains):
        ''' Set domain ranges. '''
        self.domains = domains
        return


    def set_multimer_n(self, multimer_n):
        ''' Set multiplicity of the model. '''
        self.multimer_n = multimer_n
        return


    def set_domain_binding_status(self):
        ''' Determine and set binding status for each submodel. '''

        subids = self.chimerax_session.get_submodel_ids(self.id)

        return


    def set_trans_vect(self, trans_vect):
        ''' Set translation vector of the model. '''
        self.trans_vect = trans_vect
        return


    def set_rot_matrix(self, rot_matrix):
        '''
        Set rotational symmetry axis of the SymPlex after superposition.
        '''
        self.rot_matrix = rot_matrix
        return    


    def set_rot_angle(self, rot_angle):
        ''' Set rotation angle of the model. '''
        self.rot_angle = rot_angle
        return    


    def set_connection(self, source, dest, mode='default'):
        '''
        Set connection between 2 submodels. 

        modes: 'default', 'flattened', 'snapin'
        '''
        source_id = self.model_reg.convert_model_id(source)
        dest_id = self.model_reg.convert_model_id(dest)

        if source_id[0] != self.id[0]:
            ctl.e(source_id[0])
            ctl.e(self.id)
            ctl.error('set_connection: starting id does not fit to model')

        # check if connection already set
        if source_id in self.connections[mode]:
            ctl.e('mode: '+mode)
            ctl.e(source_id)
            ctl.e('existing dest:')
            ctl.e(self.connections[mode][source_id])
            ctl.e('new dest:')
            ctl.e(dest_id)
            ctl.error('set_connection: '+ \
                      'there is already a connection from this source')

        else:
            self.connections[mode][source_id] = dest_id

            # set reverse connection (back from dest mol)
            reverse_model = self.model_reg.get_model(dest_id[0])

            if source_id in reverse_model.connections[mode] and \
               source_id != dest_id:
                    # source != dest: allow self-connection to overwrite
                    # inherited connections

                ctl.e('mode: '+mode)
                ctl.e(source_id)
                ctl.e('existing dest:')
                ctl.e(self.connections[mode][source_id])
                ctl.e('new dest:')
                ctl.e(dest_id)
                ctl.error('set_connection: reverse: '+ \
                          'there is already a connection from this source')

            else:
                reverse_model.connections[mode][dest_id] = source_id

        return


    def get_connection(self, source, mode='default'):
        '''
        Get connection for source in given mode, with fallback to standard if
        no connection available for given mode.

        modes: 'default', 'flattened', 'snapin'
        '''
        if mode == 'snapin':
            modes = ['snapin', 'flattened', 'default']

        if mode == 'flattened':
            modes = ['flattened', 'default']

        if mode == 'default':
            modes = ['default']

        source_id = self.model_reg.convert_model_id(source)

        for mode in modes:

            if source_id in self.connections[mode]:
                return self.connections[mode][source_id]

        return False


    def get_connections(self, mode='default'):
        '''
        Get all connections of given mode, but no fallback to standard.

        modes: 'default', 'flattened', 'snapin'
        '''
        conn = {}
        
        for c in self.connections[mode]:
            conn[c] = self.connections[mode][c]

        return conn


    def get_connections_for_modes(self, modes=['default']):
        '''
        Get all connections that are available in one of the given modes.
        If one connection is available in several modes, choose the first mode
        in the list.

        modes: 'default', 'flattened', 'snapin'
        '''
        conn = {}
        for mode in modes:

            for c in self.connections[mode]:
                if c not in conn:
                    conn[c] = self.connections[mode][c]

        return conn


    def delete_connection(self, model_id, mode='default'):
        '''
        Delete connection with given model_id as origin, but do not delete the
        reverse connection.

        modes: 'default', 'flattened', 'snapin'
        '''
        model_id = self.model_reg.convert_model_id(model_id)

        if model_id in self.connections:
            del self.connections[mode][model_id]

        else:
            ctl.error('delete_connection: '+ \
                      'connection to delete is not in connections')

        return


    def delete_connections(self):
        '''
        Delete all connections of all modes.

        modes: 'default', 'flattened', 'snapin'
        '''
        self.connections = { 'default': {}, 'flattened': {}, 'snapin': {} }

        return


    def set_surface(self, surface):
        ''' Set residue range that will be included in the overall model. '''

        self.surface = surface
        _surface_resids = []

        for i in range(0, 2000): # 2000: max value of resids
            if surface == None:
                _surface_resids.append(i)
            else:
                for s in surface:
                    if s[0] <= i <= s[1]:
                        _surface_resids.append(i)
    
        self.surface_resids = list(set(_surface_resids))

        return  


    def delete_termini(self):
        '''
        Delete N-terminus and C-terminus in relation to the full length
        sequence.
        '''
        if self.conf.domains[0][0] in self.resids_initial:
            self.chimerax_session.delete_res_range(self.id, \
                                                   [-1, self.termini[0]])

        if self.conf.domains[-1][1] in self.resids_initial:
            self.chimerax_session.delete_res_range(self.id, \
                                                   [self.termini[1], 10000])

        return


    def get_center(self):
        ''' Get center of whole axis molecule complex. '''

        center = [0,0,0]
        sum_vect = 0
        sum_n = 0        
        
        for submodel_id in range(1, self.multimer_n+1):
            resids = self.chimerax_session.resids((self.id[0], submodel_id))

            for resid in resids:
                if self.termini[0] <= resid <= self.termini[1]:
                    sum_vect = sum_vect+self.chimerax_session. \
                               get_xyz((self.id[0], submodel_id), resid)
                    sum_n += 1

        center = sum_vect/sum_n

        return center


    def get_center_cylinder():
        '''
        Get center of residues inside a cylinder with midpoint of the center
        of all residues.
        '''
        center = [0, 0, 0]
        center1 = self.get_center()

        sum_vect = 0
        sum_n = 0

        cylinder_r = 50
        cylinder_h = 200
        
        for submodel_id in range(1,self.multimer_n+1):
            resids = self.chimerax_session.resids((self.id[0], submodel_id))

            for resid in resids:
                p = sess.get_xyz((self.id[0], submodel_id), resid)
                d = geometry.dist(p, [center1[0], center1[1], p[2]])

                if d <= cylinder_r and abs(p[2]-center1[2]) <= cylinder_h/2:
                    sum_vect += sess.get_xyz((self.id[0], submodel_id), resid)
                    sum_n += 1

        center = sum_vect/sum_n

        return center


    def delete_passive(self):
        '''
        Delete residues ranges that will not be included in the overall model.
        '''
        for i, sur in enumerate(self.surface):
            if i == 0:
                self.chimerax_session.delete_res_range(self.id, \
                                [-1, self.surface[0][0]-1])
                self.chimerax_session.delete_res_range(self.id, \
                                [self.surface[-1][1]+1, 10000])

            if i < len(self.surface)-1:
                self.chimerax_session.delete_res_range(self.id, \
                                [self.surface[i][1]+1, self.surface[i+1][0]-1])

        return


    def combine_chains_all(self, mode='default'):
        '''
        Combine submodels with connections between them.

        modes: 'default', 'flattened', 'snapin'
        '''
        _connected = []

        conn = self.get_connections(mode)
        
        for c in conn:
            bib.combine_chains(c, conn[c], self.chimerax_session)

        return
