import ctl
import control_utils
import bibfasta
import filesystem

import glob
import os


class Config():
    '''
    This class represents the configuration.
    '''

    def __init__(self, execution_file):
        ''' Initialization of the Config class. '''

        self.export_ax_predictions = True
        self.version = 'v0-84'

        self.main_dir = filesystem.clean_path( \
                os.path.dirname(os.path.realpath(__file__))+'/../')
        self.execution_dir = os.path.dirname(execution_file)+'/'

        self.path_ax_predictions = filesystem.clean_path( \
                os.path.dirname(os.path.realpath(__file__))+ \
                '/../preassemblies/')

        self.struct_coll_path = ''
        self.struct_coll_meta_path = ''
        self.symplex_path = ''
        self.export_path = ''
        self.export_path_postfix = ''
        self.export_path_snapshot = ''
        self.export_file_prefix = ''
        self.export_path_filters = []
        self.meta_file_prefix = ''
        self.layer_raw_path = ''

        self.species = ''
        self.species_name = ''
        self.gene = ''
        self.domains = []

        self.performance = False

        self.fn_infix_termini = ['_FL', '']

        # get runlevel
        ret = 0
        try:
           if run_level:
              ctl.d('run_level set')
              ret = run_level
        except NameError:
          ctl.d('run_level not set')
          ret = 0

        self.runlevel = ret

        # runlevel 5, full
        self.conformations =        [0] # standard

        # flatten_modes 0,1,2,3,4,5:
        # 0: pure superposition, 
        # 1: flattened, 
        # 2: snapin/tile,
        # 3: completed chains,
        # 4: primitive unit cell,
        # 5: assembly of 3x3 primitive unit cells
        self.flatten_modes =        [0,1,2,4,5]

        self.delete_termini_modes = [0,1]
        self.filters_for_export =   [0] # standard, only one option because
                                        # mostly no multiaxis domain
        self.snapshot_modes =       [0,1]

        self.cif_postprocess = True

        return


    def set_conformations(self, conformations, runlevel):
        ''' Set conformation for given runlevel. '''
        if runlevel == self.runlevel:
            self.conformations = conformations
        return


    def set_flatten_modes(self, flatten_modes, runlevel):
        ''' Set flatten modes for given runlevel. '''
        if runlevel == self.runlevel:
            self.flatten_modes = flatten_modes
        return


    def set_delete_termini_modes(self, delete_termini_modes, runlevel):
        ''' Set "delete termini" modes for given runlevel. '''
        if runlevel == self.runlevel:
            self.delete_termini_modes = delete_termini_modes        
        return


    def set_filters_for_export(self, filters_for_export, runlevel):
        ''' Set "filters for export" modes for given runlevel. '''
        if runlevel == self.runlevel:
            self.filters_for_export = filters_for_export        
        return


    def set_snapshot_modes(self, snapshot_modes, runlevel):
        ''' Set "snapshot" modes for given runlevel. '''
        if runlevel == self.runlevel:
            self.snapshot_modes = snapshot_modes
        return


    def set_species(self, species, species_name=''):
        ''' Set name for species. '''

        self.species = species
        self.species_name = species_name

        self.export_file_prefix = ''
        self.export_path_filters = ['']

        self.update_export_paths()

        return


    def set_symplex_path(self, symplex_path):
        ''' Set import path for SymPlexes. '''

        self.symplex_path = symplex_path

        return


    def set_export_path_filters(self, export_path_filters):
        ''' Set export path filters. '''

        self.export_path_filters = []

        for p in export_path_filters:
            self.export_path_filters.append(p)

        return


    def set_export_path_postfix(self, export_path_postfix):
        ''' Set export path postfix. '''

        self.export_path_postfix = export_path_postfix
        self.update_export_paths()

        return


    def set_meta_file_prefix(self, meta_file_prefix):
        ''' Set meta file prefix. '''

        self.meta_file_prefix = meta_file_prefix

        return


    def update_export_paths(self):
        '''
        Update/reset export paths.
        The function can also be used for initial set of export paths.
        '''
        self.export_path = self.get_struct_coll_path()+self.species+'_'+ \
                           self.species_name+self.export_path_postfix+'/'
        filesystem.create_folder([self.export_path])


        self.export_path_ax_predictions = self.export_path+ \
                                    'snapshot1_ax_predictions/'

        self.export_path_snapshot = self.export_path+ \
                                    'snapshot2_separated_chains/'

        self.export_path_pure_superposition = self.export_path+ \
                                    'snapshot3_pure_superposition/'

        self.export_path_flattened = self.export_path+ \
                                    'snapshot4_flattened/'

        self.export_path_tile = self.export_path+ \
                                     'tile/'

        self.export_path_primitive_unit_cell = self.export_path+ \
                                     'primitive_unit_cell/'

        self.export_path_assembly = self.export_path+ \
                                     'assembly/'

        self.export_path_complete_chains = self.export_path+ \
                                     'complete_chains/'


        self.layer_raw_path = self.export_path+ \
                self.export_file_prefix+'_'+self.version+'.cxs'
        self.layer_aligned_raw_path = self.export_path+ \
                self.export_file_prefix+'_aligned_'+self.version+'.cxs'
        self.layer_snapin_raw_path = self.export_path+ \
                self.export_file_prefix+'_snapin_'+self.version+'.cxs'
        self.layer_complete_chains_raw_path = self.export_path+ \
                self.export_file_prefix+'_complete_chains_'+self.version+'.cxs'
        self.layer_primitive_unit_cell_raw_path = self.export_path+ \
                self.export_file_prefix+'_primitive_unit_cell_'+ \
                self.version+'.cxs'


        self.layer_raw_path = self.export_path+self.export_file_prefix+'_'+ \
                              self.version+'.cxs'

        return


    def set_gene(self, gene):
        ''' Set gene name. '''

        self.gene = gene
        self.export_file_prefix = self.species+'_'+self.gene
        self.export_path_filters = [self.export_file_prefix]

        preassemblies_dir = control_utils.get_preassemblies_dir( \
                                                self, self.execution_dir)
        self.set_preassemblies_dir(preassemblies_dir)
        self.update_export_paths()

        return


    def set_domains(self):
        ''' Set domain boundaries. '''

        self.domains = domains

        return


    def import_domains(self):
        ''' Import domain boundaries from '_d.fa'-file. '''

        path_fasta = self.path_ax_predictions+ \
                                 self.species+'/'+self.symplex_path
        fasta_file = sorted(glob.glob(path_fasta+self.gene+'_d.fa'))

        if len(fasta_file) != 1:
            ctl.e(fasta_file)
            ctl.e(self.gene+'_d.fa')
            ctl.error('import_domains: no unique fasta file (..._d.fa).')

        domains, domain_boundaries = bibfasta.get_domains(fasta_file[0])
        self.domains = domain_boundaries

        if os.path.exists(self.path_ax_predictions+self.species+'/'+ \
                                        self.symplex_path+'setting_meta.txt'):
            f = open(self.get_struct_coll_meta_path()+ \
                    self.export_file_prefix+'_domain_boundaries.txt', 'w')
            f.write(str(domain_boundaries)+"\r\n")  
            f.close()

        return domain_boundaries


    def set_preassemblies_dir(self, assemblies_dir):
        ''' Set preassemblies directory. '''

        self.path_ax_predictions = assemblies_dir
        self.update_export_paths()

        return


    def set_struct_coll_path(self, struct_coll_path, struct_coll_meta_path=''):
        ''' Set assemblies directory (for output structures). '''

        self.struct_coll_path = struct_coll_path
        self.update_export_paths()

        if struct_coll_meta_path != '':
            self.struct_coll_meta_path = struct_coll_meta_path

        return


    def get_struct_coll_path(self):
        ''' Get assemblies directory. '''

        if self.struct_coll_path != '':
            path = self.struct_coll_path
        else:
            if os.path.exists( \
                self.execution_dir+'setting_dir_superordinate.txt'):
                path = os.path.dirname(os.path.realpath(__file__))+ \
                       '/../../Structures_'+self.version+'/'
                path = filesystem.clean_path(path)
            else:
                path = self.execution_dir

        filesystem.create_folder([path])

        return path


    def get_struct_coll_meta_path(self):
        ''' Get path for meta data of created layer models. '''

        if self.struct_coll_meta_path != '':
            path = self.struct_coll_meta_path
        else:
            path = self.get_struct_coll_path()[0:-1]+'_meta/'

        if os.path.exists(self.path_ax_predictions+self.species+'/'+ \
                                        self.symplex_path+'setting_meta.txt'):
            filesystem.create_folder([path])

        return path


    def set_run_level(self, runlevel):
        ''' Set runlevel. '''

        self.runlevel = runlevel

        return


    def get_run_level():
        ''' Get runlevel. '''

        ret = 0

        try:
            if run_level:
                ctl.d('run_level set')
                ret = run_level
        except NameError:
           ctl.d('run_level not set')
           ret = 0

        return ret


    def delete_rawfiles(self):
        ''' Delete rawfiles. '''

        try:
            os.unlink(self.layer_raw_path)
        except IOError:
            pass
        try:
            os.unlink(self.layer_aligned_raw_path)
        except IOError:
            pass
        try:
            os.unlink(self.layer_snapin_raw_path)
        except IOError:
            pass
        try:
            os.unlink(self.layer_complete_chains_raw_path)
        except IOError:
            pass

        return
