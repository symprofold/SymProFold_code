import copy
import glob

import ctl
import bibfasta
import subchain_set


class Assemblies():
    '''
    This class describes the list of all assemblies.
    Functions to handle the data structure of all assemblies are provided.

    list for each assembly
          0:species name, 1:gene id, 2:symm group
          3:main folder of predictions
          4:include in combined plot [0:no, 1:yes,part1, 2:yes,part2]
          5:large plot for full ensemble w/o interface clustering
                  [0:no, 1:yes]
          6:filter params [lowest ranking score]
          7:check params [max number of models per prediction scenario]
          8:display params [number of clusters in reduced mode,
                  [ranks of clusters to exclude]]
          9:clustering_mode, 0:no clustering, 1:clustering
    '''


    def __init__(self, assemblies):
        ''' Initialization of the Assemblies class. '''

        self.assemblies = assemblies
        self.filter_params_init()

        return


    def init_lists(self, assembly):
        '''
        Create lists for the data structure of a assembly.
        '''
        assembly['subchains'] = {}

        return assembly


    def filter_params_init(self):
        '''
        Initialize filter parameters for filtering by score.
        '''
        for i,s in enumerate(self.assemblies):
            if s[6] == []:
                s[6] = [0]

        return


    def add_domain_info(self):
        '''
        Add domain info for each assembly.
        '''
        for i,s in enumerate(self.assemblies):
            files0 = sorted(glob.glob(self.assemblies[i][3]+ \
                                      self.assemblies[i][1]+'*_d.fa'))

            if len(files0) == 0:
                ctl.e(assemblies[i][3]+assemblies[i][1]+'*_d.fa')
                ctl.e(s)
                ctl.e(files0)
                ctl.error('Assemblies: add_domain_info: '+ \
                          'domain fasta file not found')

            if len(files0) > 1:
                ctl.e(assemblies[i][3]+assemblies[i][1]+'*_d.fa')
                ctl.e(s)
                ctl.e(files0)
                ctl.error('Assemblies: add_domain_info: '+ \
                          'more than 1 domain fasta file')

            domains, domain_boundaries = bibfasta.get_domains(files0[0], False)
            self.assemblies[i]['dom_n'] = len(domains)
            self.assemblies[i]['domains'] = domains
            self.assemblies[i]['domain_boundaries'] = domain_boundaries

        return self.assemblies


    def set_assembly_labels(self, reduced_mode):
        '''
        Set labels for each assembly.
        '''
        for i_,s in enumerate(self.assemblies):
            if reduced_mode == 1:
                s['assembly_label'] = s[0].replace(' ', '\n')+'\n'+ \
                                     '('+s[1]+')'+'\n'+s[2]
            else:
                s['assembly_label'] = s[0].replace('- ', '')+' '+ \
                                     '('+s[1]+')'+'\n'+s[2]+'\n'

        return self.assemblies


    def add_subchains(self):
        '''
        Add standard set of subchains for each assembly.
        '''
        for i,s in enumerate(self.assemblies):
            s['subchain_set'] = subchain_set.get_subchain_set(\
                                    self.assemblies[i][3]+'msa_gen/'+ \
                                    self.assemblies[i][1]+'_'+'*.fa')

            if s['subchain_set'] == []:
                ctl.error('Assemblies: add_subchains: '+ \
                          'standard subchain set not found')

        return self.assemblies


    def reduce_cluster_n(self, max_cluster_n, verbous=False):
        '''
        Reduce number of interface clusters to upper limit.
        '''
        assemblies_red = copy.deepcopy(self.assemblies)
        filled_rows = [[] for i in self.assemblies]

        for i,a in enumerate(assemblies_red):
            if assemblies_red[i][8] != []:
                max_cluster_n_ = assemblies_red[i][8][0]
            else:
                max_cluster_n_ = max_cluster_n

            # iterate through subchains_clustered
            for j,sc in enumerate(a['subchains_clustered']):
                if verbous:
                    ctl.p(a['subchains_clustered'][sc])

                # iterate through clusters
                for k,cluster in enumerate(a['subchains_clustered'][sc]):

                    # exclude clusters as specified by user
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    if len(assemblies_red[i][8]) >= 2:
                        if k in assemblies_red[i][8][1]:
                            continue

                        if not k < max_cluster_n_+len(assemblies_red[i][8][1]):
                            continue

                    else:
                        if k >= max_cluster_n_:
                            continue

                    cluster_cleaned = []

                    # iterate through SymPlexes of one cluster
                    for l,symplex in enumerate(cluster):
                        filled_rows[i].append(k)
                        cluster_cleaned.append(symplex)

                        if verbous:
                            ctl.p('cluster_cleaned')
                            ctl.p(cluster_cleaned)

                    a['subchains_clustered'][sc][k] = cluster_cleaned


        # clean empty cluster rows
        # ~~~~~~~~~~~~~~~~~~~~~~~~
        for i,a in enumerate(assemblies_red):

            # iterate through subchains_clustered
            for j,sc in enumerate(a['subchains_clustered']):
                sc_cleaned = []

                # iterate through clusters
                for k,cluster in enumerate(a['subchains_clustered'][sc]):

                    if k in filled_rows[i]:
                        sc_cleaned.append(cluster)

                a['subchains_clustered'][sc] = sc_cleaned

      
        self.assemblies = assemblies_red

        return self.assemblies


    def reduced_mode(self):
        '''
        Reduce list for reduced mode.
        reduced mode: only one prediction per prediction scenario
        '''
        assemblies_reduced = copy.deepcopy(self.assemblies)

        # empty list
        for i,a in enumerate(self.assemblies):

            # iterate through subchains_clustered
            for j,sc in enumerate(a['subchains_clustered']):

                # iterate through rot axes
                for k,rot_ax in enumerate(a['subchains_clustered'][sc]):
                    assemblies_reduced[i]['subchains_clustered'][sc][k] = []


        # iterate through assemblies
        for i,a in enumerate(self.assemblies):

            # iterate through subchains_clustered
            for j,sc in enumerate(a['subchains_clustered']):

                # iterate through rot axes
                for k,rot_ax in enumerate(a['subchains_clustered'][sc]):
                    rot_ax = sorted(rot_ax, key=lambda x: x[1], reverse=True)

                    if rot_ax != []:
                        rot_ax_ = [rot_ax[0]]
                        assemblies_reduced[i]['subchains_clustered'][sc] \
                                                                [k] += rot_ax_

        self.assemblies = assemblies_reduced

        return self.assemblies


    def reduce_to_plottype(self, path, mode):
        '''
        Reduce list dependent on type of plot.
        Types of plot: specific species, general plot with many species
        '''
        assemblies = []

        if mode.plottype == 0:
            # Check if the script was called from a specific species path that
            # is registered in config.
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            for s in self.assemblies:
                if s[3] == path:
                    assemblies.append({})

                    for i,data in enumerate(s):
                        assemblies[0][i] = s[i]

                    assemblies[0] = self.init_lists(assemblies[0])
                    break


            # If not called from a path registered in config, try automatic
            # assignment to specific species.
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if assemblies == []:
                spec_shortname = path.split('/')[-2]
                files0 = sorted(glob.glob(path+'*_d.fa'))
                gene_id = files0[0].split('/')[-1].split('_d.fa')[0]
                
                s =  [spec_shortname, gene_id, 'p?', \
                     path, 0, 0, [0], [5], [], 1]
                assemblies.append({})

                for i,data in enumerate(s):
                    assemblies[0][i] = s[i]
                assemblies[0] = self.init_lists(assemblies[0])


        elif mode.plottype == 1:
            # Prepare data for general plot with many species.
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            for s in self.assemblies:
                if s[4] == mode.reduced_part or \
                   (mode.reduced_showall and s[4] != 0):
                    assemblies.append({})

                    for i,data in enumerate(s):
                        assemblies[-1][i] = s[i]
                assemblies[-1] = self.init_lists(assemblies[-1])


        else:
            ctl.error('Assemblies: reduce_to_plottype')

        self.assemblies = assemblies

        return self.assemblies


    def get_cluster_col_widths(self, max_multiplicity_n=10, max_cluster_n=100):
        '''
        Determine maximal width of each column  of a table with clusters as
        columns (x axis) and multiplicities as rows (y axis).
        '''
        col_widths = []
        cluster_cell_widths = self.get_cluster_cell_widths( \
                                    max_multiplicity_n, max_cluster_n)

        for i,c in enumerate(cluster_cell_widths):
            col_widths.append(max(c))

        return col_widths


    def get_cluster_cell_widths(self, max_multiplicity_n=10, max_cluster_n=100):
        '''
        Determine width of each cell of a table with clusters as columns
        (x axis) and multiplicities as rows (y axis).
        '''
        cell_widths = [[ 0 for multiplicity in range(0, max_multiplicity_n) ] \
                      for cluster in range(0, max_cluster_n)]
            # multiplicity: order of rotational symm. axis
            # max_cluster_n: maximal number of clusters
        col_num_assembly = 0

        # iterate through assemblies
        for i,s in enumerate(self.assemblies):

            # iterate through subchains_clustered
            for j,sc in enumerate(s['subchains_clustered']):
                col_num = col_num_assembly

                # iterate through clusters
                for k,cluster in enumerate(s['subchains_clustered'][sc]):

                    # iterate through SymPlexes of a cluster
                    for l,symplex in enumerate(cluster):
                        cell_widths[col_num][symplex[0]] += 1

                    col_num += 1

            col_num_assembly = col_num

        return cell_widths


    def get_multiplicity(self, pred_scenario):
        '''
        Get multiplicity of a complex model from prediction scenario name.
        '''
        mult = int(pred_scenario.split('_')[-1].split('x')[1])

        return mult
