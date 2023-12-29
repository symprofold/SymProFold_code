import glob
import os

import ctl
import bibpdb
import symplot.clash_dist
import symplot.checks
import interface_cluster
import interface_matrix
import interface_matrix_signed
import symplot.prediction_scenario
import metadata
import model_filter.model_filter
import molmodel
import params
import symplot
import symplot.predictions_sort


'''
Module providing functions to determine interface clusters for each assembly.
'''


def determine_clusters(mode, subchain_mode, root_path, assemb, verbous):
    '''
    Determine interface clusters for each assembly.

    Steps: filter predicted models, load interface matrices,
    calculate correlation coefficients,
    determine interface clusters (cluster by interfaces).

    Returns:
        assemblies_interface_clusters: list of interface clusters for each
                assembly
        species_interface_distograms: list of interface distograms for each
                assembly
        assemblies_prediction_scernario_n: number of predictions for each
                prediction scenario (species/subchain/multiplicity),
                list for each assembly
        g0: graph
    '''

    # creation of lists
    # ~~~~~~~~~~~~~~~~~
    assemblies_interface_clusters = \
                            [[[], [], [], []] for a in assemb.assemblies]
    species_interface_distograms = [{} for a in assemb.assemblies]

    assemblies_prediction_scernario_n = [{} for a in assemb.assemblies]
        # number of subchain predictions for each prediction scenario
        # (species/subchain/multiplicity)
        # assemblies_prediction_scernario_n[species][subchain][multiplicity]


    for i_,a in enumerate(assemb.assemblies):

        # dictionary for interface distograms
        interface_distograms = {}
        interface_distograms2 = {}

        # get folders of all prediction scenarios
        scenario_folders = symplot.prediction_scenario.get_folders(a[3]+a[1])

        for sc_folder in scenario_folders:
            multiplicities = params.get_multiplicities(root_path, \
                                                       mode.microorganism_type)

            prediction_scenario = metadata.get_prediction_scenario(sc_folder)
            a['subchains'][prediction_scenario] = \
                                [[] for m in range(0, max(multiplicities)+1)]

            all_predictions = sorted(glob.glob(sc_folder+'/*.pdb'))
            all_predictions = symplot.predictions_sort. \
                                                sort_by_score(all_predictions)

            # check maximal models per prediction scenario
            symplot.checks.check_max_models_per_prediction_scenario( \
                                                            a, all_predictions)

            if subchain_mode == 1: # 1:standard subchain set
                for j,f in enumerate(all_predictions):
                    score = metadata.get_score(f)

                    if score != -1:
                        sc_abbr = metadata.get_subchain_abbr(f.split('/')[-2])

                        if sc_abbr in a['subchain_set']:
                            if not sc_abbr in \
                               assemblies_prediction_scernario_n[i_]:
                                assemblies_prediction_scernario_n[i_] \
                                    [sc_abbr] = \
                                    [0 for m in \
                                    range(0, max(multiplicities)+1)]

                            if verbous:
                                ctl.p(a['subchains'])

                            if assemb.get_multiplicity(prediction_scenario) \
                               in multiplicities:
                                # check if multiplicity is in multiplicities
                                # list (multiplicities to be analyzed,
                                # e.g. for this microorganism_type)
                                assemblies_prediction_scernario_n[i_] \
                                    [sc_abbr][assemb.get_multiplicity( \
                                    prediction_scenario)] += 1
                            else:
                                ctl.e('determine_clusters: '+ \
                                    'multiplicity not in multiplicities list')


            symmetries = sorted(glob.glob(sc_folder+'/'+'*symm_*'))

            for symm in symmetries:
                files = sorted(glob.glob(symm+'/clashes/*.pdb'))
                files = symplot.predictions_sort.sort_by_score(files)

                if subchain_mode == 1:
                    for j,f in enumerate(files):
                        score = metadata.get_score(f)
                        if score != -1:
                            sc_abbr = metadata.get_subchain_abbr( \
                                                            f.split('/')[-4])


                # filter
                # ------

                # filter by clashes
                # ~~~~~~~~~~~~~~~~~
                if mode.filter_clashes == 1:
                    files = model_filter.model_filter.filter_clashes(files)


                # filter by rotational symmetry
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if mode.filter_angles:
                    files = model_filter.model_filter.filter_rotang( \
                                files, multiplicities, mode.microorganism_type)


                # filter by the property that parts of the protein chain run
                # erroneously through the sphere of a another monomer
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if mode.filter_intermol_chain:
                    files, filtered_out = model_filter.model_filter. \
                                            filter_intermol_betasheet(files)

                    for f in filtered_out:
                        if verbous:
                            ctl.p('filtered_out:')
                            ctl.p(f)

                        try:
                            os.mkdir(root_path+'filter_intermol_sheet/')
                        except IOError:
                            pass

                        filename_meta = '__'.join(f.split('/')[-4:])
                        f = open(root_path+'filter_intermol_sheet/'+ \
                                filename_meta, 'w')
                        f.write('')
                        f.close()


                for j,f in enumerate(files):
                    symm_order = metadata.get_symm_order(f)
                        # order of rotational symmetry in the SymPlex
                    if symm_order == -1:
                        continue

                    score = metadata.get_score(f)

                    # filter by score: filter out if the score is < 0.2
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    if score < 0.20:
                        continue

                    clashes = model_filter.model_filter.clashes_converted(f)

                    filename_cl = f.replace('clashes/', ''). \
                                  replace('symm_180/', ''). \
                                  replace('symm_120/', ''). \
                                  replace('symm_090/', ''). \
                                  replace('symm_072/', ''). \
                                  replace('symm_060/', ''). \
                                  replace('symm_051/', '')


                    # get optional score from clash distribution file if
                    # available
                    if 'cla' in filename_cl:
                        filename_cl = filename_cl.split('_cla')[0]+'.pdb'

                    cl_dist = symplot.clash_dist.load( \
                                                    filename_cl+'_cldist.txt')
                    cl_score = symplot.clash_dist.get_score(cl_dist)


                    # load interface matrices
                    # -----------------------

                    # determine termini
                    rmsds = bibpdb.get_rmsds(f)
                    termini = molmodel.get_termini(rmsds)

                    # create list of terminus residues to exclude
                    terminus_res = [r for r in range(0, termini[0])]+ \
                                   [r for r in range(termini[1],5000)]

                    if mode.interface_matrix_type == 0:
                        distogram = interface_matrix.load(f, terminus_res)
                    elif mode.interface_matrix_type == 1:
                        distogram = interface_matrix_signed.load(f, \
                                                                 terminus_res)

                    if distogram == {}:
                        if verbous:
                            ctl.p(f)
                            ctl.p('distogram == {}')
                    else:
                        filename_ = f.split('/')[-1]
                        interface_distogram_key = \
                            (prediction_scenario, symm_order, score, filename_)

                        if interface_distogram_key in interface_distograms:

                            # check if distograms are identical
                            if interface_distograms[ \
                                    interface_distogram_key] != distogram:
                                ctl.e(f)
                                ctl.e(distogram)
                                ctl.e(interface_distograms[ \
                                        interface_distogram_key])
                                ctl.error('get_symplex_data: '+ \
                                        'different distograms')

                        else:
                            subchain_abbr = metadata.get_subchain_abbr( \
                                                    prediction_scenario)

                            # cluster if subchain_mode == 0 or
                            # if subchain in subchain set
                            if subchain_mode == 0 or \
                               subchain_abbr in a['subchain_set']:
                                interface_distograms[ \
                                        interface_distogram_key] = distogram
                                interface_distograms2[ \
                                        interface_distogram_key] = distogram


                    if score != -1 and score >= a[6][0]:
                        if subchain_mode == 1 or subchain_mode == 0:
                            cluster_name = None

                            # add SymPlex to SymPlot if subchain_mode == 0 or
                            # if subchain in subchain set
                            if subchain_mode == 0 or \
                                sc_abbr in a['subchain_set']:
                                a['subchains'][prediction_scenario] \
                                    [symm_order]. \
                                    append([score, cluster_name, 0, \
                                            cl_score, clashes, f])


        species_interface_distograms[i_] = interface_distograms


        # calculate correlation coefficients
        # ----------------------------------
        if mode.interface_matrix_type == 0:
            interface_correlations_ = interface_matrix. \
                                    corr_coefficients(interface_distograms2)
        elif mode.interface_matrix_type == 1:
            interface_correlations_ = interface_matrix_signed. \
                                    corr_coefficients(interface_distograms2)

        f = open(root_path+'interface_correlation.txt', 'w')
        f.write('')
        f.close()

        f = open(root_path+'interface_correlation.txt', 'a')

        for ic in interface_correlations_:
            f.write(str(ic)+"\r\n")
            f.write(str(interface_correlations_[ic])+"\r\n")

        f.close()


        # determine interface clusters (cluster by interfaces)
        # ----------------------------------------------------
        if verbous:
            print('determine interface clusters')

        for interface_cluster_mode in mode.interface_cluster_modes:
            interface_cluster0, g0, cluster_n0 = interface_cluster.cluster( \
                        interface_correlations_, interface_cluster_mode)

            # sort clusters by highest score
            interface_cluster0 = interface_cluster.sort_clusters( \
                                            interface_cluster0, cluster_n0)

            if a[9] == 0:
                for ic in interface_cluster0:
                    interface_cluster0[ic] = 0

                cluster_n0 = 1

                if verbous:
                    ctl.p(interface_cluster0)


            assemblies_interface_clusters[i_][interface_cluster_mode] = \
                        [interface_cluster0, g0, cluster_n0]


    return assemblies_interface_clusters, species_interface_distograms, \
                assemblies_prediction_scernario_n, g0
