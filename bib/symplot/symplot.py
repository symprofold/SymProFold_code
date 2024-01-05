import glob
import os
import sys
import time

import matplotlib.pyplot as plt

import ctl
import bibpdb
import diagram.marker
import filesystem
import interface.binding_domains
import interface_cluster
import metadata
import symplot.assemblies
import symplot.assemblies_clusters
import symplot.bib
import symplot.checks
import symplot.nnptm
import diagram.colors
import params
import symplot.config
import symplot.export
import symplot.layout
import symplot.legend
import symplot.prediction_scenario
import symplot.top_labels
import symplot.x_axis_labels


'''
Module providing functions to calculate SymPlots.
'''

def calc_symplot(root_path, mode, assemb, symplot_):
    '''
    Calculate SymPlot.
    '''
    verbous0 = mode.verbous

    # interface correlations
    # ----------------------

    for subchain_mode in mode.subchain_modes:

        # determine interface clusters for each assembly
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # steps: filter predicted models, load interface matrices,
        # calculate correlation coefficients,
        # determine interface clusters (cluster by interfaces)
        assemblies_interface_clusters, species_interface_distograms, \
            assemblies_prediction_scenario_n, g0 = \
                    symplot.assemblies_clusters.determine_clusters( \
                            mode, subchain_mode, root_path, assemb, verbous0)


        # create list of all occuring interface resids
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        assemblies_interfaces_txt = [ ['', ['' for i in range(0,100)]] \
                                        for s in assemb.assemblies]

        if mode.verbous:
            ctl.p('clusters')
            ctl.p(assemblies_interface_clusters[0])

        plot_cluster_fn = root_path+assemb.assemblies[0][1]+'_cluster_'+ \
                                            str(round(time.time()))+'.png'
        interface_cluster.plot_clusters(g0, plot_cluster_fn, plt)


        # check completeness of prediction scenarios
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if subchain_mode == 1:
            output_prediction_completeness = \
                symplot.checks.check_prediction_completeness(assemb, \
                                            assemblies_prediction_scenario_n)


        # generate subchain info box
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~
        assembly_id = 0
            # no iterations, because single plot assumed, index 0 is used as
            # assembly index

        subchain_info = ''

        if subchain_mode == 1:
            for sp in assemblies_prediction_scenario_n[assembly_id]:
                subchain_info += sp+'\u2006·\u2006('

                for i,m in enumerate(assemblies_prediction_scenario_n \
                                                         [assembly_id][sp]):
                    if m > 0:
                        subchain_info += str(i)+','

                subchain_info = subchain_info[:-1]
                subchain_info += ')\n'


        for show_main_bd in mode.show_main_bd_modes:
            if show_main_bd == 1:
                symplot_.plot_height = mode.plot_height_main_bd
                symplot_.subplots_adjust_top = 0.90
            else:
                symplot_.subplots_adjust_top = None

            if mode.verbous:
                ctl.p('subchain_info')
                ctl.p(subchain_info)


            for interface_cluster_mode in mode.interface_cluster_modes:

                # clustering by interfaces
                # ~~~~~~~~~~~~~~~~~~~~~~~~
                for i,a in enumerate(assemb.assemblies):
                    subchains_clustered = {}
                    cluster_n = assemblies_interface_clusters[i] \
                                        [interface_cluster_mode][2]

                    for sc in a['subchains']:
                        subchains_clustered[sc] = \
                                            [[] for c in range(0, cluster_n)]

                    for sc in a['subchains']:

                        # iterate through multiplicities
                        for j,mult in enumerate(a['subchains'][sc]):
                            for symplex in mult:
                                filename_ = symplex[5].split('/')[-1]

                                rank_id = int(filename_. \
                                        replace('unrelaxed_rank', ''). \
                                        replace('rank', '').split('_')[0])

                                node_id = sc.split('/')[-1]+'_'+ \
                                              str(symplex[0])+'_'+str(rank_id)

                                if node_id in assemblies_interface_clusters \
                                                [i][interface_cluster_mode][0]:
                                    cluster_id = \
                                        assemblies_interface_clusters[i] \
                                        [interface_cluster_mode][0][node_id]

                                    if mode.verbous:
                                        ctl.p(cluster_id)

                                    if mode.verbous:
                                        ctl.p(subchains_clustered[sc])
                                        ctl.p(subchains_clustered[sc] \
                                                        [cluster_id])
                                        ctl.p(symplex)

                                    subchains_clustered[sc][cluster_id]. \
                                        append([j, symplex[0], symplex[1], \
                                            symplex[2], symplex[3], \
                                            symplex[4], symplex[5]])


                    a['subchains_clustered'] = subchains_clustered

                    if mode.verbous:
                        ctl.p('a[subchains_clustered]')
                        ctl.p(a['subchains_clustered'])


                fig, ax = plt.subplots(figsize=(symplot_.plot_width, \
                                                symplot_.plot_height))

                # set margin parameters in Figure object
                fig = symplot_.subplots_adjust(fig)

                cluster_total = assemblies_interface_clusters[0] \
                                                [interface_cluster_mode][2]

                if mode.reduced == 1:
                    # reduce number of interface clusters to upper limit
                    assemb.reduce_cluster_n(mode.reduced_cluster_n)
                    cluster_total = len(assemb.assemblies)* \
                                    mode.reduced_cluster_n

                # set font sizes
                symplot_.set_ax_fontsizes(cluster_total)


                # generate x axis (bottom) labels for each cluster, w/o spacer
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if mode.reduced == 1:
                    col_labels0 = symplot.x_axis_labels. \
                                        x_axis_labels_interfacecluster( \
                                                assemb, mode.reduced_cluster_n)
                else:
                    col_labels0 = symplot.x_axis_labels. \
                                        x_axis_labels_interfacecluster( \
                                                assemb, cluster_n)


                # reduce assembly list for reduced mode
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if mode.reduced == 1:
                    assemb.reduced_mode()

                data_points_per_col = 1
                    # layout precision: number of datapoints per fig column

                # determine column labels (cluster labels) with spacer
                col_labels, col_w, col_centerpos = symplot.x_axis_labels. \
                                x_axis_labels_w_spacher_interfacecluster( \
                                    assemb, col_labels0, data_points_per_col)

                if mode.verbous:
                    ctl.p('col_w')
                    ctl.p(col_w)
                    ctl.p('cluster_n')
                    ctl.p(cluster_n)


                # determine list assembly_label to use as label on top of plot
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if mode.reduced == 1:
                    assembly_label = symplot.top_labels.assembly_labels( \
                        assemb, mode.reduced_cluster_n, col_w)
                else:
                    assembly_label = symplot.top_labels.assembly_labels( \
                        assemb, cluster_n, col_w)


                # determine label font sizes and circle size for data points
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                col_widths = assemb.get_cluster_col_widths()
                        # get maximal width for each column

                diagram_col_n = sum(col_widths)
                symplot_.set_datapoint_fontsizes(diagram_col_n)


                # prepare data for plot (plot_data and plot_data_sort)
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # plot_data, if mode.sort_scores == 0
                # plot_data_sort, if mode.sort_scores == 1
                plot_data = [{} for i in range(0, 100)]
                plot_data_sort = [[] for i in range(0, 100)]

                step_fact = (1/data_points_per_col)*0.95
                vline = [0]

                col_num_sp = 0
                col_count_sp = 0

                # iterate through assemblies
                for assembly_id, s in enumerate(assemb.assemblies):

                    # pos_register for whole plot with all displayed assemblies
                    pos_reg = [[ 0 for j in range(0,10) ] \
                                                        for i in range(0,100)]


                    # fill list of all occuring interface resids
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    species_txt = s['assembly_label'].split('\n')

                    if mode.verbous:
                        ctl.p(s['assembly_label'])
                        ctl.p(species_txt)

                    assemblies_interfaces_txt[assembly_id][0] = \
                                species_txt[0].replace('- ', '')+"\r\n". \
                                        join(['=' for s in species_txt[0]. \
                                        replace('- ', '')])


                    # iterate through subchains_clustered
                    for i,sc in enumerate(s['subchains_clustered']):
                        col_num = col_num_sp
                        col_count = col_count_sp

                        # iterate through subchain clusters
                        for k,cluster in enumerate( \
                                                s['subchains_clustered'][sc]):
                            if mode.verbous:
                                ctl.p(s['subchains_clustered'][sc])

                            col_num = col_count+col_centerpos[col_count]
                            col_count += 1

                            if k == 0 and col_num not in vline:
                                ax.axvline( \
                                    x=(col_num-(1/data_points_per_col)), \
                                    color='black', lw=0.75)

                            elif col_num not in vline:
                                ax.axvline( \
                                    x=(col_num-(1/data_points_per_col)), \
                                    color='black', lw=0.7, linestyle='dotted')
                                vline.append(col_num)


                            # create lists of datapoint data to include in plot
                            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            x = [] # x coords of datapoints in diagram
                            y = [] # y coords of datapoints in diagram
                            mol_count = []
                                # molecule counts of datapoints in diagram
                            score = []
                                # scores of datapoints in diagram
                            cl_score = []
                                # clash scores of datapoints in diagram
                            filenames = []
                                # filenames of datapoints in diagram


                            # iterate through SymPlexes of one subchain cluster
                            for symplex_numb,symplex in enumerate(cluster):
                                symplex_score = symplex[1]
                                clash_score = symplex[4]

                                if k >= len(pos_reg)-1:
                                    ctl.error('calc_symplot: '+ \
                                            'max len of pos_reg reached')

                                pos_reg[k][symplex[0]] += 1

                                # mode show_nnptm
                                infix_nnptm = ''
                                if mode.show_nnptm == 1:
                                    nnptm = symplot.nnptm.get_nnptm( \
                                                        symplex[6], root_path)
                                    if nnptm != -1:
                                        infix_nnptm = str(nnptm)
                                        symplex_score = nnptm
                                    else:
                                        symplex_score = 0.0
                                    

                                x.append(col_num+step_fact*( \
                                                pos_reg[k][symplex[0]]-1) )
                                y.append(symplex[0])
                                score.append((symplex_score**2)* \
                                                symplot_.datapoint_circle_size)
                                filenames.append(symplex[6])

                                v_dist = 0.08+symplex_score*0.20* \
                                        (symplot_.datapoint_circle_size/500)

                                infix_clash_score = ''

                                if mode.show_clash_score == 1:
                                    infix_clash_score = '\n'+str(clash_score)
                                elif mode.show_clash_score == 2:
                                    clashes = symplex[5]
                                    if clashes > 3:
                                        infix_clash_score = '\n'+'⚡'


                                if show_main_bd == 1:
                                    if verbous0:
                                        ctl.p(s['domains'])
                                        ctl.p(s['domain_boundaries'])

                                    interface_distogram_found = 0
                                    interface_distogram_ = []

                                    for interface_distogram in \
                                                species_interface_distograms \
                                                                [assembly_id]:
                                        symplex_filename_part1 = \
                                                    symplex[6].split('/')[-1]

                                        if interface_distogram[0] == sc and \
                                           interface_distogram[1] == \
                                                        symplex[0] and \
                                           interface_distogram[2] == \
                                                        symplex[1] and \
                                           interface_distogram[3] == \
                                                        symplex_filename_part1:

                                            if interface_distogram_found > 0:
                                                ctl.e(interface_distogram)
                                                ctl.e(symplex)
                                                ctl.error('calc_symplot: '+ \
                                                    'interface_distogram_'+ \
                                                    'found > 0')

                                            interface_distogram_found += 1

                                            if verbous0:
                                                ctl.p(symplex)
                                                ctl.p(interface_distogram)

                                            interface_distogram_ = \
                                                species_interface_distograms \
                                                [assembly_id] \
                                                [interface_distogram]
                                            zero, intres = \
                                                interface.binding_domains. \
                                                get_interface_residue_ranges( \
                                                        interface_distogram_)


                                    if interface_distogram_found != 1:
                                        ctl.e(interface_distogram_found)
                                        ctl.e(interface_distogram)
                                        ctl.e(symplex)
                                        ctl.error('calc_symplot: '+ \
                                            'interface_distogram_found != 1')

                                    if len(interface_distogram_) < 1:
                                        ctl.error('calc_symplot: '+ \
                                            'len(interface_distogram_) < 1')

                                    domain_bindings = \
                                            interface.binding_domains. \
                                                get_domain_bindings( \
                                                    s['domain_boundaries'], \
                                                interface_distogram_)

                                    if verbous0:
                                        ctl.p('domain_bindings')
                                        ctl.p(domain_bindings)


                                    infix_mainbd = ''
                                    parts = 0

                                    for db_i,db in enumerate(domain_bindings):
                                        if parts > 0 and db[2] < 10:
                                            continue

                                        if parts >= 2:
                                            break

                                        parts += 1
                                        infix_mainbd += str(db[0]+1)+'·'+ \
                                                        str(db[1]+1)+',  '


                                    infix_mainbd = infix_mainbd[:-3]

                                    if len(domain_bindings) > 15:
                                        ctl.p('more than 15 domain_bindings')


                                    subchaindef = sc.split('_')[-1].split('x')[0]
                                    subchains_ = metadata. \
                                                subchainabbr_to_subchains( \
                                                        subchaindef)


                                    if len(subchains_) == 0:
                                        subchain_short = \
                                                'FL·'+sc.split('_')[-1]. \
                                                split('x')[1]
                                    elif len(subchains_) == 1:
                                        subchain_short = \
                                                str(subchains_[0])+'·'+ \
                                                sc.split('_')[-1].split('x')[1]
                                    elif len(subchains_) == 2:
                                        subchain_short = \
                                                str(subchains_[0])+''+ \
                                                str(subchains_[1])+'·'+ \
                                                sc.split('_')[-1].split('x')[1]


                                    assemblies_interfaces_txt \
                                        [assembly_id][1][k] += \
                                                'symmetry complex '+ \
                                                    subchain_short+', '+ \
                                                'ipTM+pTM: '+ \
                                                    str(symplex[1])+", "+ \
                                                'interface residues: '+ \
                                                    intres+"\r\n\r\n"


                                if symplex_score >= 0.20:
                                    if mode.show_intraplot_labels:
                                        if show_main_bd != 1:
                                            ax.text( \
                                                col_num+step_fact* \
                                                    (pos_reg[k][symplex[0]]- \
                                                    1), \
                                                symplex[0]-v_dist, \
                                                '.'+ \
                                                    (str(round( \
                                                    symplex_score, 2)). \
                                                    split('.')[1]+'00')[:2]+ \
                                                    infix_clash_score, \
                                                ha='center', va='top', \
                                                fontsize=symplot_. \
                                                    datapoint_fontsize_bottom)
                                elif symplex_score == 0.0:
                                    ax.text(col_num+step_fact* \
                                        (pos_reg[k][symplex[0]]-1), \
                                        symplex[0]-v_dist, \
                                        '.'+'\n'+ \
                                            infix_clash_score+infix_mainbd, \
                                        ha='center', va='top', \
                                        fontsize= \
                                            symplot_.datapoint_fontsize_bottom)

                                subchaindef = sc.split('_')[-1].split('x')[0]
                                subchains_ = metadata. \
                                        subchainabbr_to_subchains(subchaindef)


                                if len(subchains_) == 0:
                                    subchain_short = \
                                        'F\nL\n·'+ \
                                            sc.split('_')[-1].split('x')[1]
                                elif len(subchains_) == 1:
                                    subchain_short = \
                                        ' \n'+str(subchains_[0])+'\n·'+ \
                                            sc.split('_')[-1].split('x')[1]
                                elif len(subchains_) == 2:
                                    subchain_short = \
                                        str(subchains_[0])+'\n'+ \
                                            str(subchains_[1])+'\n·'+ \
                                            sc.split('_')[-1].split('x')[1]

                                mol_count.append( \
                                        int(sc.split('_')[-1].split('x')[1]))

                                if mode.show_intraplot_labels and \
                                                            show_main_bd != 1:
                                    ax.text(col_num+step_fact* \
                                        (pos_reg[k][symplex[0]]-1), \
                                        symplex[0]+v_dist, \
                                        subchain_short, \
                                        ha='center', va='bottom', \
                                        fontsize=symplot_.datapoint_fontsize, \
                                        linespacing= \
                                            symplot_.datapoint_linespacing)
                                elif mode.show_intraplot_labels and \
                                                            show_main_bd == 1:
                                    ax.text(col_num+step_fact* \
                                        (pos_reg[k][symplex[0]]-1), \
                                        symplex[0]+v_dist, \
                                        infix_mainbd, \
                                        ha='center', va='bottom', \
                                        fontsize=symplot_.datapoint_fontsize, \
                                        linespacing= \
                                            symplot_.datapoint_linespacing, \
                                        rotation=90)


                                if sc.split('_')[-1].split('x')[0] == '':
                                    subchain_name = 'full length'            
                                else:
                                    subchain_name = 'subchain '+ \
                                                sc.split('_')[-1].split('x')[0]

                                subchain_name_mult = 'subchain '+ \
                                                sc.split('_')[-1]

                                subchain_mult = \
                                        int(subchain_name_mult.split('x')[-1])


                                # insert plot entries for each symplex_numb
                                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                if symplex_numb >= len(plot_data):
                                    ctl.e(symplex_numb)
                                    ctl.e(len(plot_data))
                                    ctl.error('calc_symplot: '+ \
                                            'symplex_numb >= len(plot_data)')

                                if subchain_name in plot_data[symplex_numb]:
                                    plot_data \
                                        [symplex_numb][subchain_name][0] = \
                                        plot_data \
                                            [symplex_numb][subchain_name][0]+x
                                    plot_data \
                                        [symplex_numb][subchain_name][1] = \
                                        plot_data \
                                            [symplex_numb][subchain_name][1]+y
                                    plot_data \
                                        [symplex_numb][subchain_name][2] = \
                                        plot_data \
                                            [symplex_numb][subchain_name][2]+ \
                                            score
                                    plot_data \
                                        [symplex_numb][subchain_name][3] = \
                                        subchain_name
                                    plot_data \
                                        [symplex_numb][subchain_name][4] = \
                                        plot_data \
                                            [symplex_numb][subchain_name][4]+ \
                                            mol_count
                                    plot_data \
                                        [symplex_numb][subchain_name][5] = \
                                        plot_data \
                                            [symplex_numb][subchain_name][5]+ \
                                            filenames
                                else:
                                    plot_data[symplex_numb][subchain_name] = \
                                        [x, y, score, subchain_name, \
                                        mol_count, filenames]


                                if len(plot_data_sort[symplex_numb]) > 0:
                                    plot_data_sort[symplex_numb][0] = \
                                        plot_data_sort[symplex_numb][0]+x
                                    plot_data_sort[symplex_numb][1] = \
                                        plot_data_sort[symplex_numb][1]+y
                                    plot_data_sort[symplex_numb][2] = \
                                        plot_data_sort[symplex_numb][2]+score
                                    plot_data_sort[symplex_numb][3] = \
                                        subchain_name
                                    plot_data_sort[symplex_numb][4] = \
                                        plot_data_sort[symplex_numb][4]+ \
                                        mol_count
                                    plot_data_sort[symplex_numb][5] = \
                                        plot_data_sort[symplex_numb][5]+ \
                                        filenames
                                else:
                                    plot_data_sort[symplex_numb] = \
                                        [x, y, score, subchain_name, \
                                        mol_count, filenames]


                    col_num_sp = col_num
                    col_count_sp = col_count


                # fill plot_data and plot_data_sort into plot
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # plot_data, if mode.sort_scores == 0
                # plot_data_sort, if mode.sort_scores == 1
                colordict = {}
                zorder = 1000000
                marker_style = 1
                size_fact = 1.0
                subplots = []

                if mode.sort_scores == 0:
                    for symplex_numb in plot_data:
                        for sc in symplex_numb:
                            rotax_dataset = symplex_numb[sc]

                            for p,symplex_dp in enumerate(rotax_dataset[0]):
                                multiplicity = rotax_dataset[1][p]

                                scenario_path = symplot.prediction_scenario. \
                                    get_path(rotax_dataset[5][p], root_path)
                                filename = symplot.prediction_scenario. \
                                    get_coord_filename(rotax_dataset[5][p])
                                mol_count = bibpdb.get_multimer_n( \
                                            scenario_path+filename)

                                path_ = '/'.join(rotax_dataset[5][p]. \
                                                split('/')[:-1])+'/'
                                filename = rotax_dataset[5][p].split('/')[-1]
                                filename_part1 = filename.split('_cla')[0]

                                rotaxorder_files = sorted(glob.glob( \
                                        filesystem.clean_path( \
                                            path_+'../../'+filename_part1+ \
                                            '_rotaxorder*.txt')))

                                if len(rotaxorder_files) > 0:
                                    rotaxorder = rotaxorder_files[0]. \
                                                        split('.txt')[0]
                                    rotaxorder = rotaxorder. \
                                                        split('_rotaxorder')[1]
                                    rotaxorder = int(rotaxorder)
                                    mol_count = rotaxorder

                                    if verbous0:
                                        ctl.p(filename)
                                        ctl.p(rotaxorder_files)                                    
                                        ctl.p(rotaxorder)


                                zorder = 1000000-rotax_dataset[0][p]
                                pos_x = rotax_dataset[0][p]
                                pos_y = rotax_dataset[1][p]
                                segments_n = multiplicity
                                segments_filled = mol_count
                                size = rotax_dataset[2][p]*size_fact

                                colordict, color_name = \
                                        diagram.colors.group_color( \
                                                colordict, symplex_numb[sc][3])

                                diagram.marker.circle(subplots, ax, \
                                                      pos_x, pos_y, zorder, \
                                                      segments_n, \
                                                      segments_filled, \
                                                      size, color_name)


                elif mode.sort_scores == 1:
                    for rotax in plot_data_sort:
                        if rotax != []:
                            for p,symplex_dp in enumerate(rotax[0]):
                                zorder -= 1
                                pos_x = rotax[0][p]
                                pos_y = rotax[1][p]
                                segments_n = rotax[1][p] # multiplicity
                                segments_filled = rotax[4][p] # mol count
                                size = rotax[2][p]*size_fact

                                colordict, color_name = \
                                        diagram.colors.group_color( \
                                                colordict, rotax[3])

                                diagram.marker.circle(subplots, ax, \
                                                      pos_x, pos_y, zorder, \
                                                      segments_n, \
                                                      segments_filled, \
                                                      size, color_name)


                ax.set_ylabel('order of rotational symm. axis (k-fold)', \
                              fontsize=symplot_.ax_fontsize)

                multiplicities = params.get_multiplicities( \
                                            root_path, mode.microorganism_type)
                plt.yticks(multiplicities, \
                           [str(m)+'-fold' for m in multiplicities], \
                           fontsize=symplot_.ax_fontsize)

                plt.xlim(-2)
                if show_main_bd != 1:
                    plt.ylim(1.5, max(6, max(multiplicities))+0.3+ \
                            symplot_.datapoint_fontsize*0.05)
                else:
                    plt.ylim(1.5, max(6, max(multiplicities))+0.4+ \
                            symplot_.datapoint_fontsize*0.05)

                plt.grid(axis='y', color='#DDDDDD', zorder=1)

                ax.set_xticks(range(len(col_labels)))
                ax.set_xticklabels(col_labels, fontsize=symplot_.ax_fontsize)

                ax2 = ax.twiny()
                ax2.set_xticks(range(len(assembly_label)))
                ax2.tick_params(axis='x', labeltop=True, labelbottom=False)
                ax2.set_xticklabels(assembly_label, \
                                    fontsize=symplot_.ax_top_fontsize)
                ax2.set_xlim(ax.get_xlim())

                ax.tick_params(axis='x', which='both', length=0)
                ax2.tick_params(axis='x', which='both', length=0)
               
                # description of x axis
                text = symplot.x_axis_labels.x_axis_description( \
                                                            mode.show_legend)

                if mode.show_legend and symplot_.plot_width <= 10:
                    fig.text(0.17*mode.font_scale, 0.093, text, \
                             ha='right', va='top', \
                             fontsize=symplot_.ax_fontsize, zorder=20)
                    pos_h = 0.04
                    pos_v = 0.7
                    circ_rad = 0.02
                elif mode.show_legend:
                    fig.text(0.01, 0.09, text, ha='left', va='top', \
                             fontsize=symplot_.ax_fontsize, zorder=20)
                    pos_h = 0.025
                    pos_v = 0.7
                    circ_rad = 0.01

                if mode.show_legend:
                    if mode.show_intraplot_labels:
                        symplot.legend.legend_datapoints( \
                                fig, pos_h, pos_v, circ_rad, \
                                symplot_.plot_width, show_main_bd, \
                                assemb.assemblies[0][6][0], \
                                symplot_.legend_score_descr)
                    else:
                        symplot.legend.legend_datapoints_reduced( \
                                fig, pos_h, pos_v, \
                                symplot_.plot_width, \
                                mode, \
                                symplot_.legend_score_descr)


                # subchain info box
                if subchain_info != '' and len(assemb.assemblies) <= 1:
                    subchain_info_list = subchain_info.split('\n')

                    if len(subchain_info_list) > 0:
                        subchain_info_list.pop()

                    legend_txt = fig.text(0.01, \
                                    pos_v-0.2-0.08*(mode.font_scale-1), \
                                    'predictions:', \
                                    ha='left', va='top', \
                                    fontsize=symplot_.legend_fontsize, \
                                    zorder=20, \
                                    linespacing=symplot_.legend_linespacing)

                    if verbous0:
                        ctl.p(subchain_info_list)

                    for subchain_info_item in subchain_info_list:
                        chain_id_ = subchain_info_item.split( \
                                                        '\u2006·\u2006(')[0]

                        if 'FL' in chain_id_:
                            chain_id_txt = 'full length'
                        else:
                            chain_id_txt = 'subchain '+chain_id_

                        colordict, color_name = diagram.colors. \
                                        group_color(colordict, chain_id_txt)

                        legend_txt = ax.annotate( \
                            '●', xycoords=legend_txt, xy=(0, -0.3), \
                            ha='left', va='top', fontsize=12*mode.font_scale, \
                            zorder=20, \
                            linespacing=1.0, color=color_name) 
                        legend_txt_description = ax.annotate(
                            ' '+subchain_info_item, xycoords=legend_txt, \
                            xy=(1, 0.1), \
                            ha='left', va='bottom', \
                            fontsize=symplot_.legend_fontsize, zorder=20, \
                            linespacing=1.0)


                missing = False

                if subchain_mode == 1:
                    if verbous0:
                        ctl.p('output_prediction_completeness')
                        ctl.p(output_prediction_completeness)
                        ctl.p('missing predictions: ')

                    for s in assemb.assemblies:
                        for scs in s['subchain_set']:
                            if scs not in assemblies_prediction_scenario_n[ \
                                                                assembly_id]:
                                ctl.p(scs+': all')
                                missing = True
                            else:
                                missing = str(scs)+': '
                                for i,m in enumerate( \
                                    assemblies_prediction_scenario_n[ \
                                    assembly_id][scs]):

                                    if i in multiplicities and m == 0:
                                        missing += str(i)+','

                                if len(missing) > len(str(scs)+': '):
                                    ctl.p(missing[:-1])
                                    missing = True


                prefix = ''
                postfix = ''

                if subchain_mode == 0:
                    prefix = 'all'
                elif subchain_mode == 1:
                    prefix = 'subchain_set1'
                    if missing == False:
                        postfix = '_complete'

                if mode.reduced == 1 and mode.reduced_showall:
                    plt.savefig(root_path+'symplot_all'+'_'+ \
                            prefix+postfix+'.png', dpi=600)
                else:
                    if show_main_bd != 1:
                        plt.savefig(root_path+assemb.assemblies[0][1]+'_'+ \
                            prefix+postfix+'_'+ \
                            str(interface_cluster_mode)+'.png', dpi=600)
                    else:
                        plt.savefig(root_path+assemb.assemblies[0][1]+'_'+ \
                            prefix+postfix+'_'+ \
                            str(interface_cluster_mode)+'_mainbd.png', dpi=600)


                ctl.p(assemblies_prediction_scenario_n[assembly_id])
                ctl.p('subchain_sets: ')

                for s in assemb.assemblies:
                    ctl.p(s['subchain_set'])


                if (len(mode.subchain_modes) == 1 or subchain_mode == 1) and \
                   show_main_bd == 0:
                   plt.show()


                # show list of all interfaces as text output
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if mode.full_interfaces_list:
                    interfaces_txt = symplot.export.get_interfaces_txt( \
                                                    assemblies_interfaces_txt)
                    ctl.p(interfaces_txt)
