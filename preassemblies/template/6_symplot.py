import os
import pathlib
import sys

sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import lib

# import SymProFold libraries
path_symprofold = lib.get_main_dir()
sys.path.append(path_symprofold+'lib/')


import ctl
import symplot.bib
import symplot.config
import symplot.symplot

from symplot.assemblies import Assemblies
from symplot.layout import PlotLayout
from symplot.runmode import RunMode


if __name__ == "__main__":
    root_path = str(pathlib.Path().absolute())+'/'
    mode = RunMode()


    assemblies0 = symplot.config.assembly_config()
    # list for each assembly
    #       0:species name, 1:gene id, 2:symm group
    #       3:main folder of predictions
    #       4:include in combined plot [0:no, 1:yes,part1, 2:yes,part2]
    #       5:large plot for full ensemble w/o interface clustering
    #               [0:no, 1:yes]
    #       6:filter params [lowest ranking score]
    #       7:check params [max number of models per prediction scenario]
    #       8:display params [number of clusters in reduced mode,
    #               [ranks of clusters to exclude]]
    #       9:clustering_mode, 0:no clustering, 1:clustering


    # primary settings
    # ================

    # reduced mode
    # ------------
    mode.reduced = 0 # 0:no, 1:yes

    mode.reduced_part = [1, 2][0] # 1:yes,part1, 2:yes,part2
        # reduced mode: only one prediction per prediction scenario

    mode.reduced_showall = [False, True][0]
    mode.reduced_cluster_n = 2 # maximal number of clusters in reduced mode


    # text output of interfaces
    # -------------------------
    mode.full_interfaces_list = [False, True][0]


    # clustering
    # ----------
    mode.clustering = 1
        # 0:no clustering, 1:clustering

    mode.interface_matrix_type = 1
        # 0:unsigned matrix, 1:signed matrix

    mode.interface_cluster_modes = range(0, 1)
        # 0: set edge weights to correlation coeff
        # 1: set all edge weights to 1


    # secondary settings
    # ==================
    mode.filter_clashes = 1
        # 0: no, 1: yes

    mode.filter_angles = 1 # stricter range for angles
        # 0: no, 1: yes

    mode.filter_intermol_chain = True
        # [False, True]
        # filter by the property that parts of the protein chain run
        # erroneously through the sphere of a another monomer

    mode.microorganism_type = None # not set
        # 0: bacteria, 1: virus

    mode.show_intraplot_labels = [True, False][0]

    mode.sort_scores = 0
        # 0: - sort plot data by rot_ax_ind and subchain_name
        #    - data point coloring: one color per subchain
        # 1: - sort plot data by rot_ax_ind
        #    - data point coloring: one color for all subchains

    mode.circ_calibration_fact = 0.031
        # calibration factor for circle size if show_intraplot_labels==False,
        # dependend on system

    mode.show_legend = [True, False][0]

    mode.verbous = False # verbous log


    # update dependencies between different modes
    # ===========================================
    mode.update()
    

    # initialization
    # ==============

    # layout params for plot
    # ----------------------
    symplot_ = PlotLayout(mode)


    # initialization
    # --------------

    # get type of microorganism (bacteria, virus)
    if mode.microorganism_type == None:
        mode.microorganism_type = symplot.bib.get_microorganism_type(root_path)

    mode.plottype = symplot.bib.get_plottype(root_path)

    assemb = Assemblies(assemblies0)
        # data structure for all assemblies

    # Reduce array dependend on type of plot.
    # Types of plot: specific species, general plot with many species
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    assemb.reduce_to_plottype(root_path, mode)


    if mode.plottype == 1:
        # plottype: 0:specific species, 1: general plot with many species

        symplot_.plot_width = 10
        symplot_.plot_height = 5
        symplot_.datapoint_linespacing = 1.0

        if mode.reduced_showall:
            symplot_.plot_width = 18
            symplot_.plot_height = 5

        if assemb.assemblies[0][5] == 1:
                # large plot for full ensemble w/o rot axis clustering
                # 0:no, 1:yes
            symplot_.plot_width = 18


    mode.update()


    # determine margin parameters for SymPlot
    symplot_.set_subplots_adjust(mode.show_legend, mode.reduced_showall)

    # add domain info for each species
    assemb.add_domain_info()

    # add standard set of subchains
    assemb.add_subchains()

    # set labels for each assembly
    assemb.set_assembly_labels(mode.reduced)


    symplot.symplot.calc_symplot(root_path, mode, assemb, symplot_)

    print('finished')
