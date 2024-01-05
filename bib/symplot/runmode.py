

class RunMode():
    '''
    This class describes the run mode of the script.
    '''


    def __init__(self):
        ''' Initialization of the Mode class. '''

        # reduced mode
        # ~~~~~~~~~~~~
        self.reduced = 0 # 0:no, 1:yes
        self.reduced_part = 0 # 1:yes,part1, 2:yes,part2
            # reduced mode: only one prediction per prediction scenario
        self.reduced_showall = [False, True][0]
        self.reduced_cluster_n = 2 # maximal number of clusters in reduced mode


        # text output of interfaces
        # -------------------------
        self.full_interfaces_list = [False, True][0]


        # clustering
        # ----------
        self.clustering = 1
            # 0:no clustering, 1:clustering

        self.interface_matrix_type = 1
            # 0:unsigned matrix, 1:signed matrix

        self.interface_cluster_modes = range(0, 1)
            # 0: set edge weights to correlation coeff
            # 1: set all edge weights to 1


        # main bd
        # -------
        self.show_main_bd_modes = [0, 1]
            # show main bd in plot, hide score and subchains
        self.plot_height_main_bd = 5.7
            # plot height for show_main_bd_modes == 1


        # other
        # -----
        self.subchain_modes = [ [1, 0], [0] ][0]
            # 0:all, 1:standard subchain set

        self.microorganism_type = None # not set
            # 0: bacteria, 1: virus

        self.filter_clashes = 1
            # 0: no, 1: yes

        self.show_clash_score = 2
            # 0:no, 1:yes, 2:indicator symbol

        self.show_nnptm = 0
            # 0: no, 1: yes

        self.filter_angles = 1 # stricter range for angles
            # 0: no, 1: yes

        self.filter_angles = 1 # stricter range for angles
            # 0: no, 1: yes

        self.show_intraplot_labels = [True, False][0]

        self.sort_scores = 0
            # 0: - sort plot data by rot_ax_ind and subchain_name
            #    - data point coloring: one color per subchain
            # 1: - sort plot data by rot_ax_ind
            #    - data point coloring: one color for all subchains

        self.circ_calibration_fact = 0.031
            # calibration factor for circle size if
            # show_intraplot_labels==False, value is system-dependent

        self.show_legend = [True, False][0]

        self.plotlayout = 0
            # 0: plot layout for print
            # 1: plot layout for presentations

        self.font_scale = 1

        self.plottype = None
            # 0:specific species
            # 1:general plot with many species

        self.verbous = False # verbous log

        return


    def update(self):
        '''
        Update dependencies.
        '''
        if self.full_interfaces_list:
            self.reduced_showall = True

        if self.plottype == 1:
            self.sort_scores = 1

            if self.full_interfaces_list:
                self.reduced = 0
                self.show_main_bd_modes = [1]
            else:
                self.reduced = 1

            if self.reduced == 1:
                self.show_legend = False
                self.show_main_bd_modes = [0]

        if self.plotlayout == 1:
            self.font_scale = 1.5

        return
