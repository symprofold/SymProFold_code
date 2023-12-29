

class PlotLayout():
    '''
    This class describes the data frame of a SymPlot layout.
    The data (points) themselves are not described by this class.
    '''


    def __init__(self, mode):
        ''' Initialization of the PlotLayout class. '''

        self.mode = mode

        self.plot_width = 10
        self.plot_height = 6

        self.ax_fontsize = None
            # label of x and y axis (bottom, left)
        self.ax_top_fontsize = None
            # label of top axis

        self.datapoint_linespacing = 1.2
        self.datapoint_fontsize = None
            # standard font size for labels of data points
        self.datapoint_fontsize_bottom = None
            # font size for bottom labels of data points

        self.datapoint_circle_size = None
            # size of circle

        self.subplots_adjust_left = None
            # parameter for subplots_adjust(left=...)
        self.subplots_adjust_right = 0.99
            # parameter for subplots_adjust(right=...)
        self.subplots_adjust_top = None
            # parameter for subplots_adjust(top=...)
        self.subplots_adjust_bottom = 0.10
            # parameter for subplots_adjust(bottom=...)

        # params for legend
        self.legend_fontsize = 8
        self.legend_linespacing = 1.4
        self.legend_score_descr = self.get_legend_score_descr()

        return


    def get_legend_score_descr(self):
        '''
        Get score description used in legend.
        '''
        if self.mode.show_nnptm == 0:
            score_descr = 'ipTM+pTM\nscore'
        elif self.mode.show_nnptm == 1:    
            score_descr = 'max(dipTM)\ndimer ipTM'

        return score_descr


    def set_subplots_adjust(self, show_legend, reduced_mode_showall):
        '''
        Determine margin parameters.
        '''
        if self.plot_width <= 10:
            self.subplots_adjust_left = 0.18
        else:
            self.subplots_adjust_left = 0.10

        if not show_legend and reduced_mode_showall:
            self.subplots_adjust_left = 0.04

        elif not show_legend:
            self.subplots_adjust_left = 0.06
            
        return


    def subplots_adjust(self, fig):
        '''
        Set margin parameters to Figure object.
        '''
        if self.subplots_adjust_left != None:
            fig.subplots_adjust(left=self.subplots_adjust_left)

        if self.subplots_adjust_right != None:
            fig.subplots_adjust(right=self.subplots_adjust_right)

        if self.subplots_adjust_top != None:
            fig.subplots_adjust(top=self.subplots_adjust_top)

        if self.subplots_adjust_bottom != None:
            fig.subplots_adjust(bottom=self.subplots_adjust_bottom)

        return fig


    def set_datapoint_fontsizes(self, column_n):
        '''
        Determine label font sizes and circle size for data points.

        Returns:
            datapoint_fontsize: standard font size for labels of
                    data points
            datapoint_fontsize_bottom: font size for bottom labels of
                    data points
            datapoint_circle_size: size of circle
        '''
        # standard params
        datapoint_fontsize = 8
        datapoint_fontsize_bottom = 8
        datapoint_circle_size = 500

        # modification of standard params e.g. for large number of columns
        # (column_n)
        if column_n > 55 and self.mode.reduced == 0:
            datapoint_fontsize = 6.5
            datapoint_fontsize_bottom = 5.0
        elif column_n > 50 and self.mode.reduced == 0:
            datapoint_fontsize = 6.5
            datapoint_fontsize_bottom = 5.8
        elif column_n > 40 and self.mode.reduced == 0:
            datapoint_fontsize = 7.0
            datapoint_fontsize_bottom = 6.0
        elif column_n > 35 and self.mode.reduced == 0:
            self.datapoint_fontsize_bottom = 6.5
        elif column_n > 30 and self.mode.reduced == 0:
            datapoint_fontsize_bottom = 7.0

        elif column_n > 60 and self.mode.reduced != 0:
            datapoint_fontsize = 6.0
            datapoint_fontsize_bottom = 4.0
            datapoint_circle_size = 100
        elif column_n > 35 and self.mode.reduced != 0:
            datapoint_fontsize = 7
            datapoint_fontsize_bottom = 5.0
            datapoint_circle_size = 100

        self.datapoint_fontsize = datapoint_fontsize
        self.datapoint_fontsize_bottom = datapoint_fontsize_bottom
        self.datapoint_circle_size = datapoint_circle_size

        return self.datapoint_fontsize, self.datapoint_fontsize_bottom, \
               self.datapoint_circle_size


    def set_ax_fontsizes(self, clusters_n):
        '''
        Determine font sizes for axis labels.

        Returns:
            ax_fontsize: label of x and y axis (bottom, left)
            ax_top_fontsize: label of top axis
        '''
        if clusters_n >= 30:
            ax_fontsize = 6.8
            ax_top_fontsize = 6.8
        elif clusters_n >= 12:
            ax_fontsize = 7.0
            ax_top_fontsize = 7.8
        elif clusters_n >= 10:
            ax_fontsize = 7.0
            ax_top_fontsize = 8.0
        elif clusters_n >= 9:
            ax_fontsize = 8.0
            ax_top_fontsize = 8.0
        elif clusters_n >= 8:
            ax_fontsize = 8.5
            ax_top_fontsize = 8.0
        else:
            ax_fontsize = 9
            ax_top_fontsize = 9

        self.ax_fontsize = ax_fontsize
        self.ax_top_fontsize = ax_top_fontsize

        return self.ax_fontsize, ax_top_fontsize
