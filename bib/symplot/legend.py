from matplotlib.patches import Ellipse


'''
Module providing functions for creating legends for SymPlots.
'''

def legend_datapoints(fig, pos_h, pos_v, circ_rad, plot_width, \
                      show_main_bd, minscore, score_descr):
    '''
    Create diagram legend with data point description.
    '''
    if show_main_bd != 1:
        fig.text(pos_h, pos_v+0.18, 'first domain\nof subchain', \
                ha='center', va='top', fontsize=8, zorder=20)
        fig.text(pos_h, pos_v+0.13, 'last domain\nof subchain', \
                ha='center', va='top', fontsize=8, zorder=20)
        fig.text(pos_h, pos_v+0.08, '·n-mer', \
                ha='center', va='top', fontsize=8, zorder=20)

    if show_main_bd == 1:
        mainbd_txt = 'main binding\ndomain\u2006·\u2006domain'
        fig.text(pos_h-0.010, pos_v+0.22, mainbd_txt, \
                ha='left', va='top', fontsize=8, zorder=20, rotation=90)

    if show_main_bd != 1:
        if minscore == 0:
            minscore_txt = '0.20'
        else:
            minscore_txt = str(minscore)

        fig.text(pos_h, pos_v+0.0, \
            score_descr+'\n(if ≥ '+minscore_txt+')', \
            ha='center', va='top', fontsize=8, zorder=20)

    circ = Ellipse((pos_h, pos_v+0.03), circ_rad, circ_rad*(plot_width/6), \
                   color='orange', alpha=0.8)
    fig.add_artist(circ)

    return


def legend_datapoints_reduced(fig, pos_h, pos_v, plot_width, \
                              circ_calibration_fact, score_descr):
    '''
    Create diagram legend with reduced data point description.
    '''
    fig.text(pos_h, pos_v+0.18, \
        score_descr, \
        ha='center', va='top', fontsize=8, zorder=20)
    fig.text(pos_h, pos_v+0.106, \
        '0.2\n\n\n\n0.5\n\n\n\n\n0.8', \
        ha='center', va='top', fontsize=8, zorder=20)

    circ_rad = circ_calibration_fact*0.2
    circ0 = Ellipse((pos_h, pos_v+0.120), circ_rad, \
                    circ_rad*(plot_width/6), \
                    color='orange', alpha=0.8, linewidth=0)

    circ_rad = circ_calibration_fact*0.5
    circ1 = Ellipse((pos_h, pos_v+0.045), circ_rad, \
                    circ_rad*(plot_width/6), \
                    color='orange', alpha=0.8, linewidth=0)

    circ_rad = circ_calibration_fact*0.8
    circ2 = Ellipse((pos_h, pos_v-0.050), circ_rad, \
                    circ_rad*(plot_width/6), \
                    color='orange', alpha=0.8, linewidth=0)

    fig.add_artist(circ0)
    fig.add_artist(circ1)
    fig.add_artist(circ2)

    return
