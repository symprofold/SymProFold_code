import ctl


'''
Module providing functions to handle labels and description of the x axis of a
SymPlot.
'''

def x_axis_description(show_legend):
    '''
    Get description of x axis.
    '''
    text = 'rotational symm. axis'+'\n'

    if not show_legend:
        text = ''

    return text


def x_axis_labels_interfacecluster(assemb, col_per_assembly):
    '''
    Generate x axis (bottom) labels for each cluster.
    '''
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    col_labels = []

    for i,a in enumerate(assemb.assemblies):
        col_labels_assembly = ['axis '+alphabet[i]+'\n' \
                                        for i in range(0, col_per_assembly)]
        col_labels += col_labels_assembly

    return col_labels


def x_axis_labels_w_spacher_interfacecluster(assemb, col_labels0, \
                                             data_points_per_col):
    '''
    Determine x axis (bottom) labels with spacer.

    Args:
        assemb: Assemblies object
        col_labels0: column labels w/o spacer
        data_points_per_col: layout precision, number of datapoints per
                fig column
    '''
    # get (maximal) width of each column
    col_widths = assemb.get_cluster_col_widths()

    col_labels = [] # column labels with spacer
    col_w = [] # column widths w/o spacer
    col_centerpos = [0, ] # reference positions for centering

    for i,col_width in enumerate(col_widths):
        if col_width < 3:
            col_width = 3 # set minimum width for displayed cluster column

        col_w.append(col_width)
        centerpos = round((col_width-round(data_points_per_col/2))/ \
                                        data_points_per_col)
        col_centerpos.append(col_centerpos[-1]+centerpos)

        if len(col_labels0) > i:
            if centerpos == 0:
                col_labels.append(col_labels0[i])
            elif centerpos == 1:
                col_labels.append(col_labels0[i])
                col_labels.append('')
            elif centerpos == 2:
                col_labels.append('')
                col_labels.append(col_labels0[i])
                col_labels.append('')
            else:
                for j in range(0, round((centerpos-1)/2)):
                    col_labels.append('')

                col_labels.append(col_labels0[i])

                for j in range(round((centerpos-1)/2), centerpos):
                    col_labels.append('')

    return col_labels, col_w, col_centerpos
