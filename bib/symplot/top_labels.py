import ctl


'''
Module providing functions to handle assembly labels to use as top labels of a
SymPlot.
'''

def assembly_labels(assemb, cluster_n, col_w):
    '''
    Determine list of assembly labels to use as top labels of a SymPlot.
    '''
    assembly_labels = []
    col_num = 0 # each cluster is a column

    # iterate through assemblies
    for a in assemb.assemblies:
        assembly_width = 0

        for i in range(0, cluster_n):
            for j in range(0, col_w[col_num]):
                assembly_width += 1

            assembly_width += 1 # insert gap between cluster columns
            col_num += 1

        if assembly_width == 0:
            assembly_labels.append(a['assembly_label'])
        elif assembly_width == 1:
            assembly_labels.append(a['assembly_label'])
            assembly_labels.append('')
        elif assembly_width == 2:
            assembly_labels.append('')
            assembly_labels.append(s['assembly_label'])
            assembly_labels.append('')
        else:
            for j in range(0, round((assembly_width-1)/2)-1):
                assembly_labels.append('')

            assembly_labels.append(a['assembly_label'])

            for j in range(round((assembly_width-1)/2), assembly_width):
                assembly_labels.append('')

    return assembly_labels
