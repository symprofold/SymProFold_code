import ctl


'''
Module providing functions to export datasets.
'''

def get_interfaces_txt(species_interfaces):
    '''
    Export interfaces as text.
    '''
    interfaces_txt = ''
    axes_labels = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    for s in species_interfaces:
        interfaces_txt += ("\r\n"+s[0]+"\r\n")
        for i,cluster in enumerate(s[1]):
            if cluster != '': 
                interfaces_txt += ('Cluster '+axes_labels[i]+':'+"\r\n"+ \
                                        '----------'+"\r\n")
                interfaces_txt += cluster

    return interfaces_txt
