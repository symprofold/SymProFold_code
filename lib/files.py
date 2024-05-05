import ctl

import symplot.prediction_scenario


'''
Module providing functions to handle and manage files, e.g. predicted
coord files.
'''

def get_model_status(conf, model_folder):
    '''
    Get status of predicted model.

    Return:
        int: 1:oriented, unrelaxed
             2:oriented, relaxed
    '''
    path = conf.path_ax_predictions+conf.species+'/'+model_folder[:-1]
    folders = symplot.prediction_scenario.get_folders(path.replace('FL', ''))

    if len(folders) != 1:
        ctl.e(path)
        ctl.e(folders)
        ctl.error('get_model_status: folder not found')

    if folders[0][-2:] == '_o':
        model_status = 1
    elif folders[0][-3:] == '_or':
        model_status = 2
    else:
        ctl.error('get_model_status: model status unclear')

    return model_status
