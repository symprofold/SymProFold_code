import ctl


'''
Module providing check functions for prediction scenarios.
'''

def check_max_models_per_prediction_scenario(assembly, predictions):
    '''
    Check maximal models per prediction scenario.
    '''
    if assembly[7] != []:

        if len(predictions) > assembly[7][0]:
            ctl.d('ERROR: check_max_models_per_prediction_scenario: '+ \
                  'max models per subchain exceeded')
            ctl.d(assembly[7])
            ctl.d(len(predictions))
            ctl.d(predictions)

    return


def check_prediction_completeness(assemb, assemblies_prediction_scenario_n):
    '''
    Check completenes of prediction scenarios.
    The desired number is 5 or a multiple of 5.
    '''
    output_prediction_completeness = ''

    for i_,a in enumerate(assemb.assemblies):
        for sc in assemblies_prediction_scenario_n[i_]:
            for i,m in enumerate(assemblies_prediction_scenario_n[i_][sc]):

                # check completeness of a prediction scenario
                if m%5 != 0:
                    if m < 20:
                        output_prediction_completeness += \
                            'subchain: '+str(sc)+', '+'mult: '+str(i)+'\n'
                        output_prediction_completeness += \
                            str(assemblies_prediction_scenario_n[i_][sc])+'\n'
                        output_prediction_completeness += \
                            str(m)+'\n'

    return output_prediction_completeness
