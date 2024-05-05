import sys

import ctl
import metadata


'''
Module providing functions to sort coord files.
'''

def sort_by_score(prediction_files, verbous=False):
    '''
    Sort coord files descending by score.
    '''
    prediction_files = sorted(prediction_files)
    prediction_files_and_scores = []

    for f in prediction_files:
        filename = f.split('/')[-1]
        score = metadata.get_score(f)
        prediction_files_and_scores.append([f, score])

    check = check_order(prediction_files_and_scores)

    if check:

        return prediction_files

    else:
        prediction_files_and_scores = sorted(prediction_files_and_scores, \
                                             key=lambda x: x[1], reverse=True)
        check = check_order(prediction_files_and_scores)

        if check:
            if verbous:
                ctl.p('sort_by_score: sorting ok')

            prediction_files_sorted = \
                                    [i[0] for i in prediction_files_and_scores]

            return prediction_files_sorted

        else:
            ctl.e(prediction_files_and_scores)
            ctl.error('sort_by_score: ERROR')


def check_order(predictions_and_scores, verbous=False):
    '''
    Check if coord files are in descending order.
    '''
    last_desc_score = 1000

    for p in predictions_and_scores:

        if p[1] > last_desc_score:
            if verbous:
                ctl.p('check_order: not in descending order')
                ctl.p(predictions_and_scores)

            return False


        if p[1] < last_desc_score:
            last_desc_score = p[1]
        
    return True
