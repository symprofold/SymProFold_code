import glob

import ctl
import metadata


'''
Module providing functions to filter oligomer model.
'''

def filter_rotang(files, multiplicities, microorganism_type=0, verbous=False):
    '''
    Filter a list of oligomer models by rotational symmetry angles within a
    tolerance.
    '''
    files_filtered = []

    for f in files:
        rotang = metadata.get_rotsymm_ang(f)

        if len(rotang) != 0:
            if microorganism_type == 0 and 5 not in multiplicities:
                    # microorganism_type: 0: bacteria, 1: virus
                t = 12 # tolerance

                if min(rotang) < 60-t or \
                   max(rotang) < 60-t or \
                   60+t < min(rotang) < 90-t or \
                   60+t < max(rotang) < 90-t or \
                   90+t < min(rotang) < 120-t or \
                   90+t < max(rotang) < 120-t or \
                   120+t < min(rotang) < 180-t or \
                   120+t < max(rotang) < 180-t or \
                   180+t < min(rotang) or \
                   180+t < max(rotang):

                    if verbous:
                        ctl.p(rotang)
                        ctl.p('filter_rotang: angle not in tolerance')

                else:
                    files_filtered.append(f)

            elif microorganism_type == 0 and 5 in multiplicities:
                t = 12 # tolerance

                if min(rotang) < 60-t or \
                   max(rotang) < 60-t or \
                   60+t < min(rotang) < 72-t or \
                   60+t < max(rotang) < 72-t or \
                   72+t < min(rotang) < 90-t or \
                   72+t < max(rotang) < 90-t or \
                   90+t < min(rotang) < 120-t or \
                   90+t < max(rotang) < 120-t or \
                   120+t < min(rotang) < 180-t or \
                   120+t < max(rotang) < 180-t or \
                   180+t < min(rotang) or \
                   180+t < max(rotang):

                    if verbous:
                        ctl.p(rotang)
                        ctl.p('filter_rotang: angle not in tolerance')

                else:
                    files_filtered.append(f)

            elif microorganism_type == 1:
                t = 5 # tolerance

                if min(rotang) < 60-t or \
                   max(rotang) < 60-t or \
                   60+t < min(rotang) < 72-t or \
                   60+t < max(rotang) < 72-t or \
                   72+t < min(rotang) < 90-t or \
                   72+t < max(rotang) < 90-t or \
                   90+t < min(rotang) < 120-t or \
                   90+t < max(rotang) < 120-t or \
                   120+t < min(rotang) < 180-t or \
                   120+t < max(rotang) < 180-t or \
                   180+t < min(rotang) or \
                   180+t < max(rotang):

                    if verbous:
                        ctl.p(rotang)
                        ctl.p('filter_rotang: angle not in tolerance')

                else:
                    files_filtered.append(f)

    return files_filtered


def filter_clashes(files, verbous=False):
    '''
    Filter a list of coord files by clashes.
    '''
    files_filtered = []

    for f in files:
        cla = clashes_converted(f)
        roll_cla = roll_clashes_converted(f)
        
        if cla <= 3 and roll_cla <= 6:
            files_filtered.append(f)

        if cla > 30:
            if verbous:
                ctl.p('filter_clashes: many clashes')
                ctl.p(cla)
                ctl.p(f)

    return files_filtered


def clashes_converted(file_path):
    '''
    Get clashes per 100aa (converted) of a complex model from path.
    '''
    filename = file_path.split('/')[-1]

    if not '_cla' in filename:
        ctl.e(file_path)
        ctl.error('clashes_converted: clashes metadata not found in path')

    elif 'unrelaxed_' in filename:
        cla = metadata.get_clashes(file_path)
        cla = round(cla*0.05, 2)
            # compensate for the higher number of clashes in unrelaxed
            # coord files using a factor of 20 between unrelaxed and relaxed

    else:
        cla = metadata.get_clashes(filename)

    return cla


def roll_clashes_converted(file_path):
    '''
    Get clashes per 100aa (rolling range of 200aa, converted) of a
    complex model from path.
    '''
    filename = file_path.split('/')[-1]

    if not '_cla' in filename:
        ctl.e(file_path)
        ctl.error('roll_clashes_converted: '+ \
                            'roll clashes metadata not found in path')

    elif 'unrelaxed_' in filename:
        cla = metadata.get_roll_clashes(file_path)
        cla = round(cla*0.05, 2)
            # compensate for the higher number of clashes in unrelaxed
            # coord files using a factor of 20 between unrelaxed and relaxed

    else:
        cla = metadata.get_roll_clashes(file_path)

    return cla


def filter_intermol_betasheet(files, verbous=False):
    '''
    Filter a list of coord files by the property that parts of the protein
    chain run erroneously through the sphere of another monomer.

    This filter evaluates the fraction of residues involved in medium-length
    intermolecular beta sheet strands.
    '''
    files_filtered = []
    filtered_out = []

    for f in files:
        fraction_minlen, fraction, betasheet_res_n = \
                                        metadata.get_intermol_betasheets(f)
        mult = int(metadata.get_prediction_scenario(f).split('x')[1])
        res_n = round(betasheet_res_n/mult)

        if not (fraction_minlen > 0.08 and fraction > 0.12 and res_n > 75) \
           or fraction >= 0.5:
            files_filtered.append(f)
        else:
            filtered_out.append(f)

    return files_filtered, filtered_out
