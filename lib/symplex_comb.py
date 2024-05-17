import ctl
import metadata

import glob
import os


'''
Module providing functions to determine potential combinations of SymPlexes and
calculate/assemble them.
'''

def symplex_comb(ax):
    '''
    Get combinations of SymPlexes.
    '''
    comb = []

    for i,f0 in enumerate(ax[0].files):
        for j,f1 in enumerate(ax[1].files):
            files0 = sorted(glob.glob(ax[0].path+'**/'+f0))
            files1 = sorted(glob.glob(ax[1].path+'**/'+f1))

            if len(files0) == 1 and len(files1) == 1:
                comb.append((i, j))
            else:
                ctl.p('not possible:')
                ctl.p(f0)
                ctl.p(files0)
                ctl.p(f1)
                ctl.p(files1)

    return comb


def get_symm_order(ax):
    '''
    Get order of rotational symmetry in a complex model from path.
    '''
    path = ax.model_active_path
    fn = path.split('/')[-1]

    files = sorted(glob.glob(ax.path+'**/'+fn))

    if len(files) != 1:
        return -1
    else:
        path = files[0]

    dir_name = path.split('/')[-2]
    mult = -1

    if 'symm_180' in dir_name:
        mult = 2
    elif 'symm_120' in dir_name:
        mult = 3       
    elif 'symm_090' in dir_name:
        mult = 4
    elif 'symm_072' in dir_name:
        mult = 5
    elif 'symm_060' in dir_name:
        mult = 6
    elif 'symm_051' in dir_name:
        mult = 7
    elif 'symm_045' in dir_name:
        mult = 8

    return mult


def domains(symplex_folder, conf):
    '''
    Get list of domains.
    '''
    domains = []

    sc = metadata.get_subchain_abbr(symplex_folder)

    if sc == 'FL':
        sc = (1, len(conf.domains))
    else:
        sc = metadata.subchainabbr_to_subchains(sc)

    if len(sc) == 1:
        sc = (sc[0], sc[0])

    for d in range(sc[0], sc[1]+1):
        if d not in domains:
            domains.append(d)

    return domains


def overlapping_domains(symplex0_folder, symplex1_folder, conf):
    '''
    Get list of overlapping domains.
    '''
    ov_domains = []

    symplex0_domains = domains(symplex0_folder, conf)
    symplex1_domains = domains(symplex1_folder, conf)

    for d0 in symplex0_domains:
        for d1 in symplex1_domains:
            if d0 == d1:
                if d0 not in ov_domains:
                    ov_domains.append(d0)

    return ov_domains


def seq_coverage(sc0, sc1, conf):
    '''
    Check if SymPlex combination covers full sequence.
    '''

    # check if N terminal end is included in sequence
    if not(sc0 == 'FL' or sc1 == 'FL' or \
           int(sc0[0]) == 1 or int(sc1[0]) == 1):
        return False

    # check if C terminal end is included in sequence
    if not(sc0 == 'FL' or sc1 == 'FL' or \
           int(sc0[-1]) == len(conf.domains) or \
           int(sc1[-1]) == len(conf.domains)):
        return False

    return True


def subchain_order(symplex0_folder, symplex1_folder, insertion_length, conf):
    '''
    Get upstream/downstream order of subchains.

    Return:
        int: SymPlex index with subchain that is located more upstream.
    '''
    sc0 = metadata.get_subchain_abbr(symplex0_folder)
    sc1 = metadata.get_subchain_abbr(symplex1_folder)

    if sc0 == 'FL':
        sc0 = (1, len(conf.domains))
    else:
        sc0 = metadata.subchainabbr_to_subchains(sc0)

    if sc1 == 'FL':
        sc1 = (1, len(conf.domains))
    else:
        sc1 = metadata.subchainabbr_to_subchains(sc1)

    if len(sc0) == 1:
        sc0 = (sc0[0], sc0[0])

    if len(sc1) == 1:
        sc1 = (sc1[0], sc1[0])

    if sc0[0] < sc1[0] or sc0[1] < sc1[1]:
        if insertion_length == 1:
            if sc0[1]-sc0[0] < len(conf.domains):

                return 1

        else:

            return 0

    if sc0[0] > sc1[0] or sc0[1] > sc1[1]:
        return 1

    if sc0[0] == sc1[0] or sc0[1] == sc1[1]:
        return -1

    ctl.error('subchain_order')

    return


def insertion_lengths(sc0, sc1):
    '''
    Determine insertion lengths compatible with combination of given subchains.

    Returns:
        [int]: list of numbers of domains
    '''
    insertion_lengths_ = [-1] # -1: all available domains

    if sc0 == 'FL' or sc1 == 'FL':
        insertion_lengths_.append(1)

    return insertion_lengths_


def surface_sections(symplex0_domains, symplex1_domains, sc_order, domains, \
                     start_domain, domain_start_pos, insertion_length=-1):
    '''
    Determine resulting surface sections (active domains) of symplex0 and
    symplex1.

    Params:
        symplex0_domains: domains of symplex0
        symplex1_domains: domains of symplex1
        sc_order: order of subchains
        domains: list of all domain boundaries
        start_domain: start domain
        start_domain_pos: 0:start position at the upstream side of the
                            start domain start_domain
                          1:start position at the downstream side of the
                            start domain start_domain
        insertion_length: number of domains to insert
                            -1: insertion of all available domains
    '''
    ov_domains = list(set(symplex0_domains) & set(symplex1_domains))

    # checks
    if len(ov_domains) == 0:
        ctl.error('surface_sections: len(ov_domains) == 0')

    if insertion_length not in (-1, 1):
        ctl.error('surface_sections: insertion length not supported')

    if start_domain not in ov_domains:
        ctl.error('surface_sections: alignment_domain not in ov_domains')

    if len(ov_domains) < insertion_length:
        ctl.error('surface_sections: len(ov_domains) < len(insertion_length)')


    if insertion_length == -1:
        surface_section0 = [ \
            [1, domains[start_domain-2+domain_start_pos][1]] \
            ]
        surface_section1 = [ \
            [domains[start_domain-1+domain_start_pos][0], 2000] \
            ]

    elif insertion_length == 1:
        surface_section0 = [ \
            [1, domains[start_domain-2+domain_start_pos][1]], \
            [domains[start_domain+domain_start_pos][0], 2000] \
            ]
        surface_section1 = [ \
            [domains[start_domain-1+domain_start_pos][0], \
             domains[start_domain-1+domain_start_pos][1]] \
            ]

    if sc_order == 0:

        return surface_section0, surface_section1

    elif sc_order == 1:

        return surface_section1, surface_section0

    ctl.error('surface_sections')

    return


def symplex_folders(symplex0_predscen, symplex1_predscen, \
                    subfolder_prediction_scenarios, conf):
    '''
    Get folders of given prediction scenarios sorted by the order of the
    rotational symmetry axis.
    The SymPlex with the highest order is ranked first.
    '''
    multiplicity0 = int(symplex0_predscen[-1])
    multiplicity1 = int(symplex1_predscen[-1])

    if multiplicity0 >= multiplicity1:
        symplex0_folder = subfolder_prediction_scenarios+conf.gene+'_'+ \
                                    symplex0_predscen.replace('FL', '')+'/'
        symplex1_folder = subfolder_prediction_scenarios+conf.gene+'_'+ \
                                    symplex1_predscen.replace('FL', '')+'/'
    else:
        symplex1_folder = subfolder_prediction_scenarios+conf.gene+'_'+ \
                                    symplex0_predscen.replace('FL', '')+'/'
        symplex0_folder = subfolder_prediction_scenarios+conf.gene+'_'+ \
                                    symplex1_predscen.replace('FL', '')+'/'

    return symplex0_folder, symplex1_folder


def assembly_dir_rename(assembler, assembly_return, ax, \
                        symplex_combination, \
                        insertion_length, \
                        alignment_domain, alignment_pivot_pos, \
                        conf):
    '''
    Rename assembly directory and include metadata in new directory name.
    '''
    result_infixes = { 2: 'axistilt', \
                       3: 'gap', \
                       4: 'notperpendicular', \
                       5: 'snapindist', \
                       6: 'rotsymmaxnotdetermined', \
                       7: 'rotsymmaxinverted', \
                       8: 'snapinrot', \
                       9: 'nosymmgroup', \
                      10: 'snapinrotstep', \
                      11: 'withintermini', \
                      12: 'nocoincidentres' }

    #  1:assembly possible
    #  2:assembly not possible because tilt >45°
    #  3:assembly not possible because snapin distance too large
    #  4:assembly not possible because rotational symmetry axis
    #    not perpendicular to xy plane
    #  6:assembly not possible because determination of
    #    rotational symmetry axis not successful
    #  7:assembly not possible because z component of
    #    rotational symmetry axis <=0
    #  8:assembly not possible because snapin rotation too large
    #  9:assembly not possible because no symmetry group determined
    # 10:assembly not possible because snapin rotation step too large
    # 11:assembly not possible because ax surface completely within termini
    # 12:assembly not possible because no coincident residues

    path_combination_part1 = '/'.join(conf.export_path.split('/')[:-2])+'/'

    if assembly_return == 1: # assembly possible
        validation_scores = \
                    assembler.layers[1].primitive_unit_cell.validation_scores

        os.rename(conf.export_path, path_combination_part1+ \
            str(ax[0].fold)+str(ax[1].fold)+'_'+ \
            str(metadata.get_subchain_abbr(ax[0].pathRaw))+'-'+ \
                    str(metadata.get_subchain_abbr(ax[1].pathRaw))+'_'+ \
            str(round(validation_scores[0], 3))+'_'+ \
            str(round(validation_scores[1], 3))+'_'+ \
            str(round(validation_scores[2], 3))+'_'+ \
            'd'+str(alignment_domain)+'_'+ \
            str(2*(insertion_length+1)+alignment_pivot_pos)+'_'+ \
            str(symplex_combination[0])+str(symplex_combination[1]))

    elif assembly_return in result_infixes: # assembly not possible
        os.rename(conf.export_path, path_combination_part1+ \
            str(ax[0].fold)+str(ax[1].fold)+'_'+ \
            str(metadata.get_subchain_abbr(ax[0].pathRaw))+'-'+ \
                    str(metadata.get_subchain_abbr(ax[1].pathRaw))+ \
                '_'+ \
            result_infixes[assembly_return]+'_'+ \
            'd'+str(alignment_domain)+'_'+ \
            str(2*(insertion_length+1)+alignment_pivot_pos)+'_'+ \
            str(symplex_combination[0])+str(symplex_combination[1]))

    else:
        ctl.e(assembly_return)
        ctl.error('unknown result type')

    return


def check_processed(comb_fn_str):
    '''
    Check if a combination of SymPlexes has already been calculated.
    '''
    folder = sorted(glob.glob(comb_fn_str))

    if len(folder) > 1:
        ctl.e(folder)
        ctl.error('check_processed: len(folder) > 1')
        
    if len(folder) == 1:
        return True

    return False
