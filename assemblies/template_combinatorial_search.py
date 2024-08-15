import os
import sys

sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import lib

# import SymProFold libraries
sys.path.append(lib.get_main_dir()+'lib/')

import ctl
import chimerax_api
import files
import metadata
import symplex_comb

import glob
import random

from structure.modelreg import ModelReg
from structure.axis import Axis
from assembler import Assembler
from config import Config


conf = Config(os.path.realpath(__file__))
try:
    conf.set_run_level(run_level)
except NameError:
    pass

sess = chimerax_api.ChimeraxSession(session)


# Configuration of the combinatorial search for a specific assembly
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
conf.set_species('Vaer', 'Vibrio aerogenes')
        # species abbreviation, species name (replace example entries)

conf.set_gene('A0A1M5ZCF8')
        # gene name (replace example entry)

conf.set_symplex_path('')
        # subdirectory in assembly directory where the SymPlex
        # (Symmetric Protein Complex) models are stored (replace example entry)

symplex0_predscenarios = ['23x2', '23x4', '23x6', '3x2', '3x4']
symplex1_predscenarios = ['12x4', '12x6', '1x3', '1x4', '1x6', 'FLx4']
        # prediction scenarios to be included in the combinatorial search
        # (replace example entries)

mode = 0 # 0:random search order, 1: sorted search order

# end of configuration section


# options for runlevel 0
conf.set_flatten_modes([5], 0) # [0,1,2,4,5]
        # 0:pure superposition, 1:flattened, 2:snapin/tile,
        # 3:completed chains, 4:primitive_unit_cell,
        # 5:assembly of 3x3 primitive unit cells
conf.set_snapshot_modes([0], 0) # snapshot_w_separated_chains
conf.set_delete_termini_modes([1], 0)


path_export_postfix_running = '_running'
            # postfix for temporary name of export folder during runtime

conf.set_export_path_postfix('_'+conf.gene+'_'+conf.version+'/'+ \
                                        conf.species+'_'+conf.species_name)
conf.import_domains()

if os.path.exists(conf.path_ax_predictions+conf.species+'/'+ \
                                    conf.symplex_path) == False:
    ctl.error('path_prediction_scenarios not found')


alignment_pivot_pos = [0, 1]
    # 0: pivot residue at the upstream side of the alignment section
    #    (e.g. domain)
    # 1: pivot residue at the downstream side of the alignment section
    #    (e.g. domain)

predscen_symplex_combinations = symplex_comb.predscen_symplex_comb( \
                    symplex0_predscenarios, symplex1_predscenarios, sess, conf)

if mode == 0:
    random.shuffle(predscen_symplex_combinations)

for predscen_symplex_combination in predscen_symplex_combinations:
    symplex0_predscen = predscen_symplex_combination[0][0]
    symplex1_predscen = predscen_symplex_combination[0][1]
    symplex0_folder, symplex1_folder = symplex_comb. \
                    symplex_folders(symplex0_predscen, symplex1_predscen, \
                                    conf.symplex_path, conf)

    ctl.p('symplex0_folder:')
    ctl.p(symplex0_folder)
    ctl.p('symplex1_folder:')
    ctl.p(symplex1_folder)

    model_reg = ModelReg()
    sess.set_model_reg(model_reg)

    ax = [Axis(model_reg), Axis(model_reg)]

    ax[0].set_session(sess, conf)
    model_status0 = files.get_model_status(conf, symplex0_folder)
    ax[0].set_folder(symplex0_folder, model_status0)

    ax[1].set_session(sess, conf)
    model_status1 = files.get_model_status(conf, symplex1_folder)
    ax[1].set_folder(symplex1_folder, model_status1)

    sc0 = metadata.get_subchain_abbr(ax[0].pathRaw)
    sc1 = metadata.get_subchain_abbr(ax[1].pathRaw)

    # check if SymPlex combination covers full sequence
    if symplex_comb.seq_coverage(sc0, sc1, conf) == False:
        ctl.e(sc0)
        ctl.e(sc1)
        ctl.error('SymPlex combination does not cover the full sequence')

    c = predscen_symplex_combination[1]
            # combination of SymPlexes

    ctl.p('combination:')
    ctl.p(c)

    ov_domains = symplex_comb.overlapping_domains( \
                                    symplex0_folder, symplex1_folder, conf)

    symplex0_domains = symplex_comb.domains(symplex0_folder, conf)
    symplex1_domains = symplex_comb.domains(symplex1_folder, conf)
    ctl.p(symplex0_domains)
    ctl.p(symplex1_domains)

    insertion_lengths_ = symplex_comb.insertion_lengths(sc0, sc1)
            # insertion_length: number of domains to insert
            # -1: insertion of all available domains

    for insertion_length in insertion_lengths_:
        sc_order = symplex_comb.subchain_order( \
                            symplex0_folder, symplex1_folder, \
                            insertion_length, conf)
        ctl.p('sc_order:')
        ctl.p(sc_order)

        if sc_order == -1:
            continue

        for alignment_domain in ov_domains:
            for p in alignment_pivot_pos:
                if alignment_domain == 1 and p == 0:
                    continue

                if alignment_domain == len(conf.domains) and p == 1:
                    continue

                if alignment_domain == 1 and insertion_length == 1:
                    continue

                if alignment_domain+p >= len(conf.domains) and \
                                                insertion_length == 1:
                    continue


                surface_section0, surface_section1 = \
                        symplex_comb.surface_sections( \
                                symplex0_domains, symplex1_domains, \
                                sc_order, conf.domains, \
                                alignment_domain, p, \
                                insertion_length)

                model_reg = ModelReg()
                sess.set_model_reg(model_reg)
                ax = [Axis(model_reg), Axis(model_reg)]


                # configure ax0
                # ~~~~~~~~~~~~~
                ax[0].set_session(sess, conf)
                ax[0].set_folder(symplex0_folder, model_status0)
                ax[0].set_domains([ conf.domains[alignment_domain-1] ])
                ax[0].set_surface(surface_section0)
                ax[0].set_model_active(c[0])
                symm_order0 = symplex_comb.get_symm_order(ax[0])

                if symm_order0 != ax[0].fold:
                    continue


                # configure ax1
                # ~~~~~~~~~~~~~
                ax[1].set_session(sess, conf)
                ax[1].set_folder(symplex1_folder, model_status1)
                ax[1].set_domains([ conf.domains[alignment_domain-1] ])
                ax[1].set_surface(surface_section1)
                ax[1].set_model_active(c[1])
                symm_order1 = symplex_comb.get_symm_order(ax[1])

                if symm_order1 != ax[1].fold:
                    continue


                # skip cases in which insertions of 1 domain in length do
                # not have an interface with each other
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if insertion_length == 1:
                    if sc_order == 0:
                        distogram = symplex_comb.interface_matrix_domain( \
                                            ax[1], alignment_domain, conf)
                    elif sc_order == 1:
                        distogram = symplex_comb.interface_matrix_domain( \
                                            ax[0], alignment_domain, conf)

                    if len(distogram) == 0:
                        continue


                path_combination_part1 = \
                        '/'.join(conf.export_path.split('/')[:-2])+'/'

                comb_fn_str = path_combination_part1+ \
                    str(ax[0].fold)+str(ax[1].fold)+'_'+ \
                    str(metadata.get_subchain_abbr(ax[0].pathRaw))+'-'+ \
                    str(metadata.get_subchain_abbr(ax[1].pathRaw))+'_*_'+ \
                    'd'+str(alignment_domain)+'_'+ \
                    str(2*(insertion_length+1)+p)+'_'+ \
                    str(c[0])+'-'+str(c[1])

                comb_fn_str_running = path_combination_part1+ \
                    '*'+ \
                    str(ax[0].fold)+str(ax[1].fold)+'_'+ \
                    str(metadata.get_subchain_abbr(ax[0].pathRaw))+'-'+ \
                    str(metadata.get_subchain_abbr(ax[1].pathRaw))+'_'+ \
                    'd'+str(alignment_domain)+'_'+ \
                    str(2*(insertion_length+1)+p)+'_'+ \
                    str(c[0])+'-'+str(c[1])+ \
                    path_export_postfix_running

                if symplex_comb.check_processed(comb_fn_str) == True or \
                   symplex_comb.check_processed(comb_fn_str_running) == True:
                    continue


                conf.set_export_path_postfix( \
                    '_'+conf.gene+'_'+conf.version+'/'+ \
                    conf.species+'_'+conf.species_name+'_'+ \
                    str(ax[0].fold)+str(ax[1].fold)+'_'+ \
                    str(metadata.get_subchain_abbr(ax[0].pathRaw))+'-'+ \
                    str(metadata.get_subchain_abbr(ax[1].pathRaw))+'_'+ \
                    'd'+str(alignment_domain)+'_'+ \
                    str(2*(insertion_length+1)+p)+'_'+ \
                    str(c[0])+'-'+str(c[1])+ \
                    path_export_postfix_running)
                conf.update_export_paths()

                assembler = Assembler(ax, conf)
                assembler.alignment_pivot_pos = p
                ctl.p('alignment_domain')
                ctl.p(alignment_domain)
                ctl.p(p)
                ctl.p(insertion_length)


                assembly_successful, assembly_return, result_message = \
                                                                assembler.run()

                # create reduced assembly model without hydrogens
                # to save storage space
                assembly_models = sorted(glob.glob( \
                                conf.export_path+'assembly/'+'*.cif'))

                if len(assembly_models) == 1:
                    assembly_model = assembly_models[0]

                    sess.init()
                    sess.open_model(assembly_model)
                    sess.run('delete H')
                    sess.save_model_id((1,), \
                            assembly_model.replace('.cif', '_woH.cif'))
                    os.unlink(assembly_model)

                symplex_comb.assembly_dir_rename(\
                                    assembler, assembly_return, ax, \
                                    c, insertion_length, \
                                    alignment_domain, p, conf)
