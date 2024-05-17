import glob
import os
import sys

sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import lib


# Parameter (can affect the separation results)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
params = {
    'min_contact_area_a': 0.1,
    'min_contact_area_b': 0.0015,
        # parameter to calculate minimum fraction of contact area for merging
        # of domain subsections
        # formula: area * ( a + b*dom_sect_len )
        #
        # Domain subsections with at least minimum fraction can be
        # merged, merging partner can be larger and have a smaller fraction.

    'maxdist_betasheet': 200,
        # parameter for crosslinking of beta sheets:
        # Between between strands of a beta sheet, external chain sections up
        # to a maximum length of maxdist_betasheet are included.

    'mincontactarea_alphahelix_fract': 0.05,
        # parameter for merging domain subsections containing alpha helices:
        # Fraction of inter domain subsection contact area of alpha helices
        # in one domain subsection.
        # Domain subsections with at least minimum fraction can be
        # merged, merging partner can be larger and have a smaller fraction.

    'mincontactarea_alphahelix': 200,
        # parameter for merging domain subsections containing alpha helices:
        # Minimal subsection contact area of alpha helices in one domain
        # subsection.
        # Domain subsections with at least minimum fraction can be
        # merged, merging partner can be larger and have a smaller fraction.

    'min_output_protparts': 1,
        # minimum number of output protein parts (protein subsections) for
        # beta sheet containing proteins that are not separatable by single
        # chain cuts
}


# General options (not affecting the result)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
options = {
    'show_plot': [False,True][0],
        # show plots

    'dom_coloring': [False,True][1],
        # color the domains in different colors

    'show_steps': [False,True][1],
        # show steps in ChimeraX

    'export_cxs': [False,True][0],
        # export ChimeraX file with graphical representation of the
        # individual steps

    'verbous': [False,True][1]
}


# import Domain_Separator libraries
path_symprofold = lib.get_main_dir()
path_domain_separator = path_symprofold+'tools/domain_separator/'
sys.path.append(path_domain_separator+'bib/')

import bibpdb
import ctl
import filesystem

work_dir = os.path.dirname(os.path.realpath(__file__))+'/'
monomer_files_relaxed = sorted(glob.glob(work_dir+ \
                                        '*_x1/*_x1/rank00_0.*.pdb'))
monomer_files_unrelaxed = sorted(glob.glob(work_dir+ \
                                        '*_x1/unrelaxed_rank00_0.*.pdb'))
monomer_files_ranked_relaxed = sorted(glob.glob(work_dir+ \
                                        '*_x1/*_x1/ranked_0.pdb'))
monomer_files_ranked_unrelaxed = sorted(glob.glob(work_dir+ \
                                        '*_x1/ranked_0.pdb'))
monomer_file = None

if len(monomer_files_relaxed) >= 1:
    monomer_file = monomer_files_relaxed[0]
elif len(monomer_files_unrelaxed) >= 1:
    monomer_file = monomer_files_unrelaxed[0]
elif len(monomer_files_ranked_relaxed) >= 1:
    monomer_file = monomer_files_ranked_relaxed[0]
elif len(monomer_files_ranked_unrelaxed) >= 1:
    monomer_file = monomer_files_ranked_unrelaxed[0]
else:
    print('ERROR: no top-ranked model found')
    sys.exit()
    
import chimerax_api
import separate_domains

from modelreg import ModelReg

input_files = [monomer_file]
sess = chimerax_api.ChimeraxSession(session, 1)

accession_code = monomer_file.split('/')[-2].split('_')[0]
f_out_prefix = work_dir+accession_code

for f in input_files:
    ctl.d('starting file '+f)
    bibpdb.clean_pdb(f, f+'_cleaned.pdb')
    separate_domains.separate_domains(f+'_cleaned.pdb', f_out_prefix, \
                                      params, options, sess)
    os.unlink(f+'_cleaned.pdb')
