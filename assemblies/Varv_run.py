import os
import sys

sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import lib

# import SymProFold libraries
sys.path.append(lib.get_main_dir()+'lib/')

import symplex_comb

from structure.modelreg import ModelReg
from structure.axis import Axis
from assembler import Assembler
import chimerax_api
from config import Config


conf = Config(os.path.realpath(__file__))
try:
    conf.set_run_level(run_level)
except NameError:
    pass

sess = chimerax_api.ChimeraxSession(session)


# options for runlevel 0
conf.set_flatten_modes([0,1,2,4,5], 0)
        # 0: pure superposition, 1: flattened, 2: snapin/tile,
        # 3:completed chains, 4:primitive_unit_cell,
        # 5: assembly of 3x3 primitive unit cells
conf.set_delete_termini_modes([0,1], 0)


conf.set_species('Varv', 'Viridibacillus arvi')
conf.set_gene('A0A0K2Z0V7')
conf.import_domains()

symplex0_folder = 'A0A0K2Z0V7/A0A0K2Z0V7_28x4/'
model_status0 = 1 # set oriented path

symplex1_folder = 'A0A0K2Z0V7/A0A0K2Z0V7_13x4/'
model_status1 = 1 # set oriented path

alignment_pivot_pos = [0, 1][0]
insertion_length = [-1, 1][0]


symplex0_domains = symplex_comb.domains(symplex0_folder, conf)
symplex1_domains = symplex_comb.domains(symplex1_folder, conf)
sc_order = symplex_comb.subchain_order( \
                    symplex0_folder, symplex1_folder, insertion_length, conf)

surface_section0, surface_section1 = \
    symplex_comb.surface_sections(symplex0_domains, symplex1_domains, \
                                  sc_order, conf.domains, \
                                  2, alignment_pivot_pos, \
                                  insertion_length)


for conformation in  conf.conformations:
    model_reg = ModelReg()
    sess.set_model_reg(model_reg)

    ax = [Axis(model_reg), Axis(model_reg)]

    ax[0].set_session(sess, conf)
    ax[0].set_folder(symplex0_folder, model_status0)
    ax[0].set_surface(surface_section0)

    ax[1].set_session(sess, conf)
    ax[1].set_folder(symplex1_folder, model_status1)
    ax[1].set_surface(surface_section1)

    # ax[0].set_model_active(0)
    # ax[1].set_model_active(0)

    assembler = Assembler(ax, conf)
    assembler.alignment_pivot_pos = 0
    assembler.run()
