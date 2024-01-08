import os
import sys
sys.path.append(os.path.dirname(__file__)+'/../bib/')

from structure.modelreg import ModelReg
from structure.axis import Axis
from assembler import Assembler
import chimerax_api
from config import Config


conf = Config()
try:
    conf.set_run_level(run_level)
except NameError:
    pass

sess = chimerax_api.ChimeraxSession(session)


#options for runlevel 0
conf.set_delete_termini_modes([1], 0)
conf.set_flatten_modes([0], 0) # not used for Cglu
        # 0: pure superposition, 1: flattened, 2: snapin/tile,
        # 3:completed chains, 4:primitive_unit_cell,
        # 5: assembly of 3x3 primitive unit cells
conf.set_filters_for_export([0,1,2], 0)
conf.set_snapshot_modes([0], 0)   # snapshot_w_separated_chains

#options for runlevel 5
conf.set_flatten_modes([0], 5) # not used for Cglu
conf.set_filters_for_export([0,1,2], 5)
conf.set_snapshot_modes([0], 5)   # snapshot_w_separated_chains


conf.set_species('Cglu', 'Corynebacterium glutamicum')
conf.set_gene('Q6QUS5')
conf.import_domains()


for conformation in  conf.conformations:
    model_reg = ModelReg()
    sess.set_model_reg(model_reg)

    ax = [Axis(model_reg), Axis(model_reg)]

    ax[0].set_session(sess, conf)
    ax[0].set_folder('Q6QUS5_x6_superposed/', 1) # set oriented path

    ax[1].set_session(sess, conf)
    ax[1].set_folder('Q6QUS5_x3_o/symm_120/', 0)

    # ax[0].set_model_active(0)
    # ax[1].set_model_active(0)

    assembler = Assembler(ax, conf)
    assembler.run()
