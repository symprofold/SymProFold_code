import os
import sys
sys.path.append(os.path.dirname(__file__)+'/../bib/')

import ctl
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


# options for runlevel 0
conf.set_flatten_modes([0,1,2,3,4], 0)
        # 0: pure superposition, 1: flattened, 2: snapin/tile,
        # 3:completed chains, 4:primitive_unit_cell
conf.set_delete_termini_modes([0,1], 0)


conf.set_species('Vqui', 'Vibrio quintilis')
conf.set_gene('A0A1M7YYM3')
conf.import_domains()


for conformation in conf.conformations:
    model_reg = ModelReg()
    sess.set_model_reg(model_reg)

    ax = [Axis(model_reg), Axis(model_reg)]

    ax[0].set_session(sess, conf)
    ax[0].set_folder('A0A1M7YYM3_12x4/', 1) #set oriented path
    ax[0].set_surface([ [1,509] ])

    ax[1].set_session(sess, conf)
    ax[1].set_folder('A0A1M7YYM3_23x4/', 2) #set oriented path
    ax[1].set_surface([ [510,1000] ])

    # ax[0].set_model_active(0)
    # ax[1].set_model_active(1)

    assembler = Assembler(ax, conf)
    assembler.run()