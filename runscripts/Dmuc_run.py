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


# options for runlevel 0
conf.set_flatten_modes([0,1,2,3,4,5], 0)
        # 0: pure superposition, 1: flattened, 2: snapin/tile,
        # 3:completed chains, 4:primitive_unit_cell,
        # 5: assembly of 3x3 primitive unit cells
conf.set_delete_termini_modes([0], 0)


conf.set_species('Dmuc', 'Desulfurococcus mucosus')
conf.set_gene('E8R795')
conf.import_domains()


for conformation in conf.conformations:
    model_reg = ModelReg()
    sess.set_model_reg(model_reg)

    ax = [Axis(model_reg), Axis(model_reg)]

    ax[0].set_session(sess, conf)
    ax[0].set_folder( '23/E8R795_28x4/', 1) # set oriented path    
    ax[0].set_surface([ [442,1000] ])

    ax[1].set_session(sess, conf)
    ax[1].set_folder( '23/E8R795_17x2/', 2) # set oriented path
    ax[1].set_surface([ [1,441] ])

    # ax[0].set_model_active(0)
    # ax[1].set_model_active(0)

    primitive_unit_cell_molecules = ['1.3', '1.4', '12.1', '13.2']

    assembler = Assembler(ax, conf, \
                          primitive_unit_cell_molecules=primitive_unit_cell_molecules)
    assembler.run()
