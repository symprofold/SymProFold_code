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
conf.set_flatten_modes([0,1,2,3,4], 0)    # 0,1,2,3,4
        # 0: pure superposition, 1: flattened, 2: snapin/tile,
        # 3:completed chains, 4:primitive_unit_cell
conf.set_snapshot_modes([0,1], 0)   #0,1     snapshot_w_separated_chains
conf.set_delete_termini_modes([0,1], 0)


conf.set_species('Tkiv', 'Thermoanaerobacter kivui')
conf.set_gene('P22258')
conf.import_domains()


for conformation in  conf.conformations:
    model_reg = ModelReg()
    sess.set_model_reg(model_reg)

    ax = [Axis(model_reg), Axis(model_reg)]

    ax[0].set_session(sess, conf)
    ax[0].set_folder('P22258_15x6/', 1) #set oriented path
    # ax[0].set_domains([ [416,478] ])
    ax[0].set_domains([ [280,346] ])
    ax[0].set_surface([ [1,345] ])

    ax[1].set_session(sess, conf)
    ax[1].set_folder('P22258_29x2/', 1) #set oriented path
    # ax[1].set_domains([ [348,478] ])
    #ax[1].set_surface([ [479,1000] ])
    ax[1].set_domains([ [280,346] ])
    ax[1].set_surface([ [346,1000] ])

    ax[0].set_model_active(2)
    ax[1].set_model_active(0)

    assembler = Assembler(ax, conf)
    assembler.run()
