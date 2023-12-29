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
conf.set_flatten_modes([0,1,2,4], 0)   #[0,1,2,3,4]
        # 0: pure superposition, 1: flattened, 2: snapin/tile,
        # 3:completed chains, 4:primitive_unit_cell
conf.set_delete_termini_modes([1], 0)


conf.set_species('Athr', 'Aneurinibacillus thermoaerophilus')
conf.set_gene('Q6TL21')
conf.import_domains()

conf.set_export_path_filters([conf.export_file_prefix])


for conformation in conf.conformations:
    model_reg = ModelReg()
    sess.set_model_reg(model_reg)

    ax = [Axis(model_reg), Axis(model_reg)]

    ax[0].set_session(sess, conf)
    ax[0].set_folder('af23/Q6TL21_13x4/', 1) #set oriented path
    ax[0].set_surface([ [1,601] ])

    ax[1].set_session(sess, conf)
    ax[1].set_folder('af23/Q6TL21_24x4/', 1) #set oriented path
    ax[1].set_surface([ [602,1000] ])                     

    # ax[0].set_model_active(0)
    # ax[1].set_model_active(0)

    assembler = Assembler(ax, conf)
    assembler.run()