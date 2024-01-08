import ctl
import export
import validation


class PrimitiveUnitCell():
    '''
    This class describes a primitive unit cell.
    '''

    def __init__(self, layer):
        ''' Initialization of the PrimitiveUnitCell class. '''

        self.layer = layer
        self.molecules = []
        self.lc_offset = 0
        self.a = None
        self.b = None
        self.c = None
        self.alpha = None
        self.beta = None
        self.gamma = None
        self.symmgroup = None

        return


    def add_molecule(self, model_id):
        ''' Add a monomer molecule to primitive unit cell. '''

        self.molecules.append(model_id)

        return


    def definition(self, lc_offset=0):
        ''' Create definition of primitive unit cell. '''

        if not(self.layer.symmgroup != '' and len(self.layer.axes) == 2):
            return False

        self.lc_offset = lc_offset
        self.a = self.layer.calc_lattice_constant(self.lc_offset)
        self.b = self.layer.calc_lattice_constant(self.lc_offset)
        self.c = 1000000
        self.alpha = 90
        self.beta = 90
        self.gamma = 90
        self.symmgroup = 'P 1'

        if self.layer.symmgroup in ('p3', 'p6'):
            self.gamma = 60

        return


    def suggest_molecules(self):
        ''' Suggest primitive unit cell molecules. '''

        if len(self.layer.axes) > 1:
            if self.layer.axes[0].fold == 3 and self.layer.axes[1].fold == 2:
                self.molecules = []
                self.add_molecule('1.1')
                self.add_molecule('1.2')
                self.add_molecule('1.3')
                self.add_molecule('8.3')
                self.add_molecule('9.1')
                self.add_molecule('10.2')

            if self.layer.axes[0].fold == 3 and self.layer.axes[1].fold == 3:
                self.molecules = []
                self.add_molecule('1.1')
                self.add_molecule('1.2')
                self.add_molecule('1.3')

            if self.layer.axes[0].fold == 4 and self.layer.axes[1].fold == 2:
                self.molecules = []
                self.add_molecule('1.3')
                self.add_molecule('1.4')
                self.add_molecule('12.1')
                self.add_molecule('12.2') 

            if self.layer.axes[0].fold == 4 and self.layer.axes[1].fold == 4:
                self.molecules = []
                self.add_molecule('1.1')
                self.add_molecule('1.2')
                self.add_molecule('11.3')
                self.add_molecule('11.4')             

        if self.layer.axes[0].fold == 6:
            self.molecules = []
            self.add_molecule('1.1')
            self.add_molecule('1.2')
            self.add_molecule('1.3')
            self.add_molecule('1.4')
            self.add_molecule('1.5')
            self.add_molecule('1.6')

        return


    def combine_to_model(self, combination_model_id):
        ''' Combine molecules of primitive unit cell. '''

        if not(self.layer.symmgroup != '' and len(self.layer.axes) == 2):
            return False

        comb_str = ''
        for m in self.molecules:
            comb_str += ' #'+str(m)

        ax0 = self.layer.axes[0]
        ax0.chimerax_session.run('combine '+comb_str+ \
                ' close false modelId #'+str(combination_model_id))

        if self.layer.symmgroup in ('p3', 'p6'):
            if not (self.layer.axes[0].fold == 3 and \
                    self.layer.axes[1].fold == 2):
                ax0.chimerax_session. \
                    run('turn z 30 models #'+str(combination_model_id))

        return True


    def export(self, combination_model_id, export_path, conf):
        ''' Export primitive unit cell. '''

        if not(self.layer.symmgroup != '' and len(self.layer.axes) == 2):
            return False

        ax0 = self.layer.axes[0]

        # validate export combined model
        validation.validation(combination_model_id, \
                              self.layer.axes, \
                              export_path, \
                              conf, \
                              ax0.chimerax_session)

        export.compatibility_cif_export_combine( \
                export_path, \
                combination_model_id, \
                ax0.chimerax_session, \
                conf.cif_postprocess)
        ctl.p('cif_add_symm')
        export.cif_add_symm(export_path, export_path, \
                            self.a, self.b, self.c, \
                            self.alpha, self.beta, self.gamma, \
                            self.symmgroup)
        ctl.p('cif_add_symm finished')

        return
