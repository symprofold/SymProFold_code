import ctl


class Monomer():
    '''
    This class describes a monomer. Monomers can be parts of larger complexes.
    '''


    def __init__(self, model_reg, model_id=-1, sequence=''):
        ''' Initialization of the Monomer class. '''

        self.model_reg = model_reg
        if model_id != -1:
            self.model_id = self.model_reg.convert_model_id(model_id)
        self.sequence = sequence
        self.modelling_completeness = 0
            # 0: unknown; 1: partly; 2: sufficient; 3: full sequence modelled

        self.missing_res_n = -1

        return


    def set_model_id(self, model_id):
        ''' Set model id. '''

        model_id = self.model_reg.convert_model_id(model_id)
        self.model_id = model_id

        return


    def get_missing_res_n(self, model_id):
        ''' Get number of missing residues. '''

        if self.missing_res_n == -1:
            res_n = self.chimerax_session.get_seq(model_id)
            self.missing_res_n = len(self.sequence)-res_n
            
        return self.missing_res_n


    def set_session(self, session):
        ''' Set Chimerax session. '''
        self.chimerax_session = session # class
        self.session = session.session
        return


    def get_surface(self):
        ''' Get residue range that will be included in an overall model. '''

        model_0 = self.model_reg.get_model((self.model_id[0],))
        surface = model_0.surface

        return surface


    def get_surface_resids(self):
        ''' Get residues that will be included in an overall model. '''

        ctl.d('self.model_id')
        ctl.d(self.model_id)

        model_0 = self.model_reg.get_model((self.model_id[0],))
        surface = model_0.surface_resids

        ctl.d('surface_resids')
        ctl.d(surface)

        return surface


    def get_passive_resids(self):
        ''' Get residues that will not be included in an overall model. '''

        ctl.d('self.model_id')
        ctl.d(self.model_id)

        ctl.d(model_0)
        model_0 = self.model_reg.get_model((self.model_id[0],))
        surface = model_0.surface_resids

        return surface


def get_center(model_id, res_range, sess):
    ''' Get center of monomer model. '''

    center = [0,0,0]
    sum_vect = [0,0,0]
    sum_n = 0        

    model_id = sess.model_reg.convert_model_id(model_id)
    resids = sess.resids(model_id)

    for resid in resids:

        # sum of all residues inside given range
        if res_range[0] <= resid <= res_range[1]:
            sum_vect = [sum_vect[i]+c
                        for i,c in enumerate(sess.get_xyz(model_id, resid))]
            sum_n += 1

    center[0] = sum_vect[0]/sum_n
    center[1] = sum_vect[1]/sum_n
    center[2] = sum_vect[2]/sum_n    

    return center
