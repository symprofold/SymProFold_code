import ctl


class Model():
    '''
    This abstract class describes a molecular model.
    '''

    def __init__(self, model_reg):
        ''' Initialization. '''

        self.model_reg = model_reg
        self.id = None

        self.chimerax_session = None

        return


    def set_id(self, model_id):
        ''' Set model id. '''

        model_id = self.model_reg.convert_model_id(model_id)
        self.id = model_id

        return


    @property
    def idstr(self):
        ''' Get model id in string representation. '''
        model_id_str = self.model_reg.convert_model_id_to_str(self.id)

        return model_id_str


    def set_session(self, session):
        ''' Set Chimerax session. '''
        self.chimerax_session = session # object
        return


    def move_model(self, vect):
        '''
        Move model.
        '''
        self.chimerax_session.move_model(self.id, vect)

        return
