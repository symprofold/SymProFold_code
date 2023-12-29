import ctl


class ModelReg():
    '''
    This class describes a model register referencing model ids to objects.
    '''


    def __init__(self):
        ''' Initialization of the ModelReg class. '''
        self.models = {}
        

    def add_model(self, model_id, model_obj, model_type, writeprotection=True):
        '''
        Add model to model register.

        model types: 'monomer', 'axis'
        '''
        model_id2 = self.convert_model_id(model_id)

        if writeprotection == True:
            try:
                r = self.models[model_id2]
                ctl.e(self.models)
                ctl.e(model_id2)
                ctl.error('add_model: model id already existing.')
            except KeyError:
                self.models[model_id2] = [model_obj, model_type, \
                                          writeprotection]
        else:
            self.remove_model(model_id2)
            self.models[model_id2] = [model_obj, model_type, writeprotection]

        return


    def remove_model(self, model_id):
        '''
        Remove model including submodels from register.
        '''
        model_id2 = self.convert_model_id(model_id)

        if len(model_id2) == 1:

            _models = []
            for m in self.models:
                _models.append(m)

            for m in _models:
                if m[0] == model_id2[0]:
                    del self.models[m]

        elif len(model_id2) == 2:

            if self.model_exists(model_id2):
                del self.models[model_id2]

        return


    def get_model(self, model_id):
        ''' Get model object from model_id. '''

        ctl.d('inp:')
        ctl.d(model_id)
        model_id = self.convert_model_id(model_id)

        try:
            r = self.models[model_id][0]
            return r
        except KeyError:
            ctl.e(model_id)
            ctl.e(self.models)
            ctl.error('ModelReg.get_model: model id not existing')

        return False


    def model_exists(self, model_id):
        ''' Check whether model_id exists in ModelReg. '''

        model_id = self.convert_model_id(model_id)

        try:
            r = self.models[model_id][0]
            return True
        except KeyError:
            return False

        return False


    def get_submodels(self, model_id1):
        ''' Get submodels from model_id of main model. '''

        submodels = []
        model_id1 = self.convert_model_id(model_id1)

        for m in self.models:
            if m[0] == model_id1[0]:
                submodels.append(self.models[m])

        return submodels


    def get_model_type(self, model_id):
        '''
        Get model type.

        model types: 'monomer', 'axis'
        '''
        model_id = self.convert_model_id(model_id1)

        try:
            r = self.models[model_id][1]
            return r
        except KeyError:
            ctl.error('ModelReg.get_model_type: model id not existing')

        return False


    def check_model_writeprotection(self, model_id):
        ''' Check if model is write-protected. '''

        model_id = self.convert_model_id(model_id)

        try:
            r = self.models[model_id][2]
            return r
        except KeyError:
            ctl.error('ModelReg.check_model_writeprotection: '+ \
                      'model id not existing')

        return


    def convert_model_id(self, model_id):
        ''' Convert model_id to tuple. '''

        if type(model_id) is tuple:
            for i,m in enumerate(model_id):
                if type(model_id[i]) is str:
                    model_id[i] = int(model_id[i])
                
                if type(model_id[i]) is not int:
                    ctl.e(type(model_id[i]))
                    ctl.e(model_id)
                    ctl.error('ModelReg.convert_model_id: '+ \
                              'number in tuple not int')
            return model_id

        if type(model_id) is list:
            for i,m in enumerate(model_id):
                if type(model_id[i]) is str:
                    model_id[i] = int(model_id[i])
                
                if type(model_id[i]) is not int:
                    ctl.e(type(model_id[i]))
                    ctl.e(model_id)
                    ctl.error('ModelReg.convert_model_id: '+ \
                              'number in tuple not int')
            return tuple(model_id)

        if type(model_id) is int:
            return (model_id,)

        if type(model_id) is str:
            if model_id[0] == '#':
                model_id = model_id[1:]
                
            ret = model_id.split('.')
            if len(ret) == 2:
                ret2 = (int(ret[0]), int(ret[1]))
                return ret2
            if len(ret) == 1:
                ret2 = (int(ret[0]),)
                return ret2
            if len(ret) == 3:
                ctl.error('ModelReg.convert_model_id: length > 2')
            

    def convert_model_id_to_str(self, model_id):
        ''' Convert model_id to string representation. '''

        model_id = self.convert_model_id(model_id)
        model_id_str = '.'.join([str(m) for m in model_id])

        return model_id_str      
