import ctl
import geometry
import math
import re

import os

from chimerax.atomic import all_atomic_structures
from chimerax.core.commands import run


class ChimeraxSession():
    '''
    This class describes an API to a ChimeraX session.
    '''


    def __init__(self, session, silent=0):
        ''' Initialization of the ChimeraxSession class. '''

        self.session = session
        self.model_reg = -1

        if silent == 0:
            self.init()

        return


    def init(self):
        '''
        Initialization of new session in the ChimeraX application.
        '''
        run(self.session, 'close session')
        run(self.session, 'camera ortho')
        run(self.session, 'set bgColor white')
        run(self.session, 'cofr 0,0,0 showPivot 10,0.3')        

        return


    def set_model_reg(self, model_reg):
        ''' Set model register. '''

        self.model_reg = model_reg

        return


    def get_structure(self, model_id):
        ''' Get structure model object of ChimeraX session. '''

        id1 = model_id[0]
        id2 = model_id[1]

        models = self.session.models.list()
        ret = []

        for m in models:
            if len(m.id) == len(model_id):

                if m.id[0] == id1 and m.id[1] == id2:
                    ret.append(m)

        if len(ret) != 1:
            ctl.d(models)
            ctl.d('input model_id')
            ctl.d(model_id)
            ctl.d('ret:')

            for r in ret:
                ctl.d(r.id)
            ctl.d(model_id)
            ctl.error('get_structure: no unique structure found to idstr')
        else:

            return ret[0]


    def open_model(self, path):
        '''
        Open model.
        '''
        opened_model_id = 0

        # check if file exists
        if not os.path.exists(path):
            ctl.e(path)
            ctl.error('ChimeraxSession: open_model: file does not exist.')

        self.run('open "'+path+'"')
        opened_model_id = (self.last_id(),)

        return opened_model_id


    def get_xyz(self, model_id, resid, atomid='CA', \
                error_when_res_not_found=True):
        ''' Get coordinates (x, y, z) of a residue (model_id, resid). '''

        model_id = self.model_reg.convert_model_id(model_id)

        struct = self.get_structure(model_id)
        residues = struct.residues

        res_found = []
        
        for r in residues:
            if r.number == resid:
                res_found.append(r)
                break

        if len(res_found) == 0 and error_when_res_not_found:
            ctl.e(model_id)
            ctl.e(resid)
            ctl.error('get_xyz: no residue found for resid')
        elif len(res_found) == 0:
            return False
        elif len(res_found) > 1:
            ctl.e(model_id)
            ctl.e(resid)
            ctl.error('get_xyz: no unique residue found for resid')

        coord_found = []

        for a in res_found[0].atoms:
            if a.name == atomid:
                coord_found.append(a.scene_coord)
                break

        if len(coord_found) != 1:
            ctl.e(idstr)
            ctl.error('get_xyz: no unique coord found for resid')

        return coord_found[0]


    def get_coord_using_getcrd_command(self, idstr):
        '''
        Get coordinates (x, y, z) of a residue (model_id, resid) using the
        "getcrd" command of ChimeraX.
        '''
        re = run(self.session, "getcrd "+str(idstr))

        return re[0]


    def get_interface_residues(self, model_id0, model_id1):
        '''
        Get interface residues between two models.
        '''
        model_id0 = self.model_reg.convert_model_id(model_id0)
        model_id1 = self.model_reg.convert_model_id(model_id1)
        model_id0_str = self.model_reg.convert_model_id_to_str(model_id0)
        model_id1_str = self.model_reg.convert_model_id_to_str(model_id1)
        
        selection = run(self.session, 'interfaces select #'+ \
                            model_id0_str+' contacting #'+ \
                            model_id1_str+' bothSides true')

        residues = []

        for s in selection:
            model = str(s).split(' #')[1].split('/')[0].split('.')
            resid = str(s).strip().split(' ')[-1]
            residues.append([(int(model[0]), int(model[1])), int(resid)])

        return residues


    def match(self, model_id, model_res_range, match_to_id, bring_id):
        ''' Superposition of model to another model. '''

        model_id = self.model_reg.convert_model_id(model_id)
        model_id_str = self.model_reg.convert_model_id_to_str(model_id)

        match_to_id_str = self.model_reg.convert_model_id_to_str(match_to_id)
        bring_id = self.model_reg.convert_model_id(bring_id)     

        model_res_range_txt = ''
        if len(model_res_range) == 2:
            model_res_range_txt = ':'+str(model_res_range[0])+'-'+ \
                                  str(model_res_range[1])

        run(self.session, 'match #'+model_id_str+model_res_range_txt+ \
            ' to #'+match_to_id_str+' bring #'+str(bring_id[0])+'.1-100')

        return


    def run(self, cmd, downgrade_errors=False):
        ''' Run command iwth the "run" command of ChimeraX. '''

        return run(self.session, cmd, downgrade_errors=downgrade_errors)


    def get_seq(self, model_id, selector=3):
        ''' Get sequence of a model. '''

        modelid = self.model_reg.convert_model_id_to_str(model_id)

        aa31 = {'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
                'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', \
                'TYR':'Y', 'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', \
                 'MET':'M', 'ALA':'A', 'GLY':'G', 'PRO':'P', 'CYS':'C'}

        aa13 = {'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', \
                'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', 'F':'PHE', \
                'Y':'TYR', 'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', \
                'M':'MET', 'A':'ALA', 'G':'GLY', 'P':'PRO', 'C':'CYS'}

        ret = []
        seq = ''

        for s in all_atomic_structures(self.session):

            #format: e.g. 1.1
            if '#'+modelid+' ' in str(s)+' ':

                for r in s.residues:
                    arr = str(r).split(' ')
                    res = arr[-2]

                    if len(res) == selector:  
                        seq += aa31[arr[-2]]

        return seq


    def get_res_n(self, model_id, selector=3):
        '''
        Get overall number of residues of a main model including submodels.
        '''
        model_id = self.model_reg.convert_model_id(model_id)

        length = 0

        for s in all_atomic_structures(self.session):
            ctl.p(str(s))

            if ('#'+str(model_id[0])+'.' in str(s)+' ') or \
               ('#'+str(model_id[0])+' ' in str(s)+' '):
                for r in s.residues:
                    arr = str(r).split(' ')
                    res = arr[-2]

                    if len(res) == selector:  
                        length += 1

        return length


    def get_model_res_n(self, model_id, selector=3):
        '''
        Get overall number of residues of a given model_id including submodels
        of this model_id.
        The model_id can be a submodel.
        '''
        model_id = self.model_reg.convert_model_id(model_id)
        model_id_str = self.model_reg.convert_model_id_to_str(model_id)

        length = 0

        for s in all_atomic_structures(self.session):
            ctl.p(str(s))

            if ('#'+model_id_str+'.' in str(s)+' ') or \
               ('#'+model_id_str+' ' in str(s)+' '):
                for r in s.residues:
                    arr = str(r).split(' ')
                    res = arr[-2]

                    if len(res) == selector:  
                        length += 1

        return length


    def last_id(self, range_max=99):
        '''
        Get highest main model id of all models in a ChimeraX session that are
        smaller than range_max.
        '''
        last_model_id = 0

        for m in all_atomic_structures(self.session):
            model_id = str(m).split('#')[1]
            model_id_0 = int(model_id.split('.')[0])
            
            if last_model_id < model_id_0 <= range_max:
                last_model_id = model_id_0

        return last_model_id


    def close_ids_txt(self, model_ids):
        '''
        Close text line with ids (separated by a space) in a ChimeraX session.
        '''
        model_ids = model_ids.split(' ')
        self.close_ids(model_ids)

        return


    def close_ids(self, model_ids):
        ''' Close list of ids in a ChimeraX session. '''

        for model_id in model_ids:
            if model_id != '':
                self.close_id(model_id)

        return


    def close_id(self, model_id, alert=True):
        ''' Close model_id in a ChimeraX session. '''

        model_id = self.model_reg.convert_model_id(model_id)
        model_id_str = self.model_reg.convert_model_id_to_str(model_id)

        ctl.d('close '+str(model_id))

        # avoid closing of whole ChimeraX session when to_close is empty
        if len(model_id_str) >= 1:
            self.run('close #'+model_id_str)


        if self.model_reg.model_exists(model_id) == True:
            if self.model_reg.check_model_writeprotection(model_id) == False:
                self.model_reg.remove_model(model_id)

        return


    def get_chainid(self, modelid):
        ''' Get chain id of a model. '''

        chain = ''

        for s in all_atomic_structures(self.session):
            if '#'+modelid+' ' in str(s)+' ':

                for r in s.residues:
                    chain = str(r).split('/')[1].split(' ')[0]
                    return chain

        return chain


    def resids(self, modelid, selector=3):
        ''' Get residues of a model. '''

        ret = []
        modelid = self.model_reg.convert_model_id_to_str(modelid)

        for s in all_atomic_structures(self.session):
            if '#'+modelid+' ' in str(s)+' ':

             for r in s.residues:
               arr = str(r).split(' ')
               res = arr[-2]

               if res == '4SU' or res == 'H2U' or res == '5MU' or res == 'PSU':
                   res = 'U'
               if res ==' OMC':
                   res = 'C'
               if res == 'G7M':
                   res = 'C'           

               if len(res) == selector or selector == 0:  
                ret.append(int(arr[-1]))

        return ret


    def resids_surface(self, model_id, selector=3):
        ''' Get "surface residues" of model as defined in Axis object. '''

        _resids = self.resids(model_id, selector)

        model = self.model_reg.get_model(model_id)
        model_surface_resids = model.get_surface_resids()

        resids_s = [ r for r in _resids if r in model_surface_resids ]

        return resids_s


    def resids_fragments(self, modelid, selector=3):
        ''' Get fragments of contiguous residues. '''

        ret = self.resids(modelid, selector)
        resid_fragments = []

        lastresid = -1000
        for r in ret:
            if r == lastresid+1:
                resid_fragments[-1].append(r)
            else:    
                resid_fragments.append([r])
            lastresid = r    

        return ret, resid_fragments


    def resids_surface_fragments(self, modelid, selector=3):
        '''
        Get fragments of contiguous "surface residues" of model as defined in
        Axis object.
        '''
        ret = self.resids_surface(modelid, selector)
        resid_fragments = []

        lastresid = -1000
        for r in ret:
            if r == lastresid+1:
                resid_fragments[-1].append(r)
            else:    
                resid_fragments.append([r])
            lastresid = r    

        return ret, resid_fragments


    def get_residue_overlap(self, id1, id2):
        ''' Get overlapping residues. '''
        overlap = []

        id1_resids = self.resids(id1)
        id2_resids = self.resids(id2)

        for res1 in id1_resids:
            for res2 in id2_resids:
                if res1 == res2:
                    overlap.append(res1)
                    break

        return overlap


    def get_submodel_ids(self, model_id):
        ''' Get ids of all submodels of given model. '''

        model_id2 = self.model_reg.convert_model_id(model_id)
        subids = []

        for s in all_atomic_structures(self.session):
            if '#'+str(model_id2[0])+'.' in str(s):
                subid = str(s).split('#')[1]
                subids.append( self.model_reg.convert_model_id(subid) )

        return subids


    def get_rot_axis_z_dist(self, model_id, model_id_ref):
        '''
        Get distance of rotation axis from molecule center.
        '''
        model_id = self.model_reg.convert_model_id(model_id)

        model_id_str = self.model_reg.convert_model_id_to_str(model_id)
        model_id_ref_str = self.model_reg.convert_model_id_to_str(model_id_ref)

        re = self.run('measure rotation #'+model_id_str+' toModel #'+ \
                      model_id_ref_str+' showAxis false')
        rot_axis_point = get_rot_axis_point(re.description())
        cen = self.model_reg.get_model((model_id[0],)).get_center()
        dist = geometry.dist(rot_axis_point, cen)

        return dist


    def get_rot_angle_z(self, model_id, model_id_ref):
        '''
        Get rotation angle around z axis of model_id regarding model_id_ref.
        '''
        model_id = self.model_reg.convert_model_id(model_id)
        model_id_ref = self.model_reg.convert_model_id(model_id_ref)

        model_id_str = self.model_reg.convert_model_id_to_str(model_id)
        model_id_ref_str = self.model_reg.convert_model_id_to_str(model_id_ref)

        angle = 0
        re = self.run('measure rotation #'+model_id_str+' toModel #'+ \
                      model_id_ref_str+' showAxis false')
        rot_axis = get_rot_axis(re.description())
        angle = re.rotation_angle()/math.pi*180

        rot_axis_point = get_rot_axis_point(re.description())
        cen = self.model_reg.get_model((model_id[0],)).get_center()
        dist = geometry.dist(rot_axis_point, cen)

        angle_z = angle*rot_axis[2]

        return angle_z, dist


    def ids_str(self, model_objects):
        '''
        Create ChimeraX string of model ids from given list of model objects.
        '''
        ids_str_ = ''

        for m in model_objects:
            ids_str_ += ('#'+m.idstr+' ')

        ids_str_ = ids_str_[:-1]

        return ids_str_


    def combine_models(self, models, combined_model_id):
        '''
        Combine models to model_id_combined.
        '''
        models_idstr = self.ids_str(models)
        combined_model_id = self.model_reg.convert_model_id(combined_model_id)

        self.run('combine '+models_idstr+ \
                        ' close false modelId #'+str(combined_model_id[0]))

        return


    def combine(self, ids, intermediate_id):
        '''
        Combine all models in ChimeraX session and code ids in chainid.

        Splitting can be done with split() function.
        '''
        idstr = ''
        ids_to_combine = []
        chainids_new = []

        for i in ids:

            for s in all_atomic_structures(self.session):
                if '#'+str(i)+'.' in str(s):
                    chainid = self.get_chainid(str(s.id[0])+'.'+str(s.id[1]))
                    chainid_new = 'mid'+str(s.id[0])+'point'+str(s.id[1])+ \
                            'cid'+str(chainid)+'end'
                    chainids_new.append(chainid_new)
                    ret = self.run('changechains '+ \
                            '#'+str(s.id[0])+'.'+str(s.id[1])+' '+chainid_new)

                    idstr = idstr+'#'+str(s.id[0])+'.'+str(s.id[1])+' '
                    ctl.d(s.name)
                    ctl.d(s.id)
                    ids_to_combine.append( s.id )
                    ctl.d(s)

        ret = self.run('combine '+idstr+' modelId #'+str(intermediate_id))
        
        return chainids_new


    def split(self, id_to_split):
        '''
        Split combined models in separate models and use ids coded in chainid.

        Combination can be done with combine() function.
        '''
        ret = self.run('split #'+str(id_to_split)+'')

        for s in all_atomic_structures(self.session):
            if '#'+str(id_to_split)+'.' in str(s):
                mid = str(s).split(' mid')[1].split('cid')[0]. \
                          replace('point', '.')
                cid = str(s).split('cid')[1].split('end')[0]

                ctl.d(mid)
                ctl.d(cid)

                self.run('rename #'+str(s.id[0])+'.'+str(s.id[1])+ \
                         ' id #'+mid+' ')
                self.run('changechains '+'#'+str(s.id[0])+'.'+str(s.id[1])+ \
                         ' '+cid)

        self.close_id(id_to_split)

        return


    def set_marker(self, points, start_id):
        ''' Set marker to each point in list. '''

        for i,p in enumerate(points):
            idnr = start_id+i
            self.run('marker #'+str(idnr)+' position '+ \
                str(p[0])+','+str(p[1])+','+str(p[2])+' color red radius 10')

        return


def get_transl_vec(txt):
    ''' Get translation vector from ChimeraX text description. '''

    vec = [0,0,0]
    txt = txt.split("\n")

    vec[0] = txt[1].strip()
    vec[1] = txt[2].strip()
    vec[2] = txt[3].strip()

    # replace multiple spaces by a single space
    vec = [re.sub(r'\s{2,}', ' ', v) for v in vec]

    vec = [v.split(' ')[3] for v in vec]
    vec = [float(v) for v in vec]

    return vec


def get_rot_axis(txt):
    ''' Get rotation axis from ChimeraX text description. '''

    vec = [0,0,0]
    txt = txt.split("\n")
    v0 = txt[4].strip()

    # replace multiple spaces by a single space
    v0 = re.sub(r'\s{2,}', ' ', v0)
    
    v0 = v0.split(' ')

    vec = [v0[1], v0[2], v0[3]]
    vec = [float(v) for v in vec]

    return vec


def get_rot_axis_point(txt):
    ''' Get rotation axis point from ChimeraX text description. '''

    vec = [0,0,0]
    txt = txt.split("\n")
    v0 = txt[5].strip()

    # replace multiple spaces by a single space
    v0 = re.sub(r'\s{2,}', ' ', v0)

    v0 = v0.split(' ')

    vec = [v0[2], v0[3], v0[4]]
    vec = [float(v) for v in vec]

    return vec
