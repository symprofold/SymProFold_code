import bibpdb
import chimerax_api
import ctl
import filesystem
import geometry
import molmodel
import proteinmerge
from structure.model import Model

import math
import os


def open_model(sess, file, model_id, meta={}, part=0):
    '''
    Open model and calculate basic parameter.
    '''
    if part == 0 or part == 1:
        sess.run('open '+str(file))  
        model_id += 1
    if part == 0 or part == 2:
        sess.run('hide #'+str(model_id)+' atoms')
        sess.run('show #'+str(model_id)+' cartoons')
        sess.run('split #'+str(model_id))
    if part == 0 or part == 1:
        rmsds = bibpdb.get_rmsds(file)
        termini = molmodel.get_termini(rmsds)

        multimer_n = get_multimer_n(file)

        meta[model_id] = [[], termini, multimer_n, [[], []]]
      
    return model_id, meta


def format_model(model_id, ax):
    '''
    Perform formatting operations on model.
    '''
    for h in ax.hide:
        ax.chimerax_session.run(session,
                'hide #'+str(model_id)+':'+str(h[0])+'-'+str(h[1])+' cartoons')

    return


def get_multimer_n(file):
    '''
    Get number of chains in a given pdb file.
    Use case: e.g. get multiplicity of a multimer in a pdb file.
    '''
    n = 0
    f = filesystem.get_file(file)

    # remove all lines except 'ATOM', 'END' and 'TER'
    code_cleaned = ['ATOM', 'END', 'TER']
    f = [l for l in f if any(l.startswith(code) for code in code_cleaned)]

    for i, l in enumerate(f):
        if i >= 1:   
            if f[i][0:3] == 'TER' or \
                    (f[i][0:3] == 'END' and f[i-1][0:3] != 'END' and \
                    f[i-1][0:3] != 'TER'):
                n += 1

    return n


def place_plane(session):
    ''' Place 3 orientation points to allow identification of xy plane. '''

    session.run('marker #100 position 0,0,0 color red radius 1')
    session.run('marker #101 position 100,0,0 color red radius 1')
    session.run('marker #102 position 0,100,0 color red radius 1')

    return


def model_2fold_to_plane(model_id, meta, session):
    ''' Align 2-fold axis model to xy plane. '''

    model = Model(session.model_reg)
    model.set_id(model_id)

    tmp_dir = filesystem.clean_path( \
            os.path.dirname(os.path.realpath(__file__))+'/../../')
    tmpfile = tmp_dir+'tmp.pdb'
    session.run('save '+tmpfile+' models #'+str(model_id))

    tmp_model = Model(session.model_reg)
    tmp_model.set_id((model_id+1,))
    session.run('open "'+tmpfile+'"')
    session.run('split #'+tmp_model.idstr)
    session.run('match #'+tmp_model.idstr+'.2 to #'+tmp_model.idstr+'.1')
    re = session.run('measure rotation #'+tmp_model.idstr+'.2 toModel #'+ \
                     tmp_model.idstr+'.1'+' showAxis false')
    rot_axis = chimerax_api.get_rot_axis(re.description())
    session.close_id(tmp_model.id)

    ctl.d(rot_axis)

    center_res = round((meta[model_id][1][0]+ \
                        meta[model_id][1][1])/2)
    center0 = session.get_coord_using_getcrd_command(
                    '#'+model.idstr+'/A:'+str(center_res)+'@CA')
    center1 = session.get_coord_using_getcrd_command(
                    '#'+model.idstr+'/B:'+str(center_res)+'@CA')

    c0c1 = [center1[0]-center0[0], center1[1]-center0[1], \
            center1[2]-center0[2]]

    normal = geometry.crossproduct(c0c1, rot_axis)
    ctl.d(center_res)
    ctl.d(normal)

    main_dir = filesystem.clean_path( \
            os.path.dirname(os.path.realpath(__file__))+'/../')

    session.run('open "'+main_dir+'/component_files/marker.pdb'+'"')

    session.run('move x '+str(center0[0]+normal[0])+ \
                ' models #'+tmp_model.idstr)
    session.run('move y '+str(center0[1]+normal[1])+ \
                ' models #'+tmp_model.idstr)
    session.run('move z '+str(center0[2]+normal[2])+ \
                ' models #'+tmp_model.idstr)    

    session.combine_models([model, tmp_model], str(tmp_model.id[0]+1))
    session.run('marker #201 position '+str(normal[0])+','+ \
                str(normal[1])+','+str(normal[2])+' color yellow radius 1')

    cmd_align = 'align #'+str(tmp_model.id[0]+1)+'/A:'+str(center_res)+'@CA '+ \
                 '#'+str(tmp_model.id[0]+1)+'/B:'+str(center_res)+'@CA '+ \
                 '#'+str(tmp_model.id[0]+1)+'/C:1@CA '+ \
                 'toAtoms #100,101,102'
    ctl.d(cmd_align)
    session.run(cmd_align)

    session.close_id(tmp_model.id)

    session.run('split #'+str(tmp_model.id[0]+1))
    session.close_id((tmp_model.id[0]+1, 3))

    ctl.d(model.idstr+'.1')
    ctl.d(str(tmp_model.id[0]+1)+'.1')

    session.match(model.id, [], tmp_model.id[0]+1, model.id, \
                  model_chainid='A', match_to_chainid='A')

    session.close_id(tmp_model.id[0]+1)

    return


def model_to_plane(meta, model_id, session):
    ''' Align ">=2"-fold axis model to xy plane. '''

    center_res = round((meta[model_id][1][0]+meta[model_id][1][1])/2)

    if meta[model_id][2] == 2:
        model_2fold_to_plane(model_id, meta, session)

    if meta[model_id][2] >= 3:
        cmd_align = 'align #'+str(model_id)+'/A:'+str(center_res)+'@CA '+ \
                          '#'+str(model_id)+'/B:'+str(center_res)+'@CA '+ \
                          '#'+str(model_id)+'/C:'+str(center_res)+'@CA '+ \
                          'toAtoms #100,101,102'

        session.run(cmd_align)

    return


def combine_chains(model1_id, model2_id, session, dest_id=500, \
                   duplicate_after_completion=False):
    '''
    Combine two chains from two models to a combined model.

    Args:
        duplicate_after_completion: The combined model is copied to both
                input model ids after combination.
    '''
    model1_id = session.model_reg.convert_model_id(model1_id)
    model2_id = session.model_reg.convert_model_id(model2_id)


    # checks
    if model1_id == model2_id:
        ctl.e(model1_id)
        ctl.e(model2_id)
        ctl.error('combine_chains: same model id')
        return False        

    if not session.model_reg.model_exists(model1_id):
        ctl.e(model1_id)
        ctl.error('combine_chains: model missing')

    if not session.model_reg.model_exists(model2_id):
        ctl.e(model2_id)
        ctl.error('combine_chains: model missing')


    model1 = session.model_reg.get_model(model1_id)
    model2 = session.model_reg.get_model(model2_id)

    model_1_surface = model1.get_surface()
    model_2_surface = model2.get_surface()

    conf = session.model_reg.get_model(model1.id[0]).conf

    # check if gene contains only one domain
    if len(conf.domains) == 1:
        ctl.d(conf)
        ctl.d(conf.domains)
        return False


    # check if models already merged
    if model1.id in session.model_reg.get_model(model1.id[0]).connected:
        ctl.d(model1.id)
        ctl.d('combine_chains: already connection existing with model1.id')
        return False

    if model2.id in session.model_reg.get_model(model2.id[0]).connected:
        ctl.d(model2.id)
        ctl.d('combine_chains: already connection existing with model2.id')
        return False


    resids1 = session.resids(model1.id)
    resids2 = session.resids(model2.id)

    if len(resids1) == 0 or len(resids2) == 0:
        ctl.d('combine_chains not possible because one chain empty')
        return False


    for i, sur in enumerate(model_1_surface):
        if i == 0:
            session.delete_res_range( \
                        model1.id, [-1, model_1_surface[0][0]-1])
            session.delete_res_range( \
                        model1.id, [model_1_surface[-1][1]+1, 10000])

        if i < len(model_1_surface)-1:
            session.delete_res_range(model1.id, \
                        [model_1_surface[i][1]+1, model_1_surface[i+1][0]-1])

    for i, sur in enumerate(model_2_surface):
        if i == 0:
            session.delete_res_range( \
                        model2.id, [-1, model_2_surface[0][0]-1])
            session.delete_res_range( \
                        model2.id, [model_2_surface[-1][1]+1, 10000])

        if i < len(model_2_surface)-1:
            session.delete_res_range(model2.id, \
                        [model_2_surface[i][1]+1, model_2_surface[i+1][0]-1])


    # checks
    res_ov = session.get_residue_overlap(model1.id, model2.id)
    if len(res_ov) != 0:
        ctl.e('combine_chains: chains overlapping, not possible to merge')
        ctl.e(model1.id)
        ctl.e(model2.id)
        ctl.e(res_ov)
        ctl.error('combine_chains: chains overlapping, not possible to merge')


    dist_check, dist = proteinmerge.check_merge_possible( \
                                            model1.id, model2.id, session, 80)

    if dist_check == 0:
        ctl.e(model1.id)
        ctl.e(model2.id)        
        ctl.e(dist_check)
        ctl.e(dist)
        raise Exception('combine_chains: dist_check failed')

    id1_chainid = session.get_chainid(model1.id)
    id2_chainid = session.get_chainid(model2.id)
    session.rename_chainid(model1.id, id1_chainid)
    session.rename_chainid(model2.id, id1_chainid)

    session.combine_models([model2, model1], dest_id)

    session.close_id(model1.id)
    session.close_id(model2.id)

    session.change_model_id(dest_id, model1.id)
    session.rename_chainid(model1.id, id1_chainid)

    # duplicate after completion
    if duplicate_after_completion:
        session.combine_models([model1], model2.idstr)
        session.rename_chainid(model2.id, id2_chainid)


    if session.model_reg.model_exists(model1.id):
        model1.modelling_completeness = 2
        session.model_reg.get_model(model1.id[0]).connected.append(model1.id)

    if session.model_reg.model_exists(model2.id):
        model2.modelling_completeness = 2
        session.model_reg.get_model(model2.id[0]).connected.append(model2.id)

    return True
