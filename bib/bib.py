import bibpdb
import chimerax_api
import ctl
import filesystem
import geometry
import molmodel
import proteinmerge

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

    tmp_dir = filesystem.clean_path( \
            os.path.dirname(os.path.realpath(__file__))+'/../../')
    tmpfile = tmp_dir+'tmp.pdb'
    session.run('save '+tmpfile+' models #'+str(model_id))

    tmp_model_id = model_id+1
    session.run('open '+tmpfile)
    session.run('split #'+str(tmp_model_id))
    session.run('match #'+str(tmp_model_id)+'.2 to #'+str(tmp_model_id)+'.1')
    re = session.run('measure rotation #'+str(tmp_model_id)+'.2 toModel #'+ \
                     str(tmp_model_id)+'.1'+' showAxis false')
    rot_axis = chimerax_api.get_rot_axis(re.description())
    session.close_id(tmp_model_id)

    ctl.d(rot_axis)

    center_res = round((meta[model_id][1][0]+ \
                        meta[model_id][1][1])/2)
    center0 = session.get_coord_using_getcrd_command(
                    '#'+str(model_id)+'/A:'+str(center_res)+'@CA')
    center1 = session.get_coord_using_getcrd_command(
                    '#'+str(model_id)+'/B:'+str(center_res)+'@CA')

    c0c1 = [center1[0]-center0[0], center1[1]-center0[1], \
            center1[2]-center0[2]]

    normal = geometry.crossproduct(c0c1, rot_axis)
    ctl.d(center_res)
    ctl.d(normal)

    main_dir = filesystem.clean_path( \
            os.path.dirname(os.path.realpath(__file__))+'/../')

    session.run('open '+main_dir+'/component_files/marker.pdb')

    session.run('move x '+str(center0[0]+normal[0])+ \
                ' models #'+str(tmp_model_id))
    session.run('move y '+str(center0[1]+normal[1])+ \
                ' models #'+str(tmp_model_id))
    session.run('move z '+str(center0[2]+normal[2])+ \
                ' models #'+str(tmp_model_id))    

    session.run('combine #'+str(model_id)+' #'+str(tmp_model_id))
    session.run('marker #201 position '+str(normal[0])+','+ \
                str(normal[1])+','+str(normal[2])+' color yellow radius 1')

    cmd_align = 'align #'+str(tmp_model_id+1)+'/A:'+str(center_res)+'@CA '+ \
                 '#'+str(tmp_model_id+1)+'/B:'+str(center_res)+'@CA '+ \
                 '#'+str(tmp_model_id+1)+'/C:1@CA '+ \
                 'toAtoms #100,101,102'
    ctl.d(cmd_align)
    session.run(cmd_align)

    session.close_id(tmp_model_id)

    session.run('split #'+str(tmp_model_id+1))
    session.close_id((tmp_model_id+1, 3))

    ctl.d(str(model_id)+'.1')
    ctl.d(str(tmp_model_id+1)+'.1')

    session.run('match #'+str(model_id)+'/A'+' to #'+ \
                str(tmp_model_id+1)+'/A'+' bring #'+str(model_id)+'.1-100')
    
    session.close_id(tmp_model_id+1)

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
                   close_model2=True):
    ''' Combine two chains to one model. '''

    ids_to_close = ''

    id1 = session.model_reg.convert_model_id_to_str(model1_id)
    id2 = session.model_reg.convert_model_id_to_str(model2_id)

    model1_id = session.model_reg.convert_model_id(model1_id)
    model2_id = session.model_reg.convert_model_id(model2_id)


    # checks
    if model1_id == model2_id:
        ctl.e(model1_id)
        ctl.e(model2_id)
        ctl.error('combine_chains: same model id')
        return False        

    if not session.model_reg.model_exists(id1):
        ctl.e(id1)
        ctl.error('combine_chains: model missing')

    if not session.model_reg.model_exists(id2):
        ctl.e(id2)
        ctl.error('combine_chains: model missing')


    model_1_surface = session.model_reg.get_model(id1).get_surface()
    model_2_surface = session.model_reg.get_model(id2).get_surface()

    conf = session.model_reg.get_model(model1_id[0]).conf

    # check if gene contains only one domain
    if len(conf.domains) == 1:
        ctl.d(conf)
        ctl.d(conf.domains)
        return False


    # check if models already merged
    if model1_id in session.model_reg.get_model(model1_id[0]).connected:
        ctl.d(model1_id)
        ctl.d('combine_chains: already connection existing with model1_id')
        return False

    if model2_id in session.model_reg.get_model(model2_id[0]).connected:
        ctl.d(model2_id)
        ctl.d('combine_chains: already connection existing with model2_id')
        return False


    resids1 = session.resids(id1)
    resids2 = session.resids(id2)

    if len(resids1) == 0 or len(resids2) == 0:
        ctl.d('combine_chains not possible because one chain empty')
        return False


    for i, sur in enumerate(model_1_surface):
        if i == 0:
            session.run('delete #'+str(id1)+ \
                ':-1-'+str(model_1_surface[0][0]-1))
            session.run('delete #'+str(id1)+ \
                ':'+str(model_1_surface[-1][1]+1)+'-10000')

        if i < len(model_1_surface)-1:
            session.run('delete #'+str(id1)+ \
                ':'+str(model_1_surface[i][1]+1)+'-'+ \
                str(model_1_surface[i+1][0]-1))

    for i, sur in enumerate(model_2_surface):
        if i == 0:
            session.run('delete #'+str(id2)+ \
                ':-1-'+str(model_2_surface[0][0]-1))
            session.run('delete #'+str(id2)+ \
                ':'+str(model_2_surface[-1][1]+1)+'-10000')

        if i < len(model_2_surface)-1:
            session.run('delete #'+str(id2)+ \
                ':'+str(model_2_surface[i][1]+1)+'-'+ \
                str(model_2_surface[i+1][0]-1))


    # checks
    res_ov = session.get_residue_overlap(id1, id2)
    if len(res_ov) != 0:
        ctl.e('combine_chains: chains overlapping, not possible to merge')
        ctl.e(id1)
        ctl.e(id2)
        ctl.e(res_ov)
        ctl.error('combine_chains: chains overlapping, not possible to merge')


    dist_check, dist = proteinmerge.check_merge_possible(id1, id2, session, 80)

    if dist_check == 0:
        ctl.e(id1)
        ctl.e(id2)        
        ctl.e(dist_check)
        ctl.e(dist)
        ctl.error('combine_chains: dist_check failed')

    id1_chainid = session.get_chainid(id1)
    session.run('changechains #'+str(id1)+' '+id1_chainid)
    session.run('changechains #'+str(id2)+' '+id1_chainid)

    session.run('combine #'+id2+' #'+id1+' modelId #'+str(dest_id))

    session.close_id(id1)

    if close_model2 == True:
        session.close_id(id2)
    else:    
        ids_to_close += ' #'+str(id2)

    session.run('rename #'+str(dest_id)+' id #'+id1)
    session.run('changechains #'+str(id1)+' '+id1_chainid)

    if session.model_reg.model_exists(id1):
        session.model_reg.get_model(id1).modelling_completeness = 2
        session.model_reg.get_model(model1_id[0]).connected.append(model1_id)

    if session.model_reg.model_exists(id2):
        session.model_reg.get_model(id2).modelling_completeness = 2
        session.model_reg.get_model(model2_id[0]).connected.append(model2_id)

    return True
