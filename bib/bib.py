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


def open_model(sess, file, model_id, meta={}, part=0, \
               termini_with_signalsequence=True):
    '''
    Open model and calculate basic parameter.
    '''
    if part == 0 or part == 1:
        sess.open_model(file)
        model_id = sess.last_id()
    if part == 0 or part == 2:
        sess.hide_atoms([model_id])
        sess.show_cartoons(model_id)
        sess.split_model(model_id)
    if part == 0 or part == 1:
        rmsds = bibpdb.get_rmsds(file)
        termini = molmodel.get_termini( \
                    rmsds, \
                    termini_with_signalsequence=termini_with_signalsequence)

        multimer_n = get_multimer_n(file)

        meta[model_id] = [[], termini, multimer_n, [[], []]]
      
    return model_id, meta


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

    session.set_marker([[0, 0, 0]], (100,), radius=1)
    session.set_marker([[100, 0, 0]], (101,), radius=1)
    session.set_marker([[0, 100, 0]], (102,), radius=1)

    return


def model_2fold_to_plane(model_id, meta, session):
    ''' Align 2-fold axis model to xy plane. '''

    model = Model(session.model_reg)
    model.set_id(model_id)

    tmp_dir = filesystem.clean_path( \
            os.path.dirname(os.path.realpath(__file__))+'/../../')
    tmpfile = tmp_dir+'tmp.pdb'
    session.save_models([model], tmpfile, 'pdb')

    tmp_model_id = session.open_model(tmpfile)
    tmp_model = Model(session.model_reg)
    tmp_model.set_id(tmp_model_id)
    session.split_model(tmp_model.id)
    session.match((tmp_model.id[0], 2), [], (tmp_model.id[0], 1))

    rot_axis = session.measure_rot_axis((tmp_model.id[0], 2), \
                                        (tmp_model.id[0], 1))
    session.close_id(tmp_model.id)

    ctl.d(rot_axis)

    center_res = round((meta[model_id][1][0]+ \
                        meta[model_id][1][1])/2)

    center0 = session.get_xyz(model.id, center_res, chainid='A')
    center1 = session.get_xyz(model.id, center_res, chainid='B')

    c0c1 = [center1[0]-center0[0], center1[1]-center0[1], \
            center1[2]-center0[2]]

    normal = geometry.crossproduct(c0c1, rot_axis)
    ctl.d(center_res)
    ctl.d(normal)

    main_dir = filesystem.clean_path( \
            os.path.dirname(os.path.realpath(__file__))+'/../')

    marker_id = session.open_model(main_dir+'/component_files/marker.pdb')
    marker = Model(session.model_reg)
    marker.set_id(marker_id)
    marker.set_session(session)

    marker.move_model([center0[0]+normal[0], \
                       center0[1]+normal[1], \
                       center0[2]+normal[2]])

    session.combine_models([model, marker], str(marker.id[0]+1))
    session.set_marker([normal], (201,), color='yellow', radius=1)

    session.align_model([((marker.id[0]+1,), 'A', center_res), \
                         ((marker.id[0]+1,), 'B', center_res), \
                         ((marker.id[0]+1,), 'C', 1)], \
                        [(100,), (101,), (102,)])

    session.close_id(tmp_model.id)

    session.split_model((tmp_model.id[0]+1,))
    session.close_id((tmp_model.id[0]+1, 3))

    ctl.d(model.id)
    ctl.d(str(model.id[0]+1))

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
        session.align_model([(model_id, 'A', center_res), \
                             (model_id, 'B', center_res), \
                             (model_id, 'C', center_res)], \
                            [(100,), (101,), (102,)])

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
