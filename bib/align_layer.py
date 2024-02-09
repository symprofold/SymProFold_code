import bib
import bibpdb
import ctl
import geometry

import math
import os


def align_layer(axes, layer, conf, sess):
    '''
    Align layer to orientation points and move all other models relative to it.
    Typically the orientation points are in the xy plane.
    '''
    ax0 = axes[0]
    models_to_align = [i for i in range(1, sess.last_id()+1)]
            # list of all models    

    ax0_models = layer.ax_models(ax0)
            # use ax0_models as reference models for alignment

    ax_reprs_first = layer.get_representations()[ \
                                    next(iter(layer.get_representations()))]

    center_res = ax_reprs_first.termini[0]+ \
            round( (ax_reprs_first.termini[1]-ax_reprs_first.termini[0])/2 )

    sess.run('open "'+str(conf.layer_raw_path)+'"')

    # combine all models to intermediate_id
    intermediate_id = 101
    chainids_new = sess.combine(models_to_align, intermediate_id)

    chainids_A = []
    for c in chainids_new:
        if 'A' in c:
            for m in ax0_models: # use ax0_models as reference models for
                                 # alignment
                if 'mid'+str(m.id[0])+'point' in c:
                    chainids_A.append(c)

    align_str = ''
    for c in chainids_A:
        align_str = align_str+"#101//chain_id='"+c+"':"+str(center_res)+'@CA '

    for m in models_to_align:
        sess.close_id((m,))

    ctl.d(layer.lattice_constant)
    ctl.d(layer.symmgroup)
    ctl.d(layer.symmgroups_compatible)
    
    points = layer.orient_points()
    ctl.d(points)

    # set marker on orientation points
    marker_start_id = 300
    sess.set_marker(points, marker_start_id)
    sess.run('rename #'+str(marker_start_id)+'-'+str(marker_start_id+50)+ \
             ' id #'+str(marker_start_id)) # rename markers


    # before performing the alignment:
    # determine residue with lowest and the highest coordinate in z direction
    sess.run('save "'+str(conf.layer_aligned_raw_path[:-4])+'.pdb" #'+ \
             str(intermediate_id))
    coords = bibpdb.open_pdb(str(conf.layer_aligned_raw_path[:-4])+'.pdb', \
                            False)
    os.unlink(conf.layer_aligned_raw_path[:-4]+'.pdb')
    pre_alignment_extrema = geometry.get_extrema(coords, 'z', 'all')


    # align layer to orientation points
    sess.run('align '+align_str+' to #'+str(marker_start_id))

    # split combine models and assign original ids
    sess.split(intermediate_id)


    # after the alignment has been performed:
    # determine coords of residues that had the lowest and the highest
    # coordinate in z direction before the alignment
    pre_alignment_extremum_lowest_coord = \
                    sess.get_xyz((1, 1), \
                                 pre_alignment_extrema[0][1], \
                                 error_when_res_not_found=False)

    if pre_alignment_extremum_lowest_coord == False:
        pre_alignment_extremum_lowest_coord = \
                     sess.get_xyz((2, 1), \
                                 pre_alignment_extrema[0][1], \
                                 error_when_res_not_found=False)


    pre_alignment_extremum_highest_coord = \
                    sess.get_xyz((1, 1), \
                                 pre_alignment_extrema[1][1], \
                                 error_when_res_not_found=False)

    if pre_alignment_extremum_highest_coord == False:
        pre_alignment_extremum_highest_coord = \
                     sess.get_xyz((2, 1), \
                                  pre_alignment_extrema[1][1], \
                                  error_when_res_not_found=False)


    if pre_alignment_extremum_lowest_coord == False:
        ctl.error('align_layer: lowest residue not found')

    if pre_alignment_extremum_highest_coord == False:
        ctl.error('align_layer: highest residue not found')


    # assign as flipped if layer model was flipped during alignment
    if pre_alignment_extremum_highest_coord[2] < \
        pre_alignment_extremum_lowest_coord[2]:

        layer.ax_models(ax0, order=0)[0].flipped = \
                            not layer.ax_models(ax0, order=0)[0].flipped


    # turn probable outside of layer model to top (estimation)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sess.run('save "'+str(conf.layer_aligned_raw_path[:-4])+'.pdb" #1')
    coord = bibpdb.open_pdb(str(conf.layer_aligned_raw_path[:-4])+'.pdb', \
                            False)
    os.unlink(conf.layer_aligned_raw_path[:-4]+'.pdb')

    cen = geometry.get_center(coord, 'all')

    # before turn: combine all models to intermediate_id
    chainids_new = sess.combine(models_to_align, intermediate_id)
    for m in models_to_align:
        sess.close_id((m,))

    if cen[2] > 0:
        if layer.ax_models(ax0, order=0)[0].flipped:
            sess.run('turn x 180 models #'+str(intermediate_id))
            layer.ax_models(ax0, order=0)[0].flipped = \
                                not layer.ax_models(ax0, order=0)[0].flipped
            layer.flipped = not layer.flipped


    # split combine models and assign original ids
    sess.split(intermediate_id)

    sess.run('save "'+str(conf.layer_aligned_raw_path)+'"')

    return
