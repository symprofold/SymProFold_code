import bib
import ctl
import geometry


def get_trans_rot_param(ax0_fold, ax1_fold, model_number):
    '''
    Get rotation angle of ax1 when placed around central ax0.

    model_number: starting with 0
    '''
    rot_angle = 0

    if ax0_fold%ax1_fold == 0:
        # e.g. ax0:3, ax1:3, ax0:4, ax1:2, ax0:4, ax1:4, ax0:6, ax1:2
        rot_angle = ( model_number%(ax0_fold/ax1_fold) )*(360/ax0_fold)

    else:
        # e.g. ax0:3, ax1:2
        if ax0_fold == 3 and ax1_fold == 2:
            rot_angle = ( model_number%(ax0_fold) )*(360/ax0_fold)

    return rot_angle


def snapin_layer(axes, layer, conf, preserve_connections=False):
    '''
    Snap axes 0 and 1 to orientation points.
    2 steps:
    - step 1: alignment of first (initial) representation/representative of
        axis 0 to [0, 0, 0], with global relative translation of
        all representations/representatives of all axes
    - step 2:
        snap-in of axes 0 and 1 to snapin points (orientation points)
    '''
    max_snapin_distance = 80
    
    ax0 = axes[0]
    sess = ax0.chimerax_session
    intermediate_id = 101

    ax0_order0_model = layer.ax_models(ax0, order=0)[0]
    ax0_center = ax0_order0_model.get_center()


    # collect representations (models) to snapin for each axis
    models_to_snapin = [[] for ax in axes]

    for i,ax in enumerate(axes):
        models_to_snapin[i] = layer.ax_models(ax)


    # step 1: align center of ax0 to [0, 0, 0]
    # ----------------------------------------
    for i,ax in enumerate(axes):
        for j,m in enumerate(models_to_snapin[i]):
            sess.run('move x '+str(-ax0_center[0])+' models #'+m.idstr)
            sess.run('move y '+str(-ax0_center[1])+' models #'+m.idstr)


    # step 2: snap axis 0 and 1 representations to snapin points
    # (orientation points)
    # --------------------------------------------------------

    # iterate  through axes
    for i, ax in enumerate(axes):

        # axis 0
        if i == 0:
            # vector from snapin point (which is origin 0, 0, 0) to
            # axis center of rep1
            vect_start = ax0_order0_model.get_center()

        # axis 1
        if i == 1:
            # initial coords of representation 1 (rep1)
            ax1_rep1 = models_to_snapin[i][ax0.preferred_bs].get_center()

            # get snapin point for representation 1 (rep1) of ax1
            # Snapin point is nearest reference point.
            snapin_p_first, d = layer.get_ref_point_ax1(ax1_rep1, ax)
            d_xy = geometry.dist_xy(ax1_rep1, snapin_p_first)

            if d_xy > max_snapin_distance:
                ctl.e('snapin_distance')
                ctl.e(d_xy)
                ctl.e('ax1_rep1')
                ctl.e(ax1_rep1)
                ctl.e('snapin_p_first')
                ctl.e(snapin_p_first)
                raise Exception('snapin_layer: snapin distance too large')

            vect_start = [ax1_rep1[0]-snapin_p_first[0], \
                          ax1_rep1[1]-snapin_p_first[1], \
                          ax1_rep1[2]-snapin_p_first[2]]

        # further axes set
        if i >= 2:
            continue

        ctl.d('vect_start')
        ctl.d(vect_start)


        # for current axis: iterate through axis representations to snapin
        for j, model_to_snapin in enumerate(models_to_snapin[i]):
            ctl.d(model_to_snapin.id)

            # initial coords of axis representation (ax_rep)
            ax_rep = model_to_snapin.get_center()


            # axis 0
            # ~~~~~~
            # Center of each ax0 representation is directly aligned to
            # snapin point.

            if i == 0:
                # coord of real molecule center
                ax_rep_center = [ax_rep[0]-vect_start[0], \
                                 ax_rep[1]-vect_start[1], \
                                 ax_rep[2]-vect_start[2]]
                snapin_p, d = layer.get_ref_point(ax_rep_center)
                d_xy = geometry.dist_xy(ax_rep_center, snapin_p)

                if d_xy > max_snapin_distance:
                    ctl.e('model_to_snapin')
                    ctl.e(model_to_snapin.id)
                    ctl.e('snapin_distance')
                    ctl.e(d_xy)
                    ctl.e('ax_rep_center')
                    ctl.e(ax_rep_center)
                    ctl.e('snapin_p')
                    ctl.e(snapin_p)
                    raise Exception('snapin_layer: snapin distance too large')


            # axis 1
            # ~~~~~~
            # Center of each ax1 representation is directly aligned to
            # snapin point.

            if i == 1:
                # get snapin point for representation (ax_rep),
                # snapin point is nearest reference point
                snapin_p, d = layer.get_ref_point_ax1(ax_rep, ax)
                d_xy = geometry.dist_xy(ax_rep, snapin_p)

                if d_xy > max_snapin_distance:
                    ctl.e('model_to_snapin')
                    ctl.e(model_to_snapin.id)
                    ctl.e('snapin_distance')
                    ctl.e(d_xy)
                    ctl.e('ax_rep')
                    ctl.e(ax_rep)
                    ctl.e('snapin_p')
                    ctl.e(snapin_p)
                    raise Exception('snapin_layer: snapin distance too large')


                # align representation (ax_rep) to model with
                # preferred binding site
                if j != ax0.preferred_bs:
                    sess.match( \
                            (model_to_snapin.id[0], 1), (), \
                            (models_to_snapin[i][ax0.preferred_bs].id[0], 1), \
                            model_to_snapin.id)


                # calculated (snapin) rotation angle of for
                # representation (ax_rep) regarding preferred_bs model
                rot_angle = get_trans_rot_param(axes[0].fold, axes[1].fold, j)


                # perform relative rotation between models, if necessary and
                # update connections
                # Rotation is necessary if rot_angle != 0
                # (not necessary e.g. ax0 and ax1 are both 4-fold).
                if rot_angle != 0:

                    # perform snapin rotation of representation (ax_rep)
                    rot_center = model_to_snapin.get_center()

                    if ax0_order0_model.flipped == True:
                        sign = -1
                    else:
                        sign = 1

                    sess.run('turn z '+str(rot_angle*sign)+' center '+ \
                             str(rot_center[0])+', '+ \
                             str(rot_center[1])+', '+ \
                             str(rot_center[2])+ \
                             ' models #'+ \
                             model_to_snapin.idstr)


                # rot angle of representation (ax_rep) (before rotation by
                # rot_angle)
                ax_rep_rot_angle = model_to_snapin.rot_angle
                ctl.d('rot_angle')
                ctl.d(rot_angle)
                ctl.d('model_to_snapin.id')
                ctl.d(model_to_snapin.id)
                ctl.d('ax_rep_rot_angle')
                ctl.d(ax_rep_rot_angle)

                # rot angle of preferred_bs model
                rot_angle_preferred = \
                                models_to_snapin[i][ax0.preferred_bs].rot_angle
                ctl.d('rot_angle_preferred:')
                ctl.d(rot_angle_preferred)

                d_ang = rot_angle-(ax_rep_rot_angle-rot_angle_preferred)
                ctl.d('d_ang')
                ctl.d(d_ang)

                steps = d_ang/(360/axes[1].fold)
                ctl.d('steps')
                ctl.d(steps)

                steps_rounded = int(round(steps))%axes[1].fold
                        # modulo (%) to include as well e.g. +359°

                if steps_rounded < 0:
                    ctl.d('steps_rounded')
                    ctl.d(steps_rounded)
                    ctl.error('snapin_layer: steps_rounded < 0')
                elif steps_rounded > 2:
                    ctl.d('steps_rounded')
                    ctl.d(steps_rounded)
                    ctl.d(steps)
                    raise Exception('snapin_layer: '+ \
                                            'snapin rotation step too large')


                # update connections only when necessary
                if steps_rounded != 0 and not preserve_connections: 

                    if steps_rounded not in [1, 2]:
                        ctl.d('steps_rounded')
                        ctl.d(steps_rounded)
                        ctl.error('snapin_layer: steps_rounded error')

                    if abs((steps%axes[1].fold)-steps_rounded) > 0.2:
                        ctl.e('steps-steps_rounded too large')
                        ctl.e(steps-steps_rounded)
                        ctl.e('steps_rounded')
                        ctl.e(steps_rounded)
                        ctl.e('steps')
                        ctl.e(steps)
                        ctl.e(steps%axes[1].fold)
                        raise Exception('snapin_layer: '+ \
                                                'snapin rotation too large')

                    # previous connections
                    conn0 = model_to_snapin. \
                            get_connections_for_modes(['flattened', 'default'])

                    # updated connections
                    conn1 = {}
                    conn1_dest = {}

                    for c in conn0:
                        c_new = (c[0], \
                            ((c[1]-steps_rounded-1)%axes[1].fold)+1)
                        conn1[(c_new[0], c_new[1])] = conn0[c]

                        # create self-connections to avoid inhertited
                        # connections
                        if c not in conn1:
                            conn1[c] = c

                        if conn0[c] not in conn1_dest:
                            conn1_dest[conn0[c]] = conn0[c]

                    # create self-connections to avoid inhertited
                    # connections
                    for c in conn1:
                        model_to_snapin.set_connection(c, conn1[c], 'snapin')

                    for c in conn1_dest:
                        conn1_dest = sess.model_reg.get_model( \
                                (c[0],)).get_connections('snapin')

                        if c not in conn1_dest:
                            sess.model_reg.get_model((c[0],)). \
                                set_connection(c, conn1_dest[c], 'snapin')


                # coords (complex center) of representation (ax_rep) after
                # match to preferred_bs model and rotation
                ax_rep_center = model_to_snapin.get_center()


            # translation vector to final position on snapin point
            transl = [snapin_p[0]-ax_rep_center[0], \
                      snapin_p[1]-ax_rep_center[1]]

            # translate model to calculated snapin position
            sess.run('move x '+str(transl[0])+ \
                                    ' models #'+str(model_to_snapin.id[0]))
            sess.run('move y '+str(transl[1])+ \
                                    ' models #'+str(model_to_snapin.id[0]))


    models_all = [i for i in range(1, sess.last_id()+1)]
            # list of all models    

    chainids_new = sess.combine(models_all, intermediate_id)
    for model in models_all:
        sess.close_id(model)

    sess.split(intermediate_id)
    sess.run('save "'+conf.layer_snapin_raw_path+'"')

    return
