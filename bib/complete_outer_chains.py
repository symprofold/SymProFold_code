import bib
import ctl
import export


def complete_outer_chains(assembler, contact_submodel_orientation):
    ''' Complete outer chains of layer model. '''

    ax0 = assembler.axes[0]
    ax1 = assembler.axes[1]

    ax0_incomplete_models = []

    for r in ax0.get_representations():
        if r[0] > 1+2*ax0.fold:

            for i_id in ax0.chimerax_session.get_submodel_ids(r):
                i = ax0.chimerax_session.model_reg. \
                    convert_model_id_to_str(i_id)
                
                conn_ = ax0.get_representation(r).get_connection(i, 'snapin')
                if conn_ != False:
                    conn = conn_

                # check if models already merged
                elif i_id not in ax0.chimerax_session. \
                        model_reg.get_model(i_id[0]).connected:
                    ax0_incomplete_models.append(i)
                    ctl.d('incomplete model ax0: '+i)                
                    ctl.d(ax0_incomplete_models)

    ctl.d(ax0_incomplete_models)

    # value set to allow for an acceptable level of small clashes
    max_clashes = 5000

    for matchto in ax0_incomplete_models:
        contact_model21 = assembler.ax1_contact_submodel_0( \
                            int(matchto.split('.')[1]), ax1.fold)

        last_model = ax1.open_model(0, False)
        bib.format_model(last_model.id, ax1)

        assembler.ax_a_model_add_ax_b(ax0, ax1, \
                                      matchto,
                                      last_model.idstr+'.'+ \
                                      str(contact_model21), 'completechains', \
                                      False)
                # param "False", because connection not permanent/universal

        for i in range(1, ax1.fold+1):
            if i != contact_model21:
                ax0.chimerax_session.close_id((last_model.id[0], i))
                ctl.d('close')
                ctl.d((last_model.id[0], i))
               

        ax0.chimerax_session.run('hide #'+str(matchto)+' models')

        # hide ax1 of unflattened layer
        ax0.chimerax_session.run('hide #'+str(2+(1+0)*ax0.fold)+'-'+ \
                                          str(1+(2+0)*ax0.fold)+' models')

        clashes = ax0.chimerax_session.run('clashes #'+ \
                            last_model.idstr+'.'+str(contact_model21)+ \
                            ' restrict cross ignoreHiddenModels true'+'')
        ax0.chimerax_session.close_id(last_model.id[0]+1)
        ax0.chimerax_session.run('show #'+str(matchto)+' models')
        ax0.chimerax_session.run('show #'+str(2+(1+0)*ax0.fold)+'-'+ \
                                          str(1+(2+0)*ax0.fold)+' models')


        if len(clashes) > max_clashes:
            ctl.d('close due to clashes')
            ax0.chimerax_session.close_id((last_model.id[0], contact_model21))
            ctl.d('complete_outer_chains: close due to clashes, pos. 1')



    ax1_incomplete_models = []

    for r in ax1.get_representations():
        if 1+0*ax0.fold < r[0] <= 1+1*ax0.fold:

            for i_id in ax0.chimerax_session.get_submodel_ids(r):
                i = ax0.chimerax_session.model_reg. \
                    convert_model_id_to_str(i_id)

                conn_ = ax1.get_representation(r).get_connection(i, 'snapin')
                if conn_ != False:
                    conn = conn_

                # check if models already merged
                elif i_id not in ax0.chimerax_session. \
                         model_reg.get_model(i_id[0]).connected:
                    ctl.d(i_id)
                    ax1_incomplete_models.append(i)
                    ctl.d('incomplete model ax1: '+i)  


    ctl.d('ax0_incomplete_models')
    ctl.d(ax0_incomplete_models)
    ctl.d('ax1_incomplete_models')
    ctl.d(ax1_incomplete_models)



    for matchto in ax1_incomplete_models:
        contact_model12 = assembler.ax1_contact_submodel_1(
                                int(matchto.split('.')[1]), \
                                contact_submodel_orientation, ax1.fold)

        last_model = ax0.open_model(0, False)
        bib.format_model(last_model.id, ax0)

        assembler.ax_a_model_add_ax_b(ax1, ax0, \
                                    matchto, \
                                    last_model.idstr+'.'+ \
                                    str(contact_model12), 'completechains', \
                                    False)
                # param "False", because connection not permanent/universal

        for i in range(1, ax1.fold+1):
            if i != contact_model12:
                ax0.chimerax_session.close_id((last_model.id[0], i))
           

        ax0.chimerax_session.run('hide #'+str(matchto)+' models')

        # hide ax1 of unflattened layer
        ax0.chimerax_session.run('hide #'+str(2+(1+0)*ax0.fold)+'-'+ \
                                 str(1+(2+0)*ax0.fold)+' models')

        clashes = ax0.chimerax_session.run('clashes #'+ \
                                 last_model.idstr+'.'+str(contact_model12)+ \
                                 ' restrict cross ignoreHiddenModels true'+'')

        ax0.chimerax_session.close_id(last_model.id[0]+1)
        ax0.chimerax_session.run('show #'+str(matchto)+' models')
        ax0.chimerax_session.run('show #'+str(2+(1+0)*ax0.fold)+'-'+ \
                                 str(1+(2+0)*ax0.fold)+' models')


        if len(clashes) > max_clashes:
            ctl.d('close due to clashes')
            ax0.chimerax_session.close_id((last_model.id[0], contact_model12))
            ctl.d('complete_outer_chains: close due to clashes, pos 2')

    ax0.chimerax_session.run('save "'+ \
                             assembler.conf.layer_complete_chains_raw_path+'"')

    return
