import ctl

import chimerax


def get_center(model_id, submodel_n, termini, sess):
    ''' Get center of complex model. '''

    model_id = sess.model_reg.convert_model_id(model_id)

    center = [0, 0, 0]
    sum_vect = 0
    sum_n = 0        

    for submodel_id in range(1, submodel_n+1):

        resids = sess.resids((model_id[0], 1))

        for resid in resids:
            if termini[0] <= resid <= termini[1]:
                sum_vect = sum_vect+ \
                        sess.get_xyz((model_id[0], submodel_id), int(resid))
                sum_n += 1

    center = sum_vect/sum_n

    return center


def get_binding_status(model_id, model_id_res_range, res_range, \
                       lddt_min, sess):
    '''
    Checks if residue range res_range is involved in inter submodel contacts of
    complex model_id with residue range model_id_res_range.
    '''

    # parameter model_id_res_range='all' stands for all residues
    if model_id_res_range == 'all':
        model_id_res_range = [1, 10000]
    
    sess.run('split #'+str(model_id[0]))
    bs = []

    try:
        co = sess.run('contacts #'+str(model_id[0])+':'+ \
                  str(model_id_res_range[0])+'-'+str(model_id_res_range[1])+ \
                  '@@bfactor>='+str(lddt_min)+ \
                  ' '+ \
                  'restrict #'+str(model_id[0])+':'+ \
                  str(res_range[0])+'-'+str(res_range[1])+ \
                  '@@bfactor>='+str(lddt_min)+ \
                  ' intraModel false interModel true interSubmodel true')
    except chimerax.core.errors.UserError:
        co = []

    if len(co) >= 10: # parameter set to 10, can be adjusted
        return 1, len(co)

    return 0, len(co)


def get_binding_status_domains(model_id, model_id_res_range, domains, \
                               lddt_min, sess):
    '''
    Checks which residue ranges in the list domains are involved in
    inter submodel contacts of complex model_id.
    '''
    dom_binding = []
    binding_domains = []

    for d in domains:
        
        bs, co = get_binding_status(model_id, model_id_res_range, d, \
                                    lddt_min, sess)
        dom_binding.append([bs, co])

        if bs == 1:
            binding_domains.append(d)

    return dom_binding, binding_domains


def get_main_binding_domain(model_id, model_id_res_range, domains, \
                            lddt_min, sess):
    '''
    Get domain with most contacts to other parts of the model.

    max_dom: domain index (starting with 0) with most contacts.
    '''
    main_binding_domain = []

    dom_binding, binding_domains = get_binding_status_domains( \
            model_id, model_id_res_range, domains, lddt_min, sess)

    max_dom = -1
    max_cl = 0

    for i,db in enumerate(dom_binding):
        if db[0] == 1 and db[1] > max_cl:
            max_dom = i
            max_cl = db[1]

    if max_dom != -1:
        main_binding_domain = [max_dom, max_cl]

    return main_binding_domain, dom_binding
