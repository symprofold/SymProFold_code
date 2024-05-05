import ctl


'''
Module providing functions for model preparation to determine clashes.
'''

def clashes_per_100(model_id, sess):
    '''
    Determine clashes per 100 residues and return a filename infix with this
    information.
    '''
    clashes = sess.run('clashes #'+str(model_id)+' & ~:cys', sess)
            # skip cysteines to avoid counting clashes caused by unrecognized
            # disulfide bonds
            #
            # notice: a ChimeraX model including the clashes is created, which
            # is deleted below with the command "close"

    clashes_n = len(clashes)
    res_n = sess.get_res_n((model_id,))

    # calculate clashes per 100 residues
    cla_per_100 = round(clashes_n/res_n*100, 2)
    cla_per_100_txt = str(cla_per_100).split('.')
    cla_per_100_txt = cla_per_100_txt[0]+'-'+(cla_per_100_txt[1]+'00')[:2]
    clash_infix = '_cla'+cla_per_100_txt
    sess.run('close #'+str(model_id+1))

    return clash_infix, cla_per_100


def clashes_per_100_rolling(model_id, rolling_len, mult, sess):
    '''
    Determine clashes per 100aa in a rolling range along the whole chain.
    Return a filename infix with this information.
    '''
    clashes = sess.run('clashes #'+str(model_id)+' & ~:cys')
            # skip cysteines to avoid counting clashes caused by unrecognized
            # disulfide bonds
            #
            # notice: a ChimeraX model including the clashes is created, which
            # is deleted below with the command "close"

    resids = sess.resids(model_id)
    res_clashes = [[resid, 0] for resid in range(0, resids[-1]+1)]

    for cl in clashes:
        fields = str(cl).split(' ')
        resid = int(fields[3])
        res_clashes[resid][1] += 1

    rolling_max = 0
    rolling_max_centerresid = 0

    if resids[-1]-resids[0]+1 < rolling_len:
        zero, rolling_max_per100 = clashes_per_100(model_id, sess)
        rolling_max_centerresid = resids[0]+round((resids[-1]-resids[0])/2)

    else:
        for r_start in range(1, len(res_clashes)-rolling_len):
            roll_vals = [v[1] for v in \
                         res_clashes[r_start:(r_start+rolling_len)]]
            roll_sum = sum(roll_vals)

            if roll_sum > rolling_max:
                rolling_max = roll_sum
                rolling_max_centerresid = round(r_start+rolling_len/2)

        rolling_max = rolling_max/mult
            # mult: number of chains (oligomeric state) of model_id

        rolling_max_per100 = round(rolling_max/rolling_len*100, 2)

    roll_cla_per_100_txt = str(rolling_max_per100).split('.')
    roll_cla_per_100_txt = roll_cla_per_100_txt[0]+'-'+ \
                               (roll_cla_per_100_txt[1]+'00')[:2]
    roll_clash_infix = '_rollcla'+roll_cla_per_100_txt+'_'+ \
                               str(rolling_max_centerresid)
    sess.run('close #'+str(model_id+1))

    return roll_clash_infix, rolling_max_per100
