import bib
import ctl
import geometry


def complex_merge(id1, ids2_, sess):
    ''' Merge molecule complex 1 with list of ids (ids2). '''

    joinwith = []

    # return [] if id1 does not exist or no submodels do exist
    if not sess.model_reg.model_exists(id1):
        return []

    sids1 = sess.get_submodel_ids(id1)

    if sids1 == []:
        return []


    # skip not existing model ids in ids2_
    ids2 = []
    for id2 in ids2_:
        if sess.model_reg.model_exists(id2):
            ids2.append(id2)
        
    if ids2 == []:
        return []


    for subid1 in sids1:
        ctl.d('subid1_'+str(subid1))

        join_candidates = protein_merge(subid1, ids2, sess)

        if join_candidates != []:

            # ensure that each partner has the other partner als
            # closest candidate
            reverse_check = protein_merge(join_candidates[0][1], [id1], sess)
            ctl.d('reverse_check')
            ctl.d(join_candidates)
            ctl.d(reverse_check)
            if join_candidates[0][0] == reverse_check[0][1] and \
               join_candidates[0][1] == reverse_check[0][0]:
                joinwith.append(join_candidates[0])

                ctl.d('reverse check ok')
            else:    
                ctl.d('reverse check failed')

    return joinwith


def protein_merge(subid1, ids2, sess):
    ''' Merge molecule 1 with list of ids (ids2). '''

    joinwith = []
    subid1 = sess.model_reg.convert_model_id_to_str(subid1)

    ctl.d('subid1_'+str(subid1))

    for id2 in ids2:
        ctl.d('id2_'+str(id2))

        sids2 = sess.get_submodel_ids(id2)

        for subid2 in sids2:

            dist_cutoff = 40
                    # dist_cutoff: to adjust (higher) when chains
                    # do not get connected
            subid2 = sess.model_reg.convert_model_id_to_str(subid2)
            
            cm, dist = check_merge_possible(subid1, subid2, sess, dist_cutoff)
            ctl.d(cm)

            if cm == 1:
                ctl.d('merge possible')
                ctl.d(subid1)
                joinwith.append([subid1, subid2, dist])
                ctl.d(joinwith)
            else:
                ctl.d('merge not possible')
                ctl.d(subid1)

    if joinwith == []:
        return joinwith
    else:
        # get candidate with lowest distance
        smallest_dist = joinwith[0]
        for j in joinwith:
            if j[2] < smallest_dist[2]:
                smallest_dist = j

        return [smallest_dist]


def check_merge_possible(id_a, id_b, sess, dist_cutoff=40):
    '''
    Check whether merging of 2 submodels is possible regarding distance.
    '''
    ctl.d(id_a)
    ctl.d(id_b)

    resids_a, resid_fragments_a = sess.resids_surface_fragments(id_a)
    resids_b, resid_fragments_b = sess.resids_surface_fragments(id_b)

    ctl.d(resids_a)
    ctl.d(resid_fragments_a)
    ctl.d(resids_b)
    ctl.d(resid_fragments_b)

    d = 1000

    if resids_a == [] or resids_b == []:
        ctl.e(id_a)
        ctl.e(id_b)
        ctl.e(resids_a)
        ctl.e(resids_b)
        ctl.error('check_merge_possible2: resids_a or resids_b empty')
        return 0, 1000


    ctl.d(resids_a)
    ctl.d(resids_b)
    ctl.d(resid_fragments_a)
    ctl.d(resid_fragments_b)

    id1 = id_a
    id2 = id_b
    resids1 = resids_a
    resids2 = resids_b
    resid_fragments1 = resid_fragments_a
    resid_fragments2 = resid_fragments_b

    if resids_a[0] > resids_b[0]:
        id1 = id_b
        id2 = id_a
        resids1 = resids_b
        resids2 = resids_a
        resid_fragments1 = resid_fragments_b
        resid_fragments2 = resid_fragments_a


    if len(resid_fragments1) == 1 and len(resid_fragments2) == 1:
        if abs(resids1[-1]-resids2[0]) <= 1:
            c1 = sess.get_xyz(id1, resids1[-1])
            c2 = sess.get_xyz(id2, resids2[0])
            d = geometry.dist(c1, c2)
            if d <= dist_cutoff:
                return 1, d


    if len(resid_fragments1) == 2 and len(resid_fragments2) == 1:

        if abs(resid_fragments1[0][-1]-resids2[0] )<=1 and \
           abs(resid_fragments1[1][0] -resids2[-1])<=1:
            c1 = sess.get_xyz(id1, resid_fragments1[0][-1])
            c2 = sess.get_xyz(id2, resids2[0])

            d1 = geometry.dist(c1, c2)

            c3 = sess.get_xyz(id1, resid_fragments1[1][0])
            c4 = sess.get_xyz(id2, resids2[-1])            
            d2 = geometry.dist(c3, c4)

            d = (d1+d2)/2

            if d1 <= dist_cutoff and d2 <= dist_cutoff:
                return 1, d

    return 0, d
