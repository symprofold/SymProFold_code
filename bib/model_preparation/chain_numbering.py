import bibpdb
import ctl
import geometry
import model_preparation.geometry


'''
Module providing functions to handle chain numbers.
'''

def realign_chain_ids(sensing_file, input_file, output_file):
    '''
    Rename all chain ids into alphabetical order.
    E.g. if first chain id is "B", rename all chain ids to start with "A".
    (Preprocessing of pdb input file.)
    '''

    # check first chainid in sensing file
    # separate sensing file to avoid problem with header (ChimeraX export)
    coord = bibpdb.open_pdb(sensing_file, False)

    chain_id_names = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    chainid_order_numeric = []

    for c in coord:
        chainid = bibpdb.get_f(c, 'chainid')      

        # get position in alphabet
        al_pos = chain_id_names.find(chainid)

        if al_pos not in chainid_order_numeric:
            chainid_order_numeric.append(al_pos)


    # reorder only if neccessary
    for i,o in enumerate(chainid_order_numeric):
        if i != o:
            bibpdb.reorder_chainids(input_file, chainid_order_numeric, \
                                    output_file)
            return

    return


def determine_subchain_order(coord, meta, mult):
    '''
    Determine counterclockwise order of subchain ids.
    '''
    chain_id_names = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    chaincenter = []

    for i in range(0, mult):
        chaincenter.append(geometry.get_center(coord, chain_id_names[i], \
                           meta[1]))

    ctl.d('chaincenter')
    ctl.d(chaincenter)

    order = [0]

    for m in range(1, mult):
        for i,v in enumerate(chaincenter):

            crossp = geometry.crossproduct(chaincenter[order[m-1]], \
                                           chaincenter[i])
            if -0.1 < crossp[2] < 0.1:
                ctl.d('continue, because cp ~= 0')
                ctl.d(crossp[2])
                continue

            if crossp[2] < 0:
                ctl.d('continue, because cp > 0')
                ctl.d(crossp[2])
                continue

            # calculate distances of v to all chaincenters
            centerDist = model_preparation.geometry. \
                                        get_distances(chaincenter, v)
            minDist = min(centerDist)
            currentDist = centerDist[order[m-1]]
            if currentDist <= minDist*1.2:
                order.append(i)

            ctl.d('pointer:'+str(i))
            ctl.d('crossp')
            ctl.d(crossp)
            ctl.d('centerDist')
            ctl.d(centerDist)
            ctl.d('minDist')
            ctl.d(minDist)
            ctl.d('currentDist')
            ctl.d(currentDist)
            ctl.d('order')
            ctl.d(order)

        # remove duplicate values
        order2 = []

        for o in order:
             if o not in order2:
                 order2.append(o)

        order = order2

        if len(order) < m+1:
           ctl.d(order)
           ctl.d(m)
           ctl.d('break, len(order) < m+1')
           return -1
           break

    return order 
