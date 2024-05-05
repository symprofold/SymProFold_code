import ctl
import geometry
import math


'''
Module providing functions for handling and analyzing molecular models.
'''

def get_sibling_rmsd(model_id0, model_id1, resranges, sess):
    '''
    Calculate rmsd of the distances between all residues of a given model and
    a given sibling model.
    '''
    distances = []
    sibling_distances = {}

    resids0 = sess.resids(model_id0)
    resids1 = sess.resids(model_id1)

    if resids0 != resids1:
        ctl.error('get_mate_rmsd: given models are no siblings.')

    for resid in resids0:
        coords0 = sess.get_xyz(model_id0, resid)
        coords1 = sess.get_xyz(model_id1, resid)

        d = geometry.dist(coords0, coords1)
        distances.append(d)
        sibling_distances[resid] = d

    distances_squared = [d**2 for d in distances]
    rmsd = math.sqrt(sum(distances_squared)/len(distances))

    return rmsd, sibling_distances


def get_termini(rmsds, start_res=1, termini_with_signalsequence=True):
    '''
    Get termini from given RMSD list.

    Args:
        termini_with_signalsequence: signal included in N terminus
    '''
    cutoff_N = 5 # cutoff in Angstrom
    cutoff_C = cutoff_N

    termini = [-1,10000]
    for r in rmsds:
        if r <= start_res:
            continue
        if termini[0] == -1 and rmsds[r] <= cutoff_N:
            termini[0] = r-1
        if rmsds[r] <= cutoff_C:
            termini[1] = r+1

    if termini_with_signalsequence:

        # signal sequence included in N terminus
        if 21 in rmsds and 25 in rmsds:
            if termini[0] <= 20 and \
               rmsds[21] > cutoff_N and rmsds[22] > cutoff_N and \
               rmsds[23] > cutoff_N and rmsds[24] > cutoff_N and \
               rmsds[25] > cutoff_N:
                return get_termini(rmsds, start_res=21)

    return termini


def interface_residues_per_model(interface_residues):
    '''
    Create a dict with a list of interface residues for each model_id.
    '''
    residues_per_model = {}

    for r in interface_residues:
        if r[0] not in residues_per_model:
            residues_per_model[r[0]] = [r[1]]
        else:
            if r[1] not in residues_per_model[r[0]]:
                residues_per_model[r[0]].append(r[1])
                residues_per_model[r[0]] = sorted(residues_per_model[r[0]])

    return residues_per_model
