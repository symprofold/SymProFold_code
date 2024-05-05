import ctl
import data


'''
Module providing functions for determination and analysis of
domain binding sites.
'''

def get_interface_residue_ranges(distogram):
    '''
    Get interface residue ranges.
    '''
    residues = []

    for d in distogram:
        if d[0] not in residues:
            residues.append(d[0])

        if d[1] not in residues:
            residues.append(d[1])

    residues = sorted(residues)
    resranges = data.resids_to_ranges(residues)

    resranges_txt = ''

    for r in resranges:
        if r[0] != r[1]:
            resranges_txt += str(r[0])+'-'+str(r[1])+', '

        if r[0] == r[1]:
            resranges_txt += str(r[0])+', '

    resranges_txt = resranges_txt[0:-2]

    return resranges, resranges_txt


def get_domain(res_id, domain_boundaries):
    '''
    Get domain name of residue id.

    domain 0 for first domain
    '''
    for d_i,d in enumerate(domain_boundaries):
        if res_id in range(d[0], d[1]+1):

            return d_i

    ctl.e(res_id)
    ctl.e(domain_boundaries)
    ctl.error('get_domain')

    return


def get_interface_partner_resids(res_range, interface_distogram):
    '''
    Get partner residues of interfaces.
    '''
    partner_resids = []
    range_resids = data.ranges_to_resid_list([res_range])

    for d in interface_distogram:

        if d[0] in range_resids and interface_distogram[d][0] <= 8:
            if d[1] not in partner_resids:
                partner_resids.append(d[1])

        if d[1] in range_resids and interface_distogram[d][0] <= 8:
            if d[0] not in partner_resids:
                partner_resids.append(d[0])
        
    partner_resids = sorted(partner_resids)

    return partner_resids


def get_domain_bindings(domain_boundaries, interface_distogram):
    '''
    Get pairs of binding domains.
    '''
    interfaces = []
    partner_domains = {}

    for d_i, domain in enumerate(domain_boundaries):
        partner_resids = get_interface_partner_resids(domain, \
                                                      interface_distogram)

        for pr in partner_resids:
            partner_dom = get_domain(pr, domain_boundaries)

            if partner_dom not in partner_domains:
                key0 = min(d_i, partner_dom)
                key1 = max(d_i, partner_dom)
                key = (key0, key1)

                if key not in partner_domains:
                    partner_domains[key] = 1
                else:
                    partner_domains[key] += 1

    partner_domains = sort_domain_bindings(partner_domains)

    return partner_domains


def sort_domain_bindings(domain_bindings):
    '''
    Sort pairs of binding domains.
    '''
    domain_bindings_sort = []
    domain_bindings_tuples = []

    for d in domain_bindings:
        domain_bindings_tuples.append((d[0], d[1], domain_bindings[d]))

    domain_bindings_sort = sorted(domain_bindings_tuples, \
                                  key=lambda x: x[2], reverse=True)

    return domain_bindings_sort
