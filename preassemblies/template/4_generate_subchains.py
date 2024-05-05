import glob
import os
import lib
import pathlib
import sys


def seqlen_dom_range(domains, domain_start, domain_end):
    '''
    Calculate sequence length of subchain defined by start domain and
    end domain.
    '''
    seqlen = 0

    for d in range(domain_start-1, domain_end):
        seqlen += len(domains[d])

    return seqlen


def write_subchains_to_fasta(subchains, domains, multiplicities, path_fasta):
    '''
    Create fasta files containing the given subchain with each of the given
    multiplicities.
    '''
    for sc in subchains:

        # generate each of the given multiplicity (e.g. x2, x3, ...)
        for m in multiplicities:
            if sc[0] == 1 and sc[1] == dom_n:
                filename = path_fasta[:-4]+'x'+str(m)+'.fa'
            elif sc[0] == sc[1]:
                filename = path_fasta[:-4]+str(sc[0])+'x'+str(m)+'.fa'
            else:    
                filename = \
                        path_fasta[:-4]+str(sc[0])+str(sc[1])+'x'+str(m)+'.fa'

            f = open(filename, 'w')

            for i in range(0, m):
                f.write(header[0]+"\r\n")
                
                for d in range(sc[0]-1, sc[1]):
                    f.write(domains[d])

            f.close()

    return


def write_subchains_to_single_fasta(subchains, domains, path_fasta):
    '''
    Create fasta file containing all subchains to use as input for hetero msa
    generation.
    '''
    subchains_sort = sorted(enumerate(subchains[1:]), key=lambda x: x[1])
    subchains_sort = [subchains[0]]+[s[1] for s in subchains_sort]

    filename = path_fasta[:-4]+'x1'
    for sc in subchains_sort[1:]:
        if sc[0] == sc[1]:
            filename += '_'+str(sc[0])
        else:
            filename += '_'+str(sc[0])+str(sc[1])

    filename += '.fa'
    f = open(filename, 'w')

    for sc in subchains_sort:
        f.write(header[0]+"\r\n")
        
        for d in range(sc[0]-1, sc[1]):
            f.write(domains[d])

    f.close()

    return


def generate_subchains(domains, sequ):
    '''
    Generate standard set of subchains.
    '''
    subchains = []
    dom_n = len(domains)

    # full length
    # ~~~~~~~~~~
    sc = [1, dom_n, 'FL']
    subchains.append(sc)


    # w/o C term
    # ~~~~~~~~~~

    # ensure that >= 10% is removed
    fraction_removed_1_dom = len(domains[dom_n-1])/len(sequ)

    if fraction_removed_1_dom >= 0.1:
        sc = [1, dom_n-1, 'w/o_c']
        if sc not in subchains:
            subchains.append(sc)

    elif len(domains) >= 2:
        fraction_removed_2_dom = \
                (len(domains[dom_n-2])+len(domains[dom_n-1]))/len(sequ)
        if fraction_removed_2_dom >= 0.1:
            sc = [1, dom_n-2, 'w/o_c']
            if sc not in subchains:
                subchains.append(sc)

    elif len(domains) >= 3:
        fraction_removed_3_dom = \
                (len(domains[dom_n-3])+len(domains[dom_n-2])+ \
                    len(domains[dom_n-1]))/len(sequ)
        if fraction_removed_3_dom >= 0.1:
            sc = [1, dom_n-3, 'w/o_c']
            if sc not in subchains:
                subchains.append(sc)

    else:
        print('ERROR: generate_subchains: '+ \
                            'w/o C term: 10% removal not reached')
        sys.exit()
        

    # w/o N term
    # ~~~~~~~~~~
    if len(domains) > 1:

        # ensure that min. 10% is removed
        fraction_removed_1_dom = len(domains[0])/len(sequ)

        if fraction_removed_1_dom >= 0.1:
            sc = [2, dom_n, 'w/o_n']
            if sc not in subchains:
                subchains.append(sc)

        elif len(domains) >= 2:
            fraction_removed_2_dom = \
                    (len(domains[0])+len(domains[1]))/len(sequ)
            sc = [3, dom_n, 'w/o_n']
            if sc not in subchains:
                subchains.append(sc)

        elif len(domains) >= 3:
            fraction_removed_3_dom = \
                    (len(domains[0])+len(domains[1])+len(domains[2]))/len(sequ)
            sc = [4, dom_n, 'w/o_n']
            if sc not in subchains:
                subchains.append(sc)

        else:
            print('ERROR: generate_subchains: '+ \
                            'w/o N term: 10% removal not reached')
            sys.exit()

        
    # first third of domains
    # ~~~~~~~~~~~~~~~~~~~~~~
    fraction = seqlen_dom_range(domains, 1, round(dom_n/3))/len(sequ)
    fraction_plus_1_dom = \
                    seqlen_dom_range(domains, 1, round(dom_n/3)+1)/len(sequ)

    if fraction >= 0.2:
        sc = [1, round(dom_n/3), 'first_1/3']
        if sc not in subchains:
            subchains.append(sc)
    elif fraction_plus_1_dom >= 0.2:
        sc = [1, round(dom_n/3)+1, 'first_1/3']
        if sc not in subchains:
            subchains.append(sc)
    else:
        print('ERROR: generate_subchains: first third: 20% not reached')
        sys.exit()


    # last third of domains
    # ~~~~~~~~~~~~~~~~~~~~~
    fraction = seqlen_dom_range(domains, \
                                round(dom_n-dom_n/3)+1, dom_n)/len(sequ)
    fraction_plus_1_dom = seqlen_dom_range(domains, \
                                round(dom_n-dom_n/3), dom_n)/len(sequ)

    if fraction >= 0.2:
        sc = [round(dom_n-dom_n/3)+1, dom_n, 'last_1/3']
        if sc not in subchains:
            subchains.append(sc)
    elif fraction_plus_1_dom >= 0.2:
        sc = [round(dom_n-dom_n/3), dom_n, 'last_1/3']
        if sc not in subchains:
            subchains.append(sc)
    else:
        print('ERROR: generate_subchains: last third: 20% not reached')
        sys.exit()


    subchains_unique = []

    for sc in subchains:
        if (sc[0], sc[1]) not in subchains_unique:
            subchains_unique.append((sc[0], sc[1]))

    return subchains_unique, subchains


def generate_subchain_ensemble(domains):
    '''
    Create full ensemble of all possible subchains.
    '''
    subchains = []
    dom_n = len(domains)

    for p0 in range(1, dom_n+1):
        for p1 in range(1, dom_n+1):
            if p1 >= p0:
                sc = [p0, p1]

                if sc not in subchains:
                    subchains.append(sc)

    return subchains


def prediction_folder(monomer_fasta, n=2):
    '''
    Create folders for predictions.
    '''
    try:
        os.mkdir(monomer_fasta[:-3]+'x'+str(n))
    except IOError:
        pass

    try:
        f = open(monomer_fasta[:-3]+'x'+str(n)+'/copy_unrelaxed_here', 'w')
        f.close()
    except IOError:
        pass

    try:
        os.mkdir(monomer_fasta[:-3]+'x'+str(n)+'/'+ \
                 monomer_fasta[:-3].split('/')[-1]+'x'+str(n))
    except IOError:
        pass

    try:
        f = open(monomer_fasta[:-3]+'x'+str(n)+'/'+ \
                 monomer_fasta[:-3].split('/')[-1]+'x'+str(n)+ \
                 '/copy_relaxed_here', 'w')
        f.close()
    except IOError:
        pass

    return


def subchains_prediction_folders(subchains, domains, multiplicities, \
                                 path_fasta):
    '''
    Create prediction folders for each given subchain and each of the given
    multiplicities.
    '''
    for sc in subchains:

        # generate each of the given multiplicity (e.g. x2, x3, ...)
        for m in multiplicities:
            if sc[0] == 1 and sc[1] == dom_n:
                monomer_fasta = path_fasta[:-4]+'.fa'
            elif sc[0] == sc[1]:
                monomer_fasta = path_fasta[:-4]+str(sc[0])+'.fa'
            else:    
                monomer_fasta = \
                        path_fasta[:-4]+str(sc[0])+str(sc[1])+'.fa'

            prediction_folder(monomer_fasta, m)

    return


if __name__ == "__main__":
    '''
    Create subchain oligomer fasta files.
    '''

    # import SymProFold libraries
    path_symprofold = lib.get_main_dir()
    sys.path.append(path_symprofold+'lib/')

    import bibfasta
    import ctl
    import filesystem

    multiplicities = [2, 3, 4, 6]

    root_path = str(pathlib.Path().absolute())+'/'
            # pre-assemblies folder


    create_predictionfolders = False

    if os.path.exists(root_path+'setting_create_predictionfolders.txt'):
        create_predictionfolders = True


    domain_fasta_files = sorted(glob.glob(root_path+'*_d.fa'))

    if len(domain_fasta_files) != 1:
        ctl.error('no unique domain fasta file found')

    domain_fasta = domain_fasta_files[0]
    filename = domain_fasta.split('/')[-1]

    sequ = bibfasta.get_seq_from_fasta(domain_fasta, False). \
                    replace("\r\n\r\n", "\r\n").replace("\n\n", "\n").strip()
    header = bibfasta.get_header(domain_fasta)

    if len(header) != 1:
        ctl.error('not uniquer header found')

    domains, zero = bibfasta.get_domains(domain_fasta, False)
    dom_n = len(domains)

    # generate standard set of subchains (standard subchains set)
    subchains, zero = generate_subchains(domains, sequ)


    # create fasta files for msa generation of standard subchains set
    # ---------------------------------------------------------------
    path_msa_gen = root_path+'msa_gen/'+filename
    filesystem.create_folder([path_msa_gen])


    # create fasta files to use as input for initial msa generation
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    write_subchains_to_fasta(subchains, domains, [2], path_msa_gen)


    # create fasta file containing all subchains to use as input for
    # hetero msa generation
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    write_subchains_to_single_fasta(subchains, domains, path_msa_gen)


    # create fasta files with all multiplicities and the corresponding folders
    # for the predictions
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    write_subchains_to_fasta(subchains, domains, multiplicities, domain_fasta)
    subchains_prediction_folders(subchains, domains, multiplicities, \
                                 domain_fasta)


    # create full ensemble of all possible subchains
    # ----------------------------------------------
    filesystem.create_folder([root_path+'subchain_ensemble/'])
    path_msa_gen = root_path+'subchain_ensemble/msa_gen/'+filename
    filesystem.create_folder([path_msa_gen])

    subchains = generate_subchain_ensemble(domains)

    path_fasta_ens = domain_fasta.split('/')
    path_fasta_ens.append(path_fasta_ens[-1])
    path_fasta_ens[-2] = 'subchain_ensemble'
    path_fasta_ens = '/'.join(path_fasta_ens)

    write_subchains_to_fasta(subchains, domains, [2], path_msa_gen)
    write_subchains_to_fasta(subchains, domains, multiplicities, \
                             path_fasta_ens)

    ctl.p('finished')
