import ctl
import filesystem


def get_seq_from_fasta(file, remove_linebreaks=True):
    ''' Get sequence from fasta file. '''

    seq = ''
    f = filesystem.get_file(file)

    for i,l in enumerate(f):
        if i == 0:
            if l[0] != '>':
                ctl.error('get_seq_from_fasta: first line is not header')
            
        if i == 0:
            continue

        if remove_linebreaks == True:
            seq += l.strip()
        else:    
            seq += l

    return seq


def get_domains(file, remove_linebreaks=True):
    '''
    Get domains from fasta file in which domains are separated by empty lines.
    '''
    domain_seq = []
    domain_boundaries = [[1,0]]
    f = filesystem.get_file(file)

    last_resid = 0

    for i,l in enumerate(f):
        if i == 0:
            if l[0] != '>':
                ctl.e('get_seq_from_fasta: first line is not header')
            domain_seq.append('')

        elif l.strip() == '':
            domain_seq.append('')
            domain_boundaries.append([domain_boundaries[-1][1]+1, \
                                      domain_boundaries[-1][1]])

        else:
            if remove_linebreaks == True:
                domain_seq[-1] += l.strip()
            else:
                domain_seq[-1] += l

            domain_boundaries[-1][1] += len(l.strip())

    if domain_seq[-1] == '':
        del domain_seq[-1]

    if domain_boundaries[-1][0]-1 == domain_boundaries[-1][1]:
        del domain_boundaries[-1]


    return domain_seq, domain_boundaries


def get_header(file):
    ''' Get all sequence header from fasta file. '''

    header = []
    f = filesystem.get_file(file)

    for i,l in enumerate(f):
        if l[0] == '>':
            header.append(l.strip())

    return header
