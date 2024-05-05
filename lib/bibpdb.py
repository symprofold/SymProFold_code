import ctl
import filesystem

import glob
import math


'''
Module providing functions for handling pdb files.
'''

def get_f(l, field):
    '''
    Get field in line of pdb file.

    fields: 'name', 'res', 'chainid', 'resnr', 'x', 'atnr', 'element', 'bfact'
    '''
    if field == 'name':
        return l[13:16].strip()
    if field == 'res':
        return l[17:20].strip()
    if field == 'chainid':
        return l[21:22]
    if field == 'resnr':
        return int(l[22:30])
    if field == 'x':
        return float(l[30:38])
    if field == 'y':
        return float(l[38:46])
    if field == 'z':
        return float(l[46:54])
    if field == 'atnr':
        return int(l[4:11])
    if field == 'element':
        return l[77:78]
    if field == 'bfact':
        return float(l[60:66])

    return false


def set_f(l, field, nr):
    '''
    Set field in line of pdb file.

    fields: 'chainid', 'resnr', 'bfact'
    '''
    r = field

    if field == 'chainid':
        if len(nr) != 1:
            ctl.error('ERROR: chainid length')

        r = l[0:21]+str(nr)+l[22:]
    if field == 'resnr':
        r = l[0:22]+('          '+str(nr))[-4:]+l[26:]
    if field == 'bfact':
        if nr > 999.99:
            nr = 999.99
        r = l[0:60]+('          '+str(round(nr,2)))[-6:]+l[66:]

    return r


def lddt_to_bfact(file_in, file_out, stop_at_TER=True): 
    ''' Convert LDDT to b factors in given pdb file. '''

    c1 = open_pdb(file_in, stop_at_TER)
    w = []
    bfact = {}

    for i,l in enumerate(c1):
        lddt = float(l[60:70])

        lRMSD = 1.5*math.exp(4*(0.7-lddt/100))*1 # regular formula
        bfact[int(l[22:30])] = round(lRMSD*lRMSD*8/3*3.1415*3.1415, 2)

    return bfact


def get_rmsds(file):
    ''' Calculate RMSDs of residues in given pdb file. '''

    rmsd = {}

    bfact = lddt_to_bfact(file, file+'rmsd.pdb', False)
    for b in bfact:
        rmsd[b] = round(math.sqrt(float(bfact[b])/(8/3*3.1415*3.1415)), 2)

    return rmsd


def trim_pdb(c, stop_at_TER):
    '''
    Trim pdb file.
    E.g. Leave only 'ATOM' entries, remove 'TER'
    '''
    r = []
     
    for i, l in enumerate(c):
        if len(l) > 10 and l[0:4] == 'ATOM':
            r.append(l)      
        if l[0:3] == 'TER' and stop_at_TER == True:
            break 

    return r


def clean_pdb(file, fileout=''):
    '''
    Clean/reformat pdb file.
    Reduce file to the entries 'ATOM', 'TER' and 'END'.
    '''
    coord0 = filesystem.get_file(file)
    if fileout == '':
        fileout = file

    # replace 'ENDMDL' code with 'TER' to ensure the presence of a
    # termination code
    coord0 = [i.replace('ENDMDL', 'TER') for i in coord0]


    # remove duplicates of 'TER' in 2 consecutive lines
    coord = []
    prev_line = ''

    for l in coord0:
        if len(l) >= 3:
            if l[0:3] == 'TER' == (prev_line+'   ')[0:3]:
                prev_line = l
                continue

        coord.append(l)
        prev_line = l


    # write lines starting with 'ATOM', 'TER' or 'END' to files
    f = open(fileout, 'w')
    for l in coord:
        if (l[0:4] == 'ATOM' or l[0:3] == 'TER' or l[0:3] == 'END'):
            f.write(l) 

    f.close()

    return


def open_pdb(fileIn, stop_at_TER=True):
    ''' Open pdb file. '''

    try:
        f = open(fileIn, 'r')
        ra = f.readable()
        if ra == True:
            c1 = f.readlines()
        f.close()
    except IOError:
        return []
         
    c1 = trim_pdb(c1, stop_at_TER)
    c2 = c1

    return c2


def reorder_chainids(file, order, fileout=''):
    '''
    Reorder chain ids in pdb file according to order given in list "order".

    list "order": Key defines the new chainid (alphabet letter converted
                  to numerical value), the value stands for the old chainid
                  (alphabetic letter converted to numerical value).
    '''
    coord = filesystem.get_file(file)
    if fileout == '':
        fileout = file

    chain_id_names = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    coord2 = []

    for l in coord:
        if l[0:4] != 'ATOM':
            coord2.append(l)
            continue

        chainid = get_f(l, 'chainid')

        # get position in alphabet
        al_pos = chain_id_names.find(chainid)

        for i,o in enumerate(order):
            if o == al_pos:
                newchainid = chain_id_names[i]
                break
         
        l2 = set_f(l, 'chainid', newchainid)
        coord2.append(l2)
          
    f = open(fileout, 'w')
    for l in coord2:
        if l[0:4] == 'ATOM' or l[0:3] == 'TER' or l[0:3] == 'END':
            f.write(l) 

    f.close()

    return


def get_multimer_n(file):
    '''
    Get number of chains in a given pdb file.
    Use case: e.g. get multiplicity of a multimer in a pdb file.
    '''
    n = 0
    f = filesystem.get_file(file)

    # remove all lines except 'ATOM', 'END' and 'TER'
    code_cleaned = ['ATOM', 'END', 'TER']
    f = [l for l in f if any(l.startswith(code) for code in code_cleaned)]

    for i, l in enumerate(f):
        if i >= 1:   
            if f[i][0:3] == 'TER' or \
                    (f[i][0:3] == 'END' and f[i-1][0:3] != 'END' and \
                    f[i-1][0:3] != 'TER'):
                n += 1

    return n
