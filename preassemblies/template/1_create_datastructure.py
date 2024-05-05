import glob
import lib
import os
import pathlib
import requests
import shutil
import sys
import time


def fetch(url, destination_paths):
    '''
    Fetch file from url.
    '''
    try:
        file = requests.get(url, stream=True)
    except requests.exceptions.SSLError:
        ctl.d('fetch: SSLError, waiting for restart')
        time.sleep(10)

        try:
            file=requests.get(url, stream=True)
        except requests.exceptions.SSLError:
            ctl.error('fetch: SSLError, second try')

    for filename in destination_paths:   
        with open(filename,'w') as f:
            f = open(filename, 'w')
            f.write(file.text)
            f.close()

    return


def fetch_fasta(gene_accesscode, path):
    '''
    Fetch FASTA files of given gene accession code from uniprot.org[1].

    [1] The UniProt Consortium, UniProt: the Universal Protein Knowledgebase in
    2023, Nucleic Acids Research, Volume 51, Issue D1, 6 January 2023,
    Pages D523–D531, https://doi.org/10.1093/nar/gkac1052
    '''
    url = 'https://rest.uniprot.org/uniprotkb/'+str(gene_accesscode)+'.fasta'
    ctl.d('fetch: '+url)
    fetch(url, [path+str(gene_accesscode)+'.fa'])

    return path+str(gene_accesscode)+'.fa'


def oligomer_fasta(monomer_fasta, n=2, file_output=None):
    '''
    Create homo oligomer fasta file.
    '''
    file = filesystem.get_file(monomer_fasta)

    if file_output == None:
        file_output = monomer_fasta

    f = open(file_output[:-3]+'_x'+str(n)+'.fa', 'w')

    for i in range(0, n):
        for l in file:
            f.write(l)

    f.close()

    return


def prediction_folder(monomer_fasta, n=2):
    '''
    Create folders for predictions.
    '''
    try:
        os.mkdir(monomer_fasta[:-3]+'_x'+str(n))
    except IOError:
        pass

    try:
        f = open(monomer_fasta[:-3]+'_x'+str(n)+'/copy_unrelaxed_here', 'w')
        f.close()
    except IOError:
        pass

    try:
        os.mkdir(monomer_fasta[:-3]+'_x'+str(n)+'/'+ \
                 monomer_fasta[:-3].split('/')[-1]+'_x'+str(n))
    except IOError:
        pass

    try:
        f = open(monomer_fasta[:-3]+'_x'+str(n)+'/'+ \
                 monomer_fasta[:-3].split('/')[-1]+'_x'+str(n)+ \
                 '/copy_relaxed_here', 'w')
        f.close()
    except IOError:
        pass

    return


def oligomer_fasta_set(monomer_fasta, create_predictionfolders, \
                       file_output=None):
    '''
    Create set of homo oligomer fasta files with different oligomeric states.
    '''
    seq = bibfasta.get_seq_from_fasta(monomer_fasta)

    if create_predictionfolders:
        prediction_folder(monomer_fasta, 1)

    oligomer_fasta(monomer_fasta, 2, file_output)

    if create_predictionfolders:
        prediction_folder(monomer_fasta, 2)

    if len(seq) <= 1500:
        oligomer_fasta(monomer_fasta, 3, file_output)

        if create_predictionfolders:
            prediction_folder(monomer_fasta, 3)

    if len(seq) <= 1200:
        oligomer_fasta(monomer_fasta, 4, file_output)

        if create_predictionfolders:
            prediction_folder(monomer_fasta, 4)

    if len(seq) <= 800:
        oligomer_fasta(monomer_fasta, 6, file_output)

        if create_predictionfolders:
            prediction_folder(monomer_fasta, 6)

    if len(seq) <= 600:
        oligomer_fasta(monomer_fasta, 8, file_output)

        if create_predictionfolders:
            prediction_folder(monomer_fasta, 8)

    return


if __name__ == "__main__":

    # import SymProFold libraries
    path_symprofold = lib.get_main_dir()
    sys.path.append(path_symprofold+'lib/')

    import ctl
    import bibfasta
    import filesystem


    root_path = str(pathlib.Path().absolute())+'/'
            # pre-assemblies folder


    create_predictionfolders = False

    if os.path.exists(root_path+'setting_create_predictionfolders.txt'):
        create_predictionfolders = True


    folder_content = sorted(glob.glob(root_path+'*'))
            # content of pre-assemblies folder

    files = [] # files in pre-assemblies folder

    for f in folder_content:
        if os.path.isfile(f):
            files.append(f)


    # create list of gene accession code files
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    gene_accesscode_files = [] # list of gene accession code files

    for f in files:
        filename = f.split('/')[-1]
        gene_accesscode = [filename]
        
        if ('.' not in filename) and len(filename) <= 15:
            gene_accesscode_files.append(f)


    # create FL homo oligomer fasta files
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if len(gene_accesscode_files) == 0:
        fasta_files = sorted(glob.glob(root_path+'*.fa*'))
        fasta_files_ = []

        for f in fasta_files:
            if os.path.isfile(f) and ('_x' not in f.split('/')[-1]):
                fasta_files_.append(f)

        fasta_files = fasta_files_

        for f in fasta_files:
            filename = f.split('/')[-1].split('_')[0]

            if len(filename) <= 15:

                # creation of FL homo oligomer fasta files
                oligomer_fasta_set(f, create_predictionfolders, \
                                   file_output=None)


    elif len(gene_accesscode_files) == 1:
        for f in gene_accesscode_files:
            filename = f.split('/')[-1]
            gene_accesscode = [filename]

            if ('.' not in filename) and len(filename) <= 15:
                fasta_filename = fetch_fasta(gene_accesscode[0], root_path)

                # creation of FL homo oligomer fasta files
                oligomer_fasta_set(fasta_filename, create_predictionfolders, \
                                   file_output=None)


    elif len(gene_accesscode_files) != 1:
        ctl.e(gene_accesscode_files)
        ctl.error('not exactly 1 gene accession code file in folder')


    ctl.p('finished')
