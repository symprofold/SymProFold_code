import ctl
import filesystem

import os
import re


def get_cif_col(l, colnumb):
    """ Get column of given column number in line of cif file. """

    col = 0

    # replace multiple spaces by a single space
    l = re.sub(r'\s{2,}', ' ', l)

    l = l.split(' ')
    col = int(l[colnumb])

    return col


def sort_lines(out):
    """ Sort lines of cif file regarding residue id. """

    out2 = [[] for i in range(0,100000)]
        # 100000: max residue id
    out3 = []
    atom_nr = 0

    res_id_max = 0
    for l in out:
        res_id = int(get_cif_col(l, 13))
        if res_id > res_id_max:  
            res_id_max = res_id


    for l in out:
        res_id = int(get_cif_col(l, 13))

        res_id_str = str(res_id)+'     '
        pos_sep1 = l.find(' . ')
        l2 = l[:pos_sep1+11]+res_id_str[:len(str(res_id_max))+1]+ \
             l[pos_sep1+10+len(str(res_id_max))+1:]        
                    # len(str(res_id_max))+1:
                    # adjust column width to the number of digits

        out2[res_id].append(l2)

    atom_nr_max = 0
    for o in out2:
        for l in o:
            atom_nr_max += 1

    for o in out2:
        for l in o:
            atom_nr += 1
            l2 = l[:5]+(str(atom_nr)+'      ')[:len(str(atom_nr_max))+1]+ \
                 l[5+len(str(atom_nr_max)):]
                    # len(str(atom_nr))+1:
                    # adjust column width to the number of digits

            out3.append(l2)

    return out3


def postprocess(import_file, export_file):
    """ Postprocess of cif file. """
    
    f = filesystem.get_file(import_file)
    current_chain_id = 0

    out = []
    out2 = []
    collect_active = 0
    
    for i,l in enumerate(f):
        if l[:4] == 'HELX':
            continue
        if l[:1] == '?':
            continue        
        
        if l[:4] == 'ATOM':
            out2.append(l)
            collect_active = 1

        else:
            if collect_active == 1:
                collect_active = 0
                out2 = sort_lines(out2)
                for l2 in out2:
                    out.append(l2)
                out2 = []
            else:
                if l[:10] == 'data_chain':
                    current_chain_id = int(l[10:])
                    ctl.d(current_chain_id)

            out.append(l)            

    fo = open(export_file, 'w')
    for l in out:
        fo.write(l)
    fo.close()

    return


def compatibility_cif_export(export_path, export_model_id, session, cif_postprecess=True):
    """ Compatibility cif export. """

    export_path = filesystem.clean_path(export_path)
    debug_mode = False

    filesystem.create_folder([export_path[:-4]+'_tmp.cif'])
    session.run('save "'+export_path[:-4]+'_tmp.cif" #'+ \
                str(export_model_id)+' ')
    session.init()
    session.run('open "'+export_path[:-4]+'_tmp.cif"')
    session.run('split #1')
    for i in range(1,200):
        session.run('rename #1.'+str(i)+' chain'+str(i))

    session.run('save "'+export_path[:-4]+'_tmp.cif'+'" #'+str(1)+' ')

    if debug_mode:
        session.run('save "'+export_path[:-4]+'_tmp0.cif'+'" #'+str(1)+' ')

    if cif_postprecess:
        if debug_mode:
            postprocess(export_path[:-4]+'_tmp.cif', export_path[:-4]+'_tmp2.cif')
        postprocess(export_path[:-4]+'_tmp.cif', export_path[:-4]+'_tmp.cif')

    session.init()
    session.run('open "'+export_path[:-4]+'_tmp.cif"')
    session.run('save "'+export_path+'" #'+str(1)+' ')

    if not debug_mode:
        os.unlink(export_path[:-4]+'_tmp.cif')

    return


def compatibility_cif_export_combine(export_path, export_model_id, session, cif_postprecess=True):
    """ Compatibility cif export with combination of models. """

    export_path = filesystem.clean_path(export_path)
    debug_mode = False

    filesystem.create_folder([export_path[:-4]+'_tmp.cif'])
    session.run('save "'+export_path[:-4]+'_tmp.cif" #'+str(export_model_id))
    session.init()
    session.run('open "'+export_path[:-4]+'_tmp.cif"')
    session.run('split #1')
    for i in range(1,200):
        session.run('rename #1.'+str(i)+' chain'+str(i))

    session.run('save "'+export_path[:-4]+'_tmp.cif'+'" #'+str(1)+' ')

    if cif_postprecess:
        if debug_mode:
            postprocess(export_path[:-4]+'_tmp.cif', export_path[:-4]+'_tmp2.cif')
        postprocess(export_path[:-4]+'_tmp.cif', export_path[:-4]+'_tmp.cif')

    session.init()
    session.run('open "'+export_path[:-4]+'_tmp.cif"')

    combination_model_id = session.last_id()+1
    session.run('combine #1-'+str(session.last_id())+ \
            ' close false modelId #'+str(combination_model_id))

    session.run('save "'+export_path+'" #'+str(combination_model_id))

    if not debug_mode:
        os.unlink(export_path[:-4]+'_tmp.cif')

    return


def cif_add_symm(source_path, export_path, a, b, c,
                 alpha, beta, gamma, symmgroup):
    """ Add symmetry information to cif file. """

    a = str(round(a, 3))
    b = str(round(b, 3))
    c = str(round(c, 3))

    alpha = str(round(alpha, 2))
    beta = str(round(beta, 2))
    gamma = str(round(gamma, 2))

    ctl.d(a)
    ctl.d(b)
    ctl.d(c)    
    ctl.d(alpha)    

    sym_lines = ["#", \
    "_cell.entry_id 01", \
    "_cell.length_a "+a, \
    "_cell.length_b "+b, \
    "_cell.length_c "+c, \
    "_cell.angle_alpha "+alpha, \
    "_cell.angle_beta  "+beta, \
    "_cell.angle_gamma "+gamma, \
    "_symmetry.entry_id 01", \
    "_symmetry.space_group_name_H-M '"+symmgroup+"'"]

    ctl.d(sym_lines)
    f = filesystem.get_file(source_path)
    ctl.p('step 2')
    fo = open(export_path, 'w')
    for l in f:
        fo.write(l)
        if l[0:len('_software.pdbx_ordinal')] == '_software.pdbx_ordinal':
            for s in sym_lines:
                fo.write(s+"\r\n")
    ctl.p('step 3')
    fo.close()
    ctl.p('step 4')

    return
