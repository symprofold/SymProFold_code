import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import pathlib
import pickle
import shutil


def get_file(file, attempt=0):
    '''
    Open file and create list with lines.
    '''
    attempt += 1
    ra = 5
 
    try:
        f = open(file, 'r')
        ra = f.readable()
        if ra == True:
         lin = f.readlines()
        else:
            ra = 5
        f.close()

    except IOError:
        print('ERROR: file: '+file+', attempt:'+str(attempt))
        print('retry')

        time.sleep(0.1*attempt)
        return get_file(file)

    if ra != True:
        ctl.d(ra) 

    return lin


def sort_models(ptm_file):
    '''
    Sort models by ipTM+pTM score from ptm_file.
    '''
    models = []
    modelssort = []
 
    started = False

    for l in ptm_file:
        l = l.strip()
         
        if l[0:7] == '"model_':
            started = True

        if l[0:7] != '"model_' and started == True:
            break

        if started == True:
            l = l.split('": ')
            models.append([l[0][1:], float(l[1][:-1])])
         
    alreadysorted = []

    for i,m in enumerate(models):
        cand = ["", 0]

        for j,n in enumerate(models):
            if not n[0] in alreadysorted:
                if n[1] > cand[1]:
                    cand = n

        modelssort.append(cand)
        alreadysorted.append(cand[0])

    return modelssort


def gen_diagram(mode, ds, de, prefix, path_part1):
    '''
    Generate and save diagram for models.
    '''
    work_dir = str(pathlib.Path().absolute())+'/'
    pkl_files = sorted(glob.glob(work_dir+'result_*.pkl'))

    if len(pkl_files) == 0:
        return

    maxr = '/'+str(de)+': '

    i = -1
    matrix = []

    for j in range(ds, de):
        if j%5 == 0:
            i += 1
            matrix.append([])
           
        matrix[i].append(['1', str(j)+maxr])

    if i == 0:
        figu, axes = plt.subplots(nrows=len(matrix), ncols=len(matrix[0]), \
                                  figsize=(15, 2.5*(i+1)))
    else:
        figu, axes = plt.subplots(nrows=len(matrix), ncols=len(matrix[0]), \
                                  figsize=(15, 2.0*(i+1)))

    plt.tight_layout()

    i = -1

    for j,m in enumerate(modelssort[ds:de]):
        if j%5 == 0:
           i += 1

        pkl_files = sorted(glob.glob(work_dir+'result_'+m[0]+'.pkl'))

        if len(pkl_files) == 0:
            continue

        pkl_file = pkl_files[0]

        with open(path_part1+pkl_file, 'rb') as f:
            d1 = pickle.load(f)

        if mode == 'stat':
            data = d1['predicted_aligned_error']

        if mode == 'dist':
            distance_bins = np.append(0, d1['distogram']['bin_edges'])
            distance_logits = d1['distogram']['logits']
            distance_matrix = distance_bins[distance_logits.argmax(-1)]
            data = distance_matrix

        ptm = float(d1['ptm'])
        iptm = float(d1['iptm'])

        mpae = d1['max_predicted_aligned_error']
        if mpae != 31.75:
            print('ERROR: gen_diagram: mpae')
            sys.exit()

        if mode == 'stat':
            vmin = 0
            vmax = 32
            title1 = matrix[i][j%5][1]+'pTM:'+str(round(ptm,2))+ \
                     ' ipTM:'+str(round(iptm, 2))+' comb.:'+str(round(m[1], 2))
        if mode == 'dist':     
            vmin = 0
            vmax = 20
            title1 = matrix[i][j%5][1]+'contact map'

        if len(matrix) == 1:
            im = axes[j].imshow(data, vmin=vmin, vmax=vmax)
            axes[j].tick_params(axis='x', direction='in')
            axes[j].tick_params(axis='y', direction='in')
            axes[j].set_title(title1, fontsize=9, fontweight='bold', pad=6)
            if i == 0 and j == 0:
                axes[j].figure.colorbar(im)
        else:
            im = axes[i, j%5].imshow(data, vmin=vmin, vmax=vmax)
            axes[i, j%5].tick_params(axis='x', direction='in')
            axes[i, j%5].tick_params(axis='y', direction='in')
            axes[i, j%5].set_title(title1, fontsize=9, fontweight='bold', \
                                   pad=6)
            if i == 0 and j == 0:
                axes[i, j%5].figure.colorbar(im)

    figu.savefig(path_part1+prefix+'_'+('00'+str(ds))[-2:]+'-'+ \
                 ('00'+str(de))[-2:]+'.png')

    return


if __name__ == "__main__":
    work_dir = os.path.dirname(os.path.realpath(__file__))+'/'
    path_part1 = './'
    ptm_file = 'ranking_debug.json'

    pkl_files = sorted(glob.glob(work_dir+'*.pkl'))
    model_files_ranked = sorted(glob.glob(work_dir+'*.pkl'))


    if os.path.isfile(work_dir+ptm_file) == True and len(pkl_files) <= 1:

        import json

        f = open(work_dir+ptm_file, 'r')
        data = json.load(f)
        f.close()

        model_files = sorted(glob.glob(work_dir+'ranked_*.pdb'))

        for i,m in enumerate(model_files):
            iptmptm = float(data['iptm+ptm'][data['order'][i]])
            model_i = int(\
                        m.split('/')[-1].split('ranked_')[1].split('.pdb')[0])

            if i != model_i:
                continue

            model_file = get_file(m)
            h_atoms = 0
            atoms = 0

            for l in model_file:
                atoms += 1

                if len(l) > 13:
                    if l[13] == 'H':
                        h_atoms += 1

            prefix = '' # relaxed

            if h_atoms/atoms < 0.1:
                prefix = 'unrelaxed_' # unrelaxed

            filename_new = 'rank'+('00'+str(i))[-2:]+'_'+ \
                                (str(round(iptmptm,3))+'000')[0:5]+'.pdb'
            shutil.copyfile(m, work_dir+prefix+filename_new)

    elif os.path.isfile(work_dir+ptm_file) == True:
        ptm_file_lines = get_file(path_part1+ptm_file)
        modelssort = sort_models(ptm_file_lines)

        mode = 0 # 0: relaxed; 1: unrelaxed
        model_files = sorted(glob.glob(path_part1+'relaxed*.pdb'))

        if len(model_files) == 0:
             mode = 1

        for i,m in enumerate(modelssort):
            if mode == 0:
                shutil.copyfile(path_part1+'relaxed_'+m[0]+'.pdb', \
                                path_part1+'rank'+('00'+str(i))[-2:]+'_'+ \
                                (str(round(m[1],3))+'000')[0:5]+'.pdb')
            if mode == 1:
                shutil.copyfile(path_part1+'unrelaxed_'+m[0]+'.pdb', \
                                path_part1+'unrelaxed_rank'+ \
                                ('00'+str(i))[-2:]+'_'+ \
                                (str(round(m[1],3))+'000')[0:5]+'.pdb')

        gen_diagram('stat', 0, 5, '', path_part1)
        gen_diagram('stat', 0, 25, '', path_part1)
        gen_diagram('dist', 0, 5, 'contactmap', path_part1)
        gen_diagram('dist', 0, 25, 'contactmap', path_part1)

    else:
        pkl_files = sorted(glob.glob(path_part1+'*.pkl'))

        for i,p in enumerate(pkl_files):
            if p == path_part1+'features.pkl':
                  continue

            with open(p, 'rb') as f:
                d1 = pickle.load(f)

            print(p+':')
            print('ranking_confidence: '+str(d1['ranking_confidence']))

            postfix = p.split('result')[1][:-4]+'.pdb'

            if os.path.isfile(path_part1+'relaxed'+postfix) == True:
                shutil.copyfile(path_part1+'relaxed'+postfix, \
                    path_part1+'rank'+('00'+str(50+i))[-2:]+'_'+ \
                    (str(round(d1['ranking_confidence'],3))+'000')[0:5] \
                    +'.pdb')

            elif os.path.isfile(path_part1+'unrelaxed'+postfix) == True:
                shutil.copyfile(path_part1+'unrelaxed'+postfix, \
                    path_part1+'unrelaxed_rank'+('00'+str(50+i))[-2:]+ \
                    '_'+ \
                    (str(round(d1['ranking_confidence'],3))+'000')[0:5]+ \
                    '.pdb')

    print('finished')  
