import ctl
import model_preparation.sstruct


'''
Module providing functions for model preparation to determine beta sheets.
'''

def intermolecular_betasheet_fraction(model_id, path_export, sess, verbous):
    '''
    Determine fraction of beta sheet residues that is involved in
    intermolecular sheets.
    '''

    # determine secondary structure data using dssp in ChimeraX
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # residues involved in beta sheet strands
    dssp = model_preparation.sstruct.get_dssp(model_id, path_export, sess)

    # residues involved in medium-length beta sheet strands
    dssp_minlen = model_preparation.sstruct.get_dssp(model_id, path_export, \
                        sess, ' minStrandLen 5')


    # restrict to intermolecular beta sheets
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # residues involved in intermolecular beta sheet strands
    sheet_intermol_res, sheet_intramol_res = \
                        intermolecular_betasheet_residues( \
                            model_id, path_export, dssp, sess)

    # residues involved in intermolecular medium-length beta sheet strands
    sheet_intermol_minlen_res, zero = \
                        intermolecular_betasheet_residues( \
                            model_id, path_export, dssp_minlen, sess)

    sheet_res = list(set(sheet_intramol_res+sheet_intermol_res))

    sheet_res_n = len(sheet_res)
    sheet_intramol_res_n = len(sheet_intramol_res)
    sheet_intermol_res_n = len(sheet_intermol_res)
    sheet_intermol_minlen_res_n = len(sheet_intermol_minlen_res)

    if sheet_res_n > 0:
        intermol_res_fract = sheet_intermol_res_n/sheet_res_n
        intermol_minlen_res_fract = sheet_intermol_minlen_res_n/sheet_res_n
    else:
        intermol_res_fract = 0
        intermol_minlen_res_fract = 0

    if verbous:
        ctl.p('sheet_res')
        ctl.p(sheet_res)
        ctl.p(sheet_res_n)
        ctl.p('sheet_intramol_res')
        ctl.p(sheet_intramol_res)
        ctl.p(sheet_intramol_res_n)
        ctl.p('sheet_intermol_res')
        ctl.p(sheet_intermol_res)
        ctl.p(sheet_intermol_res_n)
        ctl.p('inermol_res_fract')
        ctl.p(intermol_res_fract)

    sheetintermol_infix = '_sheetintermol'+ \
                                str(round(intermol_minlen_res_fract, 3))+'_'+ \
                                str(round(intermol_res_fract, 3))+'_'+ \
                                str(sheet_res_n)

    return sheetintermol_infix, intermol_res_fract, sheet_intermol_res_n


def intermolecular_betasheet_residues(model_id, path_export, dssp, sess):
    '''
    Determine residues involved in intermolecular beta sheets.
    '''
    sheet_intramol_res = []
    sheet_intermol_res = []

    for d in dssp:

        # antiparallel beta sheets
        # ~~~~~~~~~~~~~~~~~~~~~~~~
        if ' antiparallel ' in d:
            ranges = d.split(' antiparallel ')
            range0 = ranges[0].split(' -> ')
            begin0 = int(range0[0].split(':')[1].strip())
            end0 = int(range0[1].split(':')[1].strip())
            chain0 = range0[0].split('/')[1].split(':')[0].strip()

            range1 = ranges[1].split(' -> ')
            begin1 = int(range1[0].split(':')[1].strip())
            end1 = int(range1[1].split(':')[1].strip())
            chain1 = range1[0].split('/')[1].split(':')[0].strip()

            r0 = [chain0+str(r) for r in range(begin0, end0+1)]
            r1 = [chain1+str(r) for r in range(begin1, end1+1)]

            if chain0 != chain1:
                for r in r0:
                    if r not in sheet_intermol_res:
                        sheet_intermol_res.append(r)

                for r in r1:
                    if r not in sheet_intermol_res:
                        sheet_intermol_res.append(r)

            elif chain0 == chain1:
                for r in r0:
                    if r not in sheet_intramol_res:
                        sheet_intramol_res.append(r)

                for r in r1:
                    if r not in sheet_intramol_res:
                        sheet_intramol_res.append(r)

            else:
                ctl.error('intermolecular_betasheet_fraction: '+ \
                            'error in chain ids')


        # parallel beta sheets
        # ~~~~~~~~~~~~~~~~~~~~
        elif ' parallel ' in d:
            ranges = d.split(' parallel ')
            range0 = ranges[0].split(' -> ')
            begin0 = int(range0[0].split(':')[1].strip())
            end0 = int(range0[1].split(':')[1].strip())
            chain0 = range0[0].split('/')[1].split(':')[0].strip()

            range1 = ranges[1].split(' -> ')
            begin1 = int(range1[0].split(':')[1].strip())
            end1 = int(range1[1].split(':')[1].strip())
            chain1 = range1[0].split('/')[1].split(':')[0].strip()

            r0 = [chain0+str(r) for r in range(begin0, end0+1)]
            r1 = [chain1+str(r) for r in range(begin1, end1+1)]

            if chain0 != chain1:
                for r in r0:
                    if r not in sheet_intermol_res:
                        sheet_intermol_res.append(r)

                for r in r1:
                    if r not in sheet_intermol_res:
                        sheet_intermol_res.append(r)

            elif chain0 == chain1:
                for r in r0:
                    if r not in sheet_intramol_res:
                        sheet_intramol_res.append(r)

                for r in r1:
                    if r not in sheet_intramol_res:
                        sheet_intramol_res.append(r)

            else:
                ctl.error('intermolecular_betasheet_fraction: '+ \
                            'error in chain ids')

    return sheet_intermol_res, sheet_intramol_res
