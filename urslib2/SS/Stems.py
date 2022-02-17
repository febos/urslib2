

def Stem(model, stem_dict, ID, Type):

    wing1   = model.wings[Type][stem_dict['L'][ID]-1]
    wing2   = model.wings[Type][stem_dict['R'][ID]-1]

    pairseq = ','.join([model.bpairs[i-1]['TYPE'] for i in wing1['PAIRS']])

    jmol = wing1['JMOL'] + ',' + wing2['JMOL']

    return {'ID'        :             ID,
            'MODEL'     :              1,
            'CHAIN1'    : wing1['CHAIN'],
            'CHAIN2'    : wing2['CHAIN'],
            'LEFT'      :    wing1['ID'],
            'RIGHT'     :    wing2['ID'],
            'PAIRS'     : wing1['PAIRS'],
            'LSEQ'      :   wing1['SEQ'],
            'RSEQ'      :   wing2['SEQ'],
            'LEN'       :   wing1['LEN'],
            'PAIRSEQ'   :        pairseq,
            'TOWER'     :           None,
            'OLDSTEM'   :           None,
            'FULLSTEM'  :           None,
            'STEMS'     :              0,
            'OLDSTEMS'  :              0,
            'LOOPTYPE'  :             '',
            'LOOPPSEUDO':             '',
            'SCHEME'    :          '\\N',
            'SIGN'      :          '\\N',
            'MULT'      :          '\\N',
            'DIAGRAM'   :          '\\N',
            'NUMINDIAG' :          '\\N',
            'ECF'       :          '\\N',
            'NUMINECF'  :          '\\N',
            'JMOL'      :           jmol,
            'LOOPJMOL'  :             ''}

def Stems(model, Type):

    stems     = []
    stem_dict = {'L':{},'R':{}} # X: {stem_id: X_id,} 

    for wing in model.wings[Type]: stem_dict[wing['TYPE']][wing['STEM']] = wing['ID']

    for stem_id in stem_dict['L'].keys(): stems.append(Stem(model,stem_dict,stem_id,Type))

    stems.sort(key = lambda x: x['ID'])

    return stems

def relation(model, stems, Type):

    # Marking model.bpairs[i][FULLSTEM/OLDSTEM/REVSTEM/LUSTEM] (stems for bpairs)

    for stem in stems:

        for i in stem['PAIRS']:

            model.bpairs[i-1][Type+'STEM'] = stem['ID']

def relation2(model): # internal relations among different types of stems

    # stem-oldstem,stem-fullstem,oldstem-fullstem,oldstem-stems,fullstem-stems,fullstem-oldstems

    for stem in model.stems:

        pair_ind = stem['PAIRS'][0] - 1

        oldstem  = model.oldmask[pair_ind]
        fullstem = model.fullmask[pair_ind]

        stem['OLDSTEM'] = oldstem
        model.oldstems[oldstem-1]['STEMS'] += 1

        stem['FULLSTEM'] = fullstem
        model.fullstems[fullstem-1]['STEMS'] += 1

    for oldstem in model.oldstems:

        pair_ind = oldstem['PAIRS'][0] - 1

        fullstem = model.fullmask[pair_ind]

        oldstem['FULLSTEM'] = fullstem
        model.fullstems[fullstem-1]['OLDSTEMS'] += 1
    

def add(model):

    stems     = Stems(model,   'LU')
    oldstems  = Stems(model,  'OLD')
    fullstems = Stems(model, 'FULL')
    revstems  = Stems(model,  'REV')

    relation(model,     stems,     '')
    relation(model,  oldstems,  'OLD')
    relation(model, fullstems, 'FULL')
    relation(model,  revstems,  'REV')

    model.stems     = stems
    model.oldstems  = oldstems
    model.fullstems = fullstems
    model.revstems  = revstems

    relation2(model)
