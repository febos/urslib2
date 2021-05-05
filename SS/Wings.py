
try:
    import urslib2.SS.Mask       as Mask
    import urslib2.SS.ChainOrder as ChainOrder
except:
    try:
        import SS.Mask       as Mask
        import SS.ChainOrder as ChainOrder
    except:
        import Mask
        import ChainOrder

def MakeJmol(dssrs):

    if all([x.endswith('.') for x in dssrs]) and not any(['-' in x for x in dssrs]):
        if len(dssrs)>1: return dssrs[0].split('.')[2]+'-'+dssrs[-1].split('.')[2]+':'+dssrs[0].split('.')[0]
        else           : return dssrs[0].split('.')[2]+':'+dssrs[0].split('.')[0]

    else:
        res = []

        for x in dssrs:
            jmolx = x.split('.')
            if not jmolx[-1]: res.append(jmolx[2]+':'+jmolx[0])
            else:             res.append(jmolx[2]+'^'+jmolx[3]+':'+jmolx[0])
        return ','.join(res)

def Wing(model, stem_id, pairs, side):

    Where = {'RES':'SEQ2','LIGANDS':'LIGSEQ'}

    if side == 'L':

        start = model.bpairs[pairs[0]-1]['NUCL1']
        end   = model.bpairs[pairs[-1]-1]['NUCL1']
        chain = start[0].split('.')[0]

    else: # == 'R'

        start = model.bpairs[pairs[-1]-1]['NUCL2']
        end   = model.bpairs[pairs[0]-1]['NUCL2']
        chain = start[0].split('.')[0]

    if not start[2] < end[2]: start,end = end,start

    seq = ','.join(model.chains[chain][Where[start[1]]][start[2]:end[2]+1])

    Len = len(pairs)

    jmol  = MakeJmol([x['DSSR'] for x in model.chains[chain][start[1]][start[2]:end[2]+1]])

    return {'ID'       :       0,
            'CHAIN'    :   chain,
            'STEM'     : stem_id,
            'ANOTHER'  :    None,
            'PREV'     :    None,
            'NEXT'     :    None,
            'PREVW'    :    None,
            'NEXTW'    :    None,
            'PAIRS'    :   pairs,
            'TYPE'     :    side,
            'START'    :   start,
            'END'      :     end,
            'SEQ'      :     seq,
            'LEN'      :     Len,
            'ECF'      :   '\\N',
            'NUMINECF' :   '\\N',
            'JMOL'     :    jmol}

def Wings(model, mask):

    wings = []
    stems = {} # {stem_id : [pair1_id,pair2_id,...]}
    Max   = max(mask+[0])
    
    if not Max: return wings

    for i in range(Max): stems[i+1] = []

    for i in range(len(mask)):

        if mask[i]:

            stems[mask[i]].append(i+1) # adding id of pair into current stem

    for stem_id in stems:

        wings.append(Wing(model,stem_id,stems[stem_id],'L'))
        wings.append(Wing(model,stem_id,stems[stem_id],'R'))

    return wings

def sort_and_relation(model, wings, Type):

    # Sorting

    where_dict = {'RES':1,'LIGANDS':2}

    wings.sort( key = lambda x: [where_dict[x['START'][1]],model.chain_order[x['CHAIN']], x['START'][2]])

    # Rebuilding wing['TYPE'] (if chain_order has mistakes)

    seen_stems = []

    for wing in wings:

        if wing['STEM'] not in seen_stems:

            seen_stems.append(wing['STEM'])
            wing['TYPE'] = 'L'

        else: wing['TYPE'] = 'R'

    # Marking ID, ANOTHER, PREVW and NEXTW

    for i in range(len(wings)): wings[i]['ID'] = i+1

    anothers = {'L':{},'R':{}}
    mirror   = {'L':'R','R':'L'}

    for w in wings: anothers[w['TYPE']][w['STEM']] = w['ID']

    for i in range(len(wings)):

        wings[i]['ANOTHER']  = anothers[mirror[wings[i]['TYPE']]][wings[i]['STEM']]

        if   i == 0:

            if wings[i]['CHAIN'] == wings[i+1]['CHAIN']: wings[i]['NEXTW'] = i+2

        elif i == len(wings) - 1:

            if wings[i]['CHAIN'] == wings[i-1]['CHAIN']: wings[i]['PREVW'] = i

        else:

            if wings[i]['CHAIN'] == wings[i+1]['CHAIN']: wings[i]['NEXTW'] = i+2
            if wings[i]['CHAIN'] == wings[i-1]['CHAIN']: wings[i]['PREVW'] = i

    # Marking model.chains[ch][RES/LIGANDS][i][WING/FULLWING/OLDWING] (wings for nucleotides)

    if Type not in ('FULLWING','REVWING'):

        for w in wings:

            for nucl in model.chains[w['CHAIN']][w['START'][1]][w['START'][2]:w['END'][2]+1]:

                nucl[Type] = w['ID']    # nucl['WING'/'OLDWING']

    if Type == 'FULLWING':

        for w in wings:

            for nucl in model.chains[w['CHAIN']][w['START'][1]][w['START'][2]:w['END'][2]+1]:

                nucl['FSTEMS'] += 1

def add(model):

    wings     = Wings(model, model.mask)
    fullwings = Wings(model, model.fullmask)
    oldwings  = Wings(model, model.oldmask)
    revwings  = Wings(model, model.revmask)

    sort_and_relation(model,     wings,     'WING')
    sort_and_relation(model, fullwings, 'FULLWING')
    sort_and_relation(model,  oldwings,  'OLDWING')
    sort_and_relation(model,  revwings,  'REVWING')

    model.wings = {}

    model.wings['LU']     = wings
    model.wings['FULL']   = fullwings
    model.wings['OLD']    = oldwings
    model.wings['REV']    = revwings

