
try:    import urslib2.SS.Graph as Graph
except:
    try:    import SS.Graph as Graph
    except: import Graph

def StemNPair(fstem1, fstem2, nucl):

    return {'ID'    :      0,
            'FSTEM1': fstem1,
            'FSTEM2': fstem2,
            'NUCL'  :   nucl,
            'MULT'  :   None,
            'TYPE'  :    'N'}

def StemLPair(model, nucl1, nucl2, fstem1, fstem2, link):

    nucl12  = [nucl1, nucl2]
    stem12  = ['\\N', '\\N']
    ostem12 = ['\\N', '\\N']
    fstem12 = [fstem1, fstem2]

    for i in (0,1):

        res = model.chains[nucl12[i][0].split('.')[0]][nucl12[i][1]][nucl12[i][2]]

        if res['WING']:

            stem  = model.wings['LU'][res['WING']-1]['STEM']

            if model.stems[stem-1]['FULLSTEM'] == fstem12[i]:     stem12[i]  = stem

        if res['OLDWING']:

            ostem = model.wings['OLD'][res['OLDWING']-1]['STEM']

            if model.oldstems[ostem-1]['FULLSTEM'] == fstem12[i]: ostem12[i] = ostem

    return {'ID'     : 0,
            'LINK'   : link,
            'FSTEM1' : fstem1,
            'FSTEM2' : fstem2,
            'OSTEM1' : ostem12[0],
            'OSTEM2' : ostem12[1],
            'STEM1'  : stem12[0],
            'STEM2'  : stem12[1],
            'MULT'   : None,
            'NUCL1'  : nucl1[0],
            'NUCL2'  : nucl2[0],
            'TYPE'   : 'L'}

def NuclMult(model, bps, ID):

    seen = []
    nts  = []
    wcwb = 0

    for bp_id in bps:

        model.bpairs[bp_id-1]['NUCLMULT'] = ID # BasePairs.nuclmult

        if model.bpairs[bp_id-1]['TYPE'] in ('WC','WB'): wcwb += 1

        for nucl in (model.bpairs[bp_id-1]['NUCL1'], model.bpairs[bp_id-1]['NUCL2']):

            if nucl[0] not in seen:

                seen.append(nucl[0])
                nts.append(nucl)

                model.chains[nucl[0].split('.')[0]][nucl[1]][nucl[2]]['MULT'] = ID # nucls.mult

    Type = Graph.Build_Type([model.bpairs[bp-1] for bp in bps], 'NUCL')

    return {'ID'   :   ID,
            'NTS'  :  nts,   # list of [dssr, where, index] of nucleotides
            'BPS'  :  bps,   # list of bp ids
            'WCWB' : wcwb,   # number of WC/WB BasePairs
            'TYPE' : Type}

def StemMult(model, stempairs, cc, ID):

    seens  = []     # for fstems
    seenn  = []     # for shared nucleotides
    seenl  = []     # for links
    seenp  = []     # for stempairs  (for edges actually)
    seenpL = []     # for stemLpairs (for ledges actually)
    seenpN = []     # for stemNpairs (for nedges actually)
    fstems = []     # list of fstems
    nts    = []     # list of shared nucleotides
    links  = []     # list of links
    edges  = 0      # number of edges
    ledges = 0      # number of link-edges
    nedges = 0      # number of nucl-edges
    cc_edges = []
    
    for p_id in cc:

        stempairs[p_id-1]['MULT'] = ID # StemNPairs.mult or StemLPairs.mult

        pair = stempairs[p_id-1]

        fstem1, fstem2 = pair['FSTEM1'], pair['FSTEM2']

        if str(fstem1)+str(fstem2) not in seenp and\
           str(fstem2)+str(fstem1) not in seenp:

            seenp.append(str(fstem1)+str(fstem2))
            edges += 1
            cc_edges.append(pair) # for Type

        for fstem in (fstem1, fstem2):

            if fstem not in seens:

                seens.append(fstem)
                fstems.append(fstem)

                model.fullstems[fstem-1]['MULT'] = ID # StemsFull.mult

        if pair['TYPE'] == 'N':

            if pair['NUCL'] not in seenn:

                seenn.append(pair['NUCL'])
                nts.append(pair['NUCL'])

            if str(fstem1)+str(fstem2) not in seenpN and\
               str(fstem2)+str(fstem1) not in seenpN:

                seenpN.append(str(fstem1)+str(fstem2))
                nedges += 1

        if pair['TYPE'] == 'L':

            if pair['LINK'] not in seenl:

                seenl.append(pair['LINK'])
                links.append(pair['LINK'])

            if str(fstem1)+str(fstem2) not in seenpL and\
               str(fstem2)+str(fstem1) not in seenpL:

                seenpL.append(str(fstem1)+str(fstem2))
                ledges += 1

            for stem  in (pair['STEM1'],  pair['STEM2']):

                if stem != '\\N' : model.stems[stem-1]['MULT'] = ID     #stems.mult

            for ostem in (pair['OSTEM1'], pair['OSTEM2']):

                if ostem != '\\N': model.oldstems[ostem-1]['MULT'] = ID # StemsOld.mult

    Type = Graph.Build_Type(cc_edges, 'FSTEM')

    return {'ID'     :     ID,
            'FSTEMS' : fstems,    # list of fullstem ids
            'NTS'    :    nts,    # list of shared nucleotides
            'LINKS'  :  links,    # list of links
            'NEDGES' : nedges,    # number of distinct edges from Npairs (distinct fstem1,fstem2)
            'LEDGES' : ledges,    # number of distinct edges from Lpairs (distinct fstem1,fstem2)
            'EDGES'  :  edges,    # number of distinct edges from all pairs (distinct fstem1,fstem2)
            'NPAIRS' :     [],    # list of stemNpair ids (will be filled in add_pairs()
            'LPAIRS' :     [],    # list of stemLpair ids (will be filled in add_pairs()
            'TYPE'   :   Type}

def NuclMults(model):

    nuclmults, ID = [], 1 
    
    # ccs contains lists of bp-ids, one element = edges of one nuclmult
    ccs = Graph.ConnectedComponents(model.bpairs, 'NUCL')

    for cc in ccs:

        if len(cc) > 1:

            nuclmults.append(NuclMult(model,cc,ID))
            ID += 1

    return nuclmults

def nucl_stems_dict(model, Type):

    nucl_stems = {}

    stems = {'':model.stems,'OLD':model.oldstems,'FULL':model.fullstems}[Type]

    for s in stems: 

        for bp in s['PAIRS']:

            nucl1, nucl2 = model.bpairs[bp-1]['NUCL1'][0], model.bpairs[bp-1]['NUCL2'][0]

            if nucl1 not in nucl_stems: nucl_stems[nucl1] = [s['ID'],]
            else:                       nucl_stems[nucl1].append(s['ID'])

            if nucl2 not in nucl_stems: nucl_stems[nucl2] = [s['ID'],]
            else:                       nucl_stems[nucl2].append(s['ID'])

    return nucl_stems

def StemPairs(model):     # StemNPairs and StemLPairs

    stempairs   = []
    nucl_fstems = nucl_stems_dict(model, 'FULL')    # source for fstempairs

    # Nucleotide Pairs
    for nucl in nucl_fstems:

        fstems = nucl_fstems[nucl]

        if len(fstems) > 1:

            for i in range(len(fstems)):
                for j in range(i+1, len(fstems)):

                    stempairs.append(StemNPair(fstems[i],fstems[j],nucl))
    # Link Pairs

    for bp in model.bpairs:

        if bp['FULLSTEM']: continue     # if it's a part of a fullstem

        if bp['NUCL1'][0] in nucl_fstems and bp['NUCL2'][0] in nucl_fstems: # if it's link between fullstems

            wing1 = model.chains[bp['CHAIN1']][bp['NUCL1'][1]][bp['NUCL1'][2]]['WING']
            wing2 = model.chains[bp['CHAIN2']][bp['NUCL2'][1]][bp['NUCL2'][2]]['WING']

            #### for case if this link between two stems which are in one fullstem
            if wing1 and wing2: 

                wing1 = model.wings['LU'][wing1-1]
                wing2 = model.wings['LU'][wing2-1]

                if wing1['STEM'] != wing2['STEM'] and\
               model.stems[wing1['STEM']-1]['FULLSTEM'] == model.stems[wing2['STEM']-1]['FULLSTEM']:

                    fs = model.stems[wing1['STEM']-1]['FULLSTEM']
                    stempairs.append(StemLPair(model,bp['NUCL1'],bp['NUCL2'],fs,fs,bp['LINK']))
            ######################################################################

            for fs1 in nucl_fstems[bp['NUCL1'][0]]:
                for fs2 in nucl_fstems[bp['NUCL2'][0]]:
                    if fs1 != fs2:  # if it's not the same fullstem

                        stempairs.append(StemLPair(model,bp['NUCL1'],bp['NUCL2'],fs1,fs2,bp['LINK']))

    for_sort = {'L':'NUCL1','N':'NUCL'}
    stempairs.sort(key = lambda x: [x['FSTEM1'],x['FSTEM2'],x[for_sort[x['TYPE']]]])  # sorting
    for i in range(len(stempairs)): stempairs[i]['ID'] = i+1    # marking ID

    return stempairs

def add_pairs(model, stemmults):

    for Type in ('N','L'):

        for pair in model.stempairs[Type]:

            stemmults[pair['MULT']-1][Type+'PAIRS'].append(pair['ID'])  # filling stemmult['NPAIRS'/'LPAIRS']

def StemMults(model):

    stemmults, ID = [], 1  
    stempairs     = StemPairs(model)
    ccs           = Graph.ConnectedComponents(stempairs, 'FSTEM')
    
    for cc in ccs:

        stemmults.append(StemMult(model,stempairs,cc,ID))
        ID += 1

    model.stempairs = {}    

    for Type in ('N','L'):      # Making resulting stemNpairs and stemLpairs

        model.stempairs[Type] = sorted([sp for sp in stempairs if sp['TYPE'] == Type],key = lambda x: x['ID'])
        for i in range(len(model.stempairs[Type])): model.stempairs[Type][i]['ID'] = i+1 # marking ID

    add_pairs(model,stemmults)
    
    return stemmults

def add(model):

    nuclmults = NuclMults(model)
    stemmults = StemMults(model)

    model.nuclmults = nuclmults
    model.stemmults = stemmults
    
