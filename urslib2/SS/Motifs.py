

def BIE_BWE(model):

    stackparam = {}

    for x in model.non_pairs:
        if x['NUCL1'] not in stackparam:
            stackparam[x['NUCL1']] = {}
        stackparam[x['NUCL1']][x['NUCL2']] = x['STACKING']
        if x['NUCL2'] not in stackparam:
            stackparam[x['NUCL2']] = {}
        if 'forward' in x['STACKING']:
            stackparam[x['NUCL2']][x['NUCL1']] = x['STACKING'].replace('pm(>>,forward)','mp(<<,backward)')
        elif 'backward' in x['STACKING']:
            stackparam[x['NUCL2']][x['NUCL1']] = x['STACKING'].replace('mp(<<,backward)','pm(>>,forward)')
        else:
            stackparam[x['NUCL2']][x['NUCL1']] = x['STACKING']

    stacks = []

    for s in model.stacks:
        if s[0]<3:
            continue

        for i in range(1,s[0]-1):
            thr = s[1][i-1:i+2]
            thrs = [thr,thr[::-1]]

            if not all([x in model.dssrnucls for x in thr]):
                continue
            
            for thr in thrs:
                if model.dssrnucls[thr[0]][0]==model.dssrnucls[thr[2]][0] and\
                   (model.dssrnucls[thr[0]][2]+1==model.dssrnucls[thr[2]][2] or\
                    model.dssrnucls[thr[0]][2]+2==model.dssrnucls[thr[2]][2]) and\
                   (model.dssrnucls[thr[0]][0]!=model.dssrnucls[thr[1]][0] or \
                    model.dssrnucls[thr[0]][2]+1!=model.dssrnucls[thr[1]][2]):
                    stacks.append(thr)

    good = []

    for st in stacks:
        stack = {}
        if   model.dssrnucls[st[0]][2]+1==model.dssrnucls[st[2]][2]: type1 = 'BIE'
        elif model.dssrnucls[st[0]][2]+2==model.dssrnucls[st[2]][2]: type1 = 'BWE'
        stack['TYPE'] = type1
        stack['NUCLS'] = st
        stack['STACK1'] = stackparam[st[0]][st[1]].split('--')
        stack['STACK2'] = stackparam[st[1]][st[2]].split('--')
                
        good.append(stack)

    # sorted by chain_length, chain_id, firstNucl_index 
    return sorted(good,key= lambda x: (model.chains[model.dssrnucls[x['NUCLS'][0]][0]]['LENGTH'],
                                       model.dssrnucls[x['NUCLS'][0]][0],
                                       model.dssrnucls[x['NUCLS'][0]][2]))

    
def HelicalStacking(model):

    stackparam = {}

    for x in model.non_pairs:

        if x['NUCL1'] not in stackparam:
            stackparam[x['NUCL1']] = {}
        stackparam[x['NUCL1']][x['NUCL2']] = x['STACKING']
        
        if x['NUCL2'] not in stackparam:
            stackparam[x['NUCL2']] = {}
        if 'forward' in x['STACKING']: 
            stackparam[x['NUCL2']][x['NUCL1']] = x['STACKING'].replace('pm(>>,forward)','mp(<<,backward)')
        elif 'backward' in x['STACKING']: 
            stackparam[x['NUCL2']][x['NUCL1']] = x['STACKING'].replace('mp(<<,backward)','pm(>>,forward)')
        else:
            stackparam[x['NUCL2']][x['NUCL1']] = x['STACKING']
        stackparam[x['NUCL2']][x['NUCL1']] = [stackparam[x['NUCL2']][x['NUCL1']], x['MINDIST'], x['ANGLE']]
        stackparam[x['NUCL1']][x['NUCL2']] = [stackparam[x['NUCL1']][x['NUCL2']], x['MINDIST'], x['ANGLE']]
    
    cstacks = []
    
    for y in model.helices:
    
        stems = []
    
        for bpid in y['PAIRS']:
            st = model.bpairs[bpid-1]['STEM']
            if st and st not in stems:
                stems.append(st)
    
        for j in range(1,len(stems)):
        
            bp11 = [model.bpairs[model.stems[stems[j-1]-1]['PAIRS'][0]-1]['NUCL1'][0],
                    model.bpairs[model.stems[stems[j-1]-1]['PAIRS'][0]-1]['NUCL2'][0]]
            bp12 = [model.bpairs[model.stems[stems[j-1]-1]['PAIRS'][-1]-1]['NUCL1'][0],
                    model.bpairs[model.stems[stems[j-1]-1]['PAIRS'][-1]-1]['NUCL2'][0]]
        
            bp21 = [model.bpairs[model.stems[stems[j]-1]['PAIRS'][0]-1]['NUCL1'][0],
                    model.bpairs[model.stems[stems[j]-1]['PAIRS'][0]-1]['NUCL2'][0]]
            bp22 = [model.bpairs[model.stems[stems[j]-1]['PAIRS'][-1]-1]['NUCL1'][0],
                    model.bpairs[model.stems[stems[j]-1]['PAIRS'][-1]-1]['NUCL2'][0]]
       
            for bp1 in (bp11,bp12):
                for bp2 in (bp21,bp22):

                    stacks = []

                    for nucl1 in bp1:
                        for nucl2 in bp2:
                            if nucl1 in stackparam and nucl2 in stackparam[nucl1]\
                               and stackparam[nucl1][nucl2][0] != '\\N':
                                stacks.append([nucl1,nucl2,*stackparam[nucl1][nucl2]])
                    if stacks:
                        cstacks.append({'BP1':bp1,
                                        'BP2':bp2,
                                        'STACKING':stacks})   
    return cstacks

def NN_Platform(model):

    plats = []

    for bp in model.bpairs:
    
        ch1,pl1,i1 = model.dssrnucls[bp['NUCL1'][0]]
        ch2,pl2,i2 = model.dssrnucls[bp['NUCL2'][0]]

        if ch1==ch2 and pl1==pl2 and abs(i1-i2)==1:
        
            plats.append({'NUCL1':bp['NUCL1'][0],
                          'NUCL2':bp['NUCL2'][0],
                          'BPID' :bp['ID']})
    return plats


def add(model): 

    biebwe           = BIE_BWE(model)
    helicalstacking  = HelicalStacking(model)
    nnplatform       = NN_Platform(model)

    
    model.biebwe              = biebwe
    model.helicalstacking     = helicalstacking
    model.nnplatform          = nnplatform





    
