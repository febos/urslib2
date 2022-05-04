

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
    

def add(model): 

    biebwe           = BIE_BWE(model)
    #helicalstacking = HelicalStacking(model)

    
    model.biebwe               = biebwe
    #model.helicalstacking     = helicalstacking
