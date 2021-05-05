
def change_desc(desc,wts): #converting desc from useful to understandable

    desc2 = ''

    if desc[0] in ('Z','T'): desc2 += str(wts[0]['LEN'])

    i = 1

    while i < len(desc):

        if desc[i] in ('Z','T'):

            if desc[i-1] in ('Z','T'): desc2 += '--' # for break between different chains
            elif desc[i-1] == 'W'    : desc2 += '.'  # . will be between wing and thread

            desc2 += str(wts[i]['LEN'])

        elif desc[i] == 'W':

            desc2 += '.w'+str(wts[i]['LEN'])

        else:

            desc2 += '-|'+desc[i]+'|-' # for stems and blocks
            i += 1

        i += 1

    return desc2

def elements(preloop,looptype,ID): # for rss.threadloop, rss.wingloop and rss.faceloop

    threads = []
    wings   = []
    faces   = []
    num     = 1
    wnum    = 1
    fnum    = 1
    stem1   = 0
    side    = 1

    tjmol  = []
    sfjmol = []
    bfjmol = []
    wjmol  = []

    flag = False
    
    for i in range(len(preloop[1])):

        if preloop[2][i] in ('B','S'):      # faces

            if preloop[2][i] == 'B': bfjmol.append(preloop[1][i]['JMOL'])
            else                   : sfjmol.append(preloop[1][i]['JMOL'])

            if not flag:    # start of face

                side += 1
                stem1 = preloop[1][i]['STEM']  
                flag  = True

            else:           # end of face

                faces.append({'MODEL'    :                      1,
                              'LOOPID'   :                     ID,
                              'LOOPTYPE' :               looptype,
                              'STEM1'    :                  stem1,
                              'STEM2'    :  preloop[1][i]['STEM'],
                              'TYPE'     :          preloop[2][i],
                              'NUM'      :                  fnum})
                fnum += 1
                flag  = False

        elif preloop[2][i] in ('Z','T'):    # threads

            if preloop[2][i] == 'T': tjmol.append(preloop[1][i]['JMOL'])

            threads.append({'MODEL'    :                      1,
                            'LOOPID'   :                     ID,
                            'LOOPTYPE' :               looptype,
                            'THREAD'   :    preloop[1][i]['ID'],
                            'NUM'      :                    num,
                            'SIDE'     :                   side,
                            'LINKS'    : preloop[1][i]['LINKS'],
                            'SEQ'      :   preloop[1][i]['SEQ'],
                            'LEN'      :  preloop[1][i]['LEN']})
            num += 1

        else:                               # wings

            wjmol.append(preloop[1][i]['JMOL'])

            wings.append({'MODEL'    :                     1,
                          'LOOPID'   :                    ID,
                          'LOOPTYPE' :              looptype,
                          'WING'     :   preloop[1][i]['ID'],
                          'NUM'      :                  wnum,
                          'SIDE'     :                  side,
                          'SEQ'      :  preloop[1][i]['SEQ'],
                          'LEN'      : preloop[1][i]['LEN']})

            wnum += 1


    tjmol, sfjmol, bfjmol, wjmol = ','.join(tjmol), ','.join(sfjmol), ','.join(bfjmol), ','.join(wjmol)

    return threads, wings, faces, tjmol, sfjmol, bfjmol, wjmol 

def Bulge(preloop, ID):

    desc = preloop[2]
    wts  = preloop[1] 

    ch1  = wts[0]['CHAIN']
    ch2  = wts[-1]['CHAIN']
    ws   = desc.count('W')

    if preloop[0]!='\\N':

        if 'SS' in desc: place = desc.find('SS')
        else:            place = desc.find('BB')

        if   desc[:place]   == 'Z': side = 'R'
        elif desc[place+2:] == 'Z': side = 'L'

    else:
        if 'SSZSS' in desc.replace('BB','SS'): side = 'R'
        else                                 : side = 'L'

    Len    = sum([wts[i]['LEN']   for i in range(len(wts)) if desc[i] in ('T','W')])
    TLen   = sum([wts[i]['LEN']   for i in range(len(wts)) if desc[i] == 'T'])
    Links  = sum([wts[i]['LINKS'] for i in range(len(wts)) if desc[i] == 'T'])

    breaks = int(preloop[0]=='\\N') + desc.count('TT') + \
             desc.count('TZ') + desc.count('ZT') + desc.count('ZZ')

    tloop, wloop, floop, tjmol, sfjmol, bfjmol, wjmol = elements(preloop,'B',ID)
    desc    = change_desc(preloop[2],preloop[1])  

    return {'ID'     :         ID,
            'STEM'   : preloop[0],
            'MODEL'  :          1,
            'CHAIN1' :        ch1,
            'CHAIN2' :        ch2,
            'WINGS'  :         ws,
            'SIDE'   :       side,
            'TYPE'   :       desc,
            'PTYPE'  :         '',
            'LEN'    :        Len,
            'TLEN'   :       TLen,
            'LINKS'  :      Links,
            'TLOOP'  :      tloop,
            'WLOOP'  :      wloop,
            'FLOOP'  :      floop,
            'BREAKS' :     breaks,
            'SJMOL'  :         '',
            'TJMOL'  :      tjmol,
            'SFJMOL' :     sfjmol,
            'BFJMOL' :     bfjmol,
            'WJMOL'  :      wjmol,
            'MISS'   :      False}

def Hairpin(preloop, ID):

    desc    = preloop[2]
    wts     = preloop[1] 

    ch      = wts[0]['CHAIN']
    ws      = desc.count('W')

    Len     = sum([wts[i]['LEN']   for i in range(len(wts)) if desc[i] in ('T','W')])
    TLen    = sum([wts[i]['LEN']   for i in range(len(wts)) if desc[i] == 'T'])
    Links   = sum([wts[i]['LINKS'] for i in range(len(wts)) if desc[i] == 'T'])

    breaks  = int(preloop[0]=='\\N') + desc.count('TT') + \
              desc.count('TZ') + desc.count('ZT') + desc.count('ZZ')
    
    tloop, wloop, floop, tjmol, sfjmol, bfjmol, wjmol = elements(preloop,'H',ID)
    desc    = change_desc(preloop[2],preloop[1])  

    return {'ID'     :         ID,
            'STEM'   : preloop[0],
            'MODEL'  :          1,
            'CHAIN'  :         ch,
            'WINGS'  :         ws,
            'TYPE'   :       desc,
            'PTYPE'  :         '',
            'LEN'    :        Len,
            'TLEN'   :       TLen,
            'LINKS'  :      Links,
            'TLOOP'  :      tloop,
            'WLOOP'  :      wloop,
            'FLOOP'  :      floop,
            'BREAKS' :     breaks,
            'SJMOL'  :         '',
            'TJMOL'  :      tjmol,
            'SFJMOL' :     sfjmol,
            'BFJMOL' :     bfjmol,
            'WJMOL'  :      wjmol,
            'MISS'   :      False}

def Internal(preloop, ID):

    desc = preloop[2]
    wts  = preloop[1] 

    ch1  = wts[0]['CHAIN']
    ch2  = wts[-1]['CHAIN']
    ws   = desc.count('W')

    if preloop[0]!='\\N':

        if 'SS' in desc: place = desc.find('SS')
        else:            place = desc.find('BB')

        Len1    = sum([wts[i]['LEN']   for i in range(0,place)])
        Len2    = sum([wts[i]['LEN']   for i in range(place+2,len(wts))])
    else:
        look1 = desc.replace('BB','SS').find('SS')
        look2 = desc.replace('BB','SS').rfind('SS')

        Len1    = sum([wts[i]['LEN']   for i in range(0,look1)]) + sum([wts[i]['LEN'] for i in range(look2+2,len(desc))])
        Len2    = sum([wts[i]['LEN']   for i in range(look1+2,look2)])

    Links   = sum([wts[i]['LINKS'] for i in range(len(wts)) if desc[i] == 'T'])
    TLen    = sum([wts[i]['LEN']   for i in range(len(wts)) if desc[i] == 'T'])

    breaks  = int(preloop[0]=='\\N') + desc.count('TT') + \
              desc.count('TZ') + desc.count('ZT') + desc.count('ZZ')
    
    tloop, wloop, floop, tjmol, sfjmol, bfjmol, wjmol = elements(preloop,'I',ID)
    desc    = change_desc(preloop[2],preloop[1])  

    return {'ID'     :           ID,
            'STEM'   :   preloop[0],
            'MODEL'  :            1,
            'CHAIN1' :          ch1,
            'CHAIN2' :          ch2,
            'WINGS'  :           ws,
            'SYM'    : Len1 == Len2,
            'TYPE'   :         desc,
            'PTYPE'  :           '',
            'LEN'    :    Len1+Len2,
            'TLEN'   :         TLen,
            'LINKS'  :        Links,
            'TLOOP'  :        tloop,
            'WLOOP'  :        wloop,
            'FLOOP'  :        floop,
            'BREAKS' :       breaks,
            'SJMOL'  :           '',
            'TJMOL'  :        tjmol,
            'SFJMOL' :       sfjmol,
            'BFJMOL' :       bfjmol,
            'WJMOL'  :        wjmol,
            'MISS'   :        False}

def Junction(preloop, ID):

    desc   = preloop[2]
    wts    = preloop[1] 

    ch1    = wts[0]['CHAIN']
    ch2    = wts[-1]['CHAIN']
    ws     = desc.count('W')
    ths    = desc.count('T')  + desc.count('Z')
    sfs    = desc.count('SS')
    bfs    = desc.count('BB')
    sds    = bfs + sfs + 1 - int(preloop[0]=='\\N')

    Len    = sum([wts[i]['LEN']   for i in range(len(wts)) if desc[i] in ('T','W')])
    TLen   = sum([wts[i]['LEN']   for i in range(len(wts)) if desc[i] == 'T'])
    Links  = sum([wts[i]['LINKS'] for i in range(len(wts)) if desc[i] == 'T'])

    breaks = int(preloop[0]=='\\N') + desc.count('TT') + \
             desc.count('TZ') + desc.count('ZT') + desc.count('ZZ')

    tloop, wloop, floop, tjmol, sfjmol, bfjmol, wjmol = elements(preloop,'J',ID)
    desc    = change_desc(preloop[2],preloop[1])  

    return {'ID'        :         ID,
            'STEM'      : preloop[0],
            'MODEL'     :          1,
            'THREADSNUM':        ths,
            'CHAIN1'    :        ch1,
            'CHAIN2'    :        ch2,
            'WINGS'     :         ws,
            'SIDES'     :        sds,
            'SFACES'    :        sfs,
            'BFACES'    :        bfs,
            'TYPE'      :       desc,
            'PTYPE'     :         '',
            'LEN'       :        Len,
            'TLEN'      :       TLen,
            'LINKS'     :      Links,
            'TLOOP'     :      tloop,
            'WLOOP'     :      wloop,
            'FLOOP'     :      floop,
            'BREAKS'    :     breaks,
            'SJMOL'     :         '',
            'TJMOL'     :      tjmol,
            'SFJMOL'    :     sfjmol,
            'BFJMOL'    :     bfjmol,
            'WJMOL'     :      wjmol,
            'MISS'      :      False}

# let's say it's an EXTERIOR PRE-LOOP, containing needed wings and threads
def Exteriors(model):

    ext_preloops = {}

    thred_len = {False:'Z',True:'T'}

    wts  = [] #wings and threads
    desc = ''

    wings = model.wings['LU']
    if not wings: return ext_preloops 

    level = 0
    in_loop = False

    for w in wings:

        if not in_loop:

            wts.append(model.threads[w['PREV']-1])
            wts.append(w)
            desc += thred_len[bool(model.threads[w['PREV']-1]['LEN'])]
            level += 1
            next_stem = w['STEM'] 
            in_loop = True

        else:

            if w['TYPE'] == 'L': level += 1
            else               : level -= 1

            if level == 0:

                wts.append(w)
                if wts[-1]['STEM'] == wts[-2]['STEM']: desc += 'SS'
                else: desc += 'BB'
                wts.append(model.threads[w['NEXT']-1])
                desc += thred_len[bool(model.threads[w['NEXT']-1]['LEN'])]

                if not w['NEXTW']:

                    ext_preloops[next_stem] = ['\\N', wts, desc]
                    in_loop, wts, desc = False, [], ''
                    del next_stem

            if level == 1 and w['TYPE'] == 'L': wts.append(w)
  
    return ext_preloops

# let's say it's a PRE-LOOP, containing needed wings and threads
def Loop(model, stem_id):

    thred_len = {False:'Z',True:'T'}

    wts  = [] #wings and threads
    desc = ''

    wings = model.wings['LU']

    L_id = model.stems[stem_id - 1]['LEFT']
    R_id = model.stems[stem_id - 1]['RIGHT']

    seen_threads = []

    if wings[L_id-1]['NEXT']:

        wts.append(model.threads[wings[L_id-1]['NEXT']-1])
        desc += thred_len[bool(model.threads[wings[L_id-1]['NEXT']-1]['LEN'])]
        seen_threads.append(wings[L_id-1]['NEXT'])

    i = L_id
    #print(i)
    #print(stem_id)
    while wings[i]['STEM'] != stem_id:
        
        # if L wing and his R inside our stem
        if wings[i]['TYPE'] == 'L' and wings[i]['ANOTHER'] < R_id:

            if wings[i]['PREV'] and wings[i]['PREV'] not in seen_threads:

                wts.append(model.threads[wings[i]['PREV']-1])
                desc += thred_len[bool(model.threads[wings[i]['PREV']-1]['LEN'])]

            another = wings[i]['ANOTHER']
            This    = i+1
            j       = i+1
            ### searching for isolated block
            while j != another - 1:

                if wings[j]['TYPE'] == 'L' and \
                   wings[j]['ANOTHER'] > another and\
                   wings[j]['ANOTHER'] < R_id:

                    another = wings[j]['ANOTHER']

                j += 1
            ###
            wts.append(wings[This-1])
            wts.append(wings[another-1])

            if wings[This-1]['STEM'] == wings[another-1]['STEM']: desc += 'SS'
            else:                                                 desc += 'BB'

            if wings[another-1]['NEXT']:

                wts.append(model.threads[wings[another-1]['NEXT']-1])
                desc += thred_len[bool(model.threads[wings[another-1]['NEXT']-1]['LEN'])]
                seen_threads.append(wings[another-1]['NEXT'])

            
            i = another 

        # if (R wing) or (L wing with his R outside our stem)
        else:

            if wings[i]['PREV'] and wings[i]['PREV'] not in seen_threads:

                wts.append(model.threads[wings[i]['PREV']-1])
                desc += thred_len[bool(model.threads[wings[i]['PREV']-1]['LEN'])]

            wts.append(wings[i])
            desc += 'W'

            if wings[i]['NEXT']:

                wts.append(model.threads[wings[i]['NEXT']-1])
                desc += thred_len[bool(model.threads[wings[i]['NEXT']-1]['LEN'])]
                seen_threads.append(wings[i]['NEXT'])

            i += 1
            

    if wings[R_id-1]['PREV'] and wings[R_id-1]['PREV'] not in seen_threads:

        wts.append(model.threads[wings[R_id-1]['PREV']-1])
        desc += thred_len[bool(model.threads[wings[R_id-1]['PREV']-1]['LEN'])]

    return [stem_id, wts, desc]

def classify(desc, for_nonloop = False):

    bfaces  = desc.count('BB')
    sfaces  = desc.count('SS')
    pseudo  = 'W'  in desc

    if   pseudo: ptype = 'P'
    elif bfaces: ptype = 'I'
    else       : ptype = 'C'

    if   sfaces + bfaces == 0:                           return ['Hairpin', ptype]

    elif sfaces + bfaces == 1:

        if sfaces: look = 'SS'
        else:      look = 'BB'

        place = desc.find(look)

        if desc[:place] == 'Z' or desc[place+2:] == 'Z': return ['Bulge',    ptype]
        else:                                            return ['Internal', ptype]
    else:                                                return ['Junction', ptype]

def classifyExt(desc, for_nonloop = False):

    bfaces  = desc.count('BB')
    sfaces  = desc.count('SS')
    pseudo  = 'W'  in desc

    if   pseudo: ptype = 'P'
    elif bfaces: ptype = 'I'
    else       : ptype = 'C'

    if   sfaces + bfaces <= 1:                                      return ['Hairpin', ptype]

    elif sfaces + bfaces == 2:

        dd = desc.replace('BB','SS')

        look1 = dd.find('SS')
        look2 = dd.rfind('SS')

        if dd[:look1]+dd[look2+2:]=='ZZ' or dd[look1+2:look2]=='Z': return ['Bulge',    ptype]
        else:                                                       return ['Internal', ptype]
    else:                                                           return ['Junction', ptype]

def Loops(model):

    loops    = {'HAIRPIN' :[],
                'BULGE'   :[],
                'INTERNAL':[],
                'JUNCTION':[]}

    functions = {'BULGE'   :    Bulge,
                 'HAIRPIN' :  Hairpin,
                 'INTERNAL': Internal,
                 'JUNCTION': Junction}
    
    preloops     = []
    ext_preloops = Exteriors(model)  # dict: key = next stem, value = preloop 

    for i in range(1,len(model.stems)+1):

        if i in ext_preloops: preloops.append(ext_preloops[i])
        preloops.append(Loop(model,i))

    ID = {'B':1,'H':1,'I':1,'J':1}

    for preloop in preloops:

        if preloop[0]=='\\N': Types = classifyExt(preloop[2])
        else:                 Types = classify(preloop[2])

        loops[Types[0].upper()].append(functions[Types[0].upper()](preloop,ID[Types[0][0]]))
        loops[Types[0].upper()][-1]['PTYPE']    = Types[1]

        # is missing nts in loop
        for i in [x for x in range(len(preloop[2])) if preloop[2][x] in ('W','T')]:
            start,end = preloop[1][i]['START'],preloop[1][i]['END']
            for nt in model.chains[start[0].split('.')[0]][start[1]][start[2]:end[2]+1]:
                if nt['MISS']:
                    loops[Types[0].upper()][-1]['MISS'] = True
                    break
            if loops[Types[0].upper()][-1]['MISS']: break
        ###

        if loops[Types[0].upper()][-1]['STEM'] != '\\N':
            model.stems[preloop[0]-1]['LOOPTYPE']   = Types[0][0]
            model.stems[preloop[0]-1]['LOOPPSEUDO'] = Types[1]
            loops[Types[0].upper()][-1]['SJMOL']    = model.stems[preloop[0]-1]['JMOL']
            model.stems[preloop[0]-1]['LOOPJMOL']   = ';'.join([loops[Types[0].upper()][-1]['SJMOL'],
                                                                loops[Types[0].upper()][-1]['TJMOL'],
                                                                loops[Types[0].upper()][-1]['SFJMOL'],
                                                                loops[Types[0].upper()][-1]['BFJMOL'],
                                                                loops[Types[0].upper()][-1]['WJMOL']])
        ID[Types[0][0]] += 1
        
    return loops

def add(model):

    loops = Loops(model)

    model.loops = loops
