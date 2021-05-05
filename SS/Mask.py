
def FullMask(model):

    bpairs = model.bpairs

    mask = [0]*len(bpairs)

    stem_id = 1

    for i in range(len(bpairs)-1):

        if mask[i]: continue # if already "stemmed"

        in_stem = False

        c = i   # current bp
        j = i+1 # next bp

        while j < len(bpairs):

            nucl1 = [bpairs[c]['CHAIN1'], bpairs[c]['NUCL1'][1], bpairs[c]['NUCL1'][2]] #chain, RES/LIGANDS, index
            nucl2 = [bpairs[c]['CHAIN2'], bpairs[c]['NUCL2'][1], bpairs[c]['NUCL2'][2]]
            nucl3 = [bpairs[j]['CHAIN1'], bpairs[j]['NUCL1'][1], bpairs[j]['NUCL1'][2]]
            nucl4 = [bpairs[j]['CHAIN2'], bpairs[j]['NUCL2'][1], bpairs[j]['NUCL2'][2]]

            if nucl1[:2] == nucl3[:2] and nucl3[2] <= nucl1[2] + 1: #preliminary conditions

                if nucl3[2]  == nucl1[2] + 1 and nucl2[:2] == nucl4[:2] and nucl4[2]  == nucl2[2] - 1 : # necessary conditions

                    if not in_stem:
                        in_stem = True
                        mask[c] = stem_id # stemming current bp

                    mask[j] = stem_id     # stemming next bp
                    c = j
            else:

                if in_stem: stem_id += 1 # if end of stem
                in_stem = False
                break
            
            j += 1

            if not j < len(bpairs): # end of bpairs

                if in_stem: stem_id += 1
                in_stem = False
                break

        if in_stem: stem_id += 1  # if cycle ended but else condition didn't work due to stem continued until end of bpairs

    return mask

def renum(mask): # 0,0,2,2,3,3,0,5,5 -> 0,0,1,1,2,2,0,3,3

    seen = []

    for s in mask:

        if s and s not in seen: seen.append(s)

    renum_dict = {0:0,}

    for i in range(len(seen)) : renum_dict[seen[i]] = i+1

    for i in range(len(mask)): mask[i] = renum_dict[mask[i]]


def OldMask(model,fullmask):

    bpairs = model.bpairs

    mask = [0]*len(bpairs)

    for i in range(len(fullmask)): mask[i] = fullmask[i]

    Max = max(fullmask+[0])

    if not Max: return mask

    for i in range(1,Max+1): # searching each fullstem for non WC/WB edges

        pairs  = []

        for j in range(len(fullmask)):

            if fullmask[j] == i: pairs.append(j)

        side = 0 # can be 0 or -1 and mean the direction (start->end or reversed)

        while True:

            while pairs:

                if bpairs[pairs[side]]['TYPE'] in ('WC','WB'): break # if we found good edge

                mask[pairs[side]] = 0

                if side: pairs = pairs[:side]
                else:    pairs = pairs[1:]

            if len(pairs) == 0: break

            if len(pairs) == 1:

                mask[pairs[0]] = 0
                break

            if side: break
            else:    side = -1

    renum(mask)

    return mask

def LuMask(model):

    bpairs = model.bpairs

    mask = [0]*len(bpairs)

    stem_id = 1
    in_stem = False

    for i in range(len(bpairs)-1):

        if mask[i] or bpairs[i]['TYPE'] not in ('WC','WB'): continue

        c = i   # current
        j = i+1 # next
    
        while j < len(bpairs):

            nucl1 = [bpairs[c]['CHAIN1'], bpairs[c]['NUCL1'][1], bpairs[c]['NUCL1'][2]]
            nucl2 = [bpairs[c]['CHAIN2'], bpairs[c]['NUCL2'][1], bpairs[c]['NUCL2'][2]]
            nucl3 = [bpairs[j]['CHAIN1'], bpairs[j]['NUCL1'][1], bpairs[j]['NUCL1'][2]]
            nucl4 = [bpairs[j]['CHAIN2'], bpairs[j]['NUCL2'][1], bpairs[j]['NUCL2'][2]]

            if nucl1[:2] == nucl3[:2] and nucl3[2] <= nucl1[2] + 1:

                if nucl3[2]  == nucl1[2] + 1 and nucl2[:2] == nucl4[:2] and nucl4[2]  == nucl2[2] - 1 and\
                   bpairs[j]['TYPE'] in ('WC','WB'):

                    if not in_stem:
                        in_stem = True
                        mask[c] = stem_id

                    mask[j] = stem_id
                    c = j
            else:

                if in_stem: stem_id += 1
                in_stem = False
                break

            j += 1

            if not j < len(bpairs): # end of bpairs

                if in_stem: stem_id += 1
                in_stem = False
                break

    return mask

def RevMask(model):

    bpairs = model.bpairs

    mask = [0]*len(bpairs)

    stem_id = 1
    in_stem = False

    for i in range(len(bpairs)-1):

        if mask[i]: continue

        c = i   # current
        j = i+1 # next
    
        while j < len(bpairs):

            nucl1 = [bpairs[c]['CHAIN1'], bpairs[c]['NUCL1'][1], bpairs[c]['NUCL1'][2]]
            nucl2 = [bpairs[c]['CHAIN2'], bpairs[c]['NUCL2'][1], bpairs[c]['NUCL2'][2]]
            nucl3 = [bpairs[j]['CHAIN1'], bpairs[j]['NUCL1'][1], bpairs[j]['NUCL1'][2]]
            nucl4 = [bpairs[j]['CHAIN2'], bpairs[j]['NUCL2'][1], bpairs[j]['NUCL2'][2]]

            if nucl1[:2] == nucl3[:2] and nucl3[2] <= nucl1[2] + 1:

                if nucl3[2]  == nucl1[2] + 1 and nucl2[:2] == nucl4[:2] and nucl4[2]  == nucl2[2] + 1:

                    if not in_stem:
                        in_stem = True
                        mask[c] = stem_id

                    mask[j] = stem_id
                    c = j

            else: # end of stem
 
                if in_stem: stem_id += 1
                in_stem = False
                break
      
            j += 1

            if not j < len(bpairs): # end of bpairs

                if in_stem: stem_id += 1
                in_stem = False
                break

    return mask

def add(model):

    fullmask = FullMask(model)
    oldmask  = OldMask(model,fullmask)
    mask     = LuMask(model)
    revmask  = RevMask(model)

    model.mask     = mask
    model.fullmask = fullmask
    model.oldmask  = oldmask
    model.revmask  = revmask
