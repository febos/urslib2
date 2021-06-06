
def Link(bp, ID):

    return {'ID'       :                ID,
            'MODEL'    :                 1,
            'CHAIN1'   :      bp['CHAIN1'],
            'CHAIN2'   :      bp['CHAIN2'],
            'BP'       :          bp['ID'],
            'BPTYPE'   :        bp['TYPE'],
            'LEFTHRD'  :              None,
            'RIGHTHRD' :              None,
            'LEFTWING' :              None,
            'RIGHTWING':              None,
            'TYPE'     :                 0,
            'DEPTH'    :                 0,
            'NUCL1'    :       bp['NUCL1'],
            'NUCL2'    :       bp['NUCL2'],
            'CLASS1'   :    bp['CLASS'][0],
            'CLASS2'   :    bp['CLASS'][1],
            'CLASS3'   :    bp['CLASS'][2],
            'DIST'     :              None,
            'REL'      :              None,
            'SS1'      :              None,
            'SS2'      :              None}
            
def classify(model):

    stem_inside = []    # [chain_order[chain1] + index(nucl1), chain_order[chain2] + index(nucl2)]
    bad_inside  = []    # same as stem_inside but for non WC/WB bps
    wcwb_inside = []    # same as bad_inside  but for WC/WB bps

    chains = sorted(model.chain_order,key= lambda x: model.chain_order[x])

    def StartEnd(model,chains,bp,Type):
        dist = 0
        for ch in chains:
            if ch==bp['CHAIN'+Type]:
                dist += bp['NUCL'+Type][2]+1
                if bp['NUCL'+Type][1]!='RES': dist += model.chains[ch]['LENGTH'] 
                break
            else: dist += model.chains[ch]['LENGTH']
        return dist

    # filling stem_inside
    for stem in model.stems:

        start = StartEnd(model,chains,model.bpairs[stem['PAIRS'][0]-1],'1')
        end   = StartEnd(model,chains,model.bpairs[stem['PAIRS'][0]-1],'2')

        stem_inside.append([min([start,end]), max([start,end])])

    # filling bad_inside and wcwb_inside
    for bp in model.bpairs:

        start = StartEnd(model,chains,bp,'1')
        end   = StartEnd(model,chains,bp,'2')

        if bp['TYPE'] not in ('WC','WB'):

            bad_inside.append([min([start,end]), max([start,end])])

        elif not model.mask[bp['ID']-1]: # WC/WB and not a part of stem

            wcwb_inside.append([min([start,end]), max([start,end])])

    # classifying
    for link in model.links:

        start = StartEnd(model,chains,link,'1')
        end   = StartEnd(model,chains,link,'2')

        if start>end: start,end = end,start

        for inside in stem_inside:
            
            if inside[0] < start < inside[1] and inside[1] < end or\
               start < inside[0] and inside[0] < end < inside[1]:   # if link intersects with stem

                link['TYPE'] = 3
                break

        if not link['TYPE']: # if type is still == 0

            for inside in wcwb_inside: 

                if inside[0] < start < inside[1] and inside[1] < end or\
                   start < inside[0] and inside[0] < end < inside[1]:   # if link intersects with lone WC/WB

                    link['TYPE'] = 2
                    break

        if not link['TYPE']:

            for inside in bad_inside:

                if inside[0] < start < inside[1] and inside[1] < end or\
                    start < inside[0] and inside[0] < end < inside[1]:  # if link intersects with bad bp

                    link['TYPE'] = 1
                    break

        # depth: (now == 0 i.e. link for stem but is a part of oldstem and fullstem)

        if not model.oldmask[link['BP']-1] : link['DEPTH'] = 1 # link is for stem and oldstem but is a part of fullstem  
        if not model.fullmask[link['BP']-1]: link['DEPTH'] = 2 # link is for all types of stems
    

def Links(model):

    chains = sorted(model.chain_order,key= lambda x: model.chain_order[x])

    links = []

    link_id = 1

    for bp in model.bpairs:

        if not bp['STEM']:

            links.append(Link(bp,link_id))
            bp['LINK'] = link_id
            link_id += 1

            dist = {'1':0,'2':0}
            for n in ('1','2'):
                for ch in chains:
                    if ch==links[-1]['CHAIN'+n]:
                        dist[n] += links[-1]['NUCL'+n][2]+1
                        if links[-1]['NUCL'+n][1]!='RES': dist[n] += model.chains[ch]['LENGTH'] 
                        break
                    else: dist[n] += model.chains[ch]['LENGTH']
            links[-1]['DIST'] = abs(dist['2'] - dist['1'] - 1)

    return links

def relation(model):

    # LEFTHRD RIGHTHRD LEFTWING RIGHTWING for links; links for threads
    for link in model.links:

        ch1   ,ch2    = link['CHAIN1']  , link['CHAIN2']
        where1,where2 = link['NUCL1'][1], link['NUCL2'][1] #RES/LIGANDS
        ind1  ,ind2   = link['NUCL1'][2], link['NUCL2'][2] # indexes in chain[RES/LIGANDS]

        link['LEFTHRD']     = model.chains[ch1][where1][ind1]['THREAD']
        link['RIGHTHRD']    = model.chains[ch2][where2][ind2]['THREAD']
        link['LEFTWING']    = model.chains[ch1][where1][ind1]['WING']
        link['RIGHTWING']   = model.chains[ch2][where2][ind2]['WING']

        if link['LEFTHRD'] : model.threads[link['LEFTHRD']-1]['LINKS']  += 1
        if link['RIGHTHRD']: model.threads[link['RIGHTHRD']-1]['LINKS'] += 1

def add(model):

    links = Links(model)

    model.links = links

    classify(model) 
    relation(model)
