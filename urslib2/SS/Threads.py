
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

def ExtThread(model, wing, Type, ID):

    seq_dict = {'RES':'SEQ2', 'LIGANDS':'LIGSEQ'}

    ch = wing['CHAIN']
    where = wing['START'][1]
    
    if Type == 'PREV':

        if wing['START'][2]: end_ind = wing['START'][2] - 1 # if wing['START'][2] != 0
        else:                end_ind = 0

        endnucl = model.chains[ch][where][end_ind]

        prev = None
        Next = wing['ID']

        start = [model.chains[ch][where][0]['DSSR'], where, 0]
        end   = [endnucl['DSSR'], where, end_ind]

        if wing['START'][2]: Len = end[2] - start[2] + 1
        else: Len = 0 # because in zero-len extthreads (end == start) but (end - start + 1) == 1 in this case

    else:

        if wing['END'][2] != len(model.chains[ch][where])-1: start_ind = wing['END'][2]+1 # if it's not the last index
            
        else:                                                start_ind = wing['END'][2]

        startnucl = model.chains[ch][where][start_ind]
        
        prev  = wing['ID']
        Next  = None

        start = [startnucl['DSSR'], where, start_ind]
        end   = [model.chains[ch][where][-1]['DSSR'], where, len(model.chains[ch][where])-1]

        if wing['END'][2] != len(model.chains[ch][where])-1: Len = end[2] - start[2] + 1
        else:                                                Len = 0



    seq = ','.join(model.chains[ch][seq_dict[where]][start[2]:end[2]+1])

    if Len == 0: jmol,seq = '',','
    else: jmol = MakeJmol([x['DSSR'] for x in model.chains[ch][start[1]][start[2]:end[2]+1]])

    return {'ID'   :    ID,
            'MODEL':     1,
            'CHAIN':    ch,
            'PREV' :  prev,
            'NEXT' :  Next,
            'START': start,
            'END'  :   end,
            'LEN'  :   Len,
            'LINKS':     0,
            'SEQ'  :   seq,
            'EXT'  :     1,
            'JMOL' :  jmol,
            'FULL' : False} #thread == chain
            

def Thread(model, wing1, wing2, ID):

    seq_dict = {'RES':'SEQ2', 'LIGANDS':'LIGSEQ'}

    ch = wing1['CHAIN']
    where = wing1['START'][1]

    startnucl = model.chains[ch][where][wing1['END'][2]+1]
    endnucl   = model.chains[ch][where][wing2['START'][2]-1]

    prev = wing1['ID']
    Next = wing2['ID']

    start = [startnucl['DSSR'], where, wing1['END'][2]+1]
    end   = [  endnucl['DSSR'], where, wing2['START'][2]-1]
    
    seq = ','.join(model.chains[ch][seq_dict[where]][start[2]:end[2]+1])
    Len = end[2] - start[2] + 1

    if Len == 0: jmol,seq = '',','
    else: jmol = MakeJmol([x['DSSR'] for x in model.chains[ch][start[1]][start[2]:end[2]+1]])

    return {'ID'   :    ID,
            'MODEL':     1,
            'CHAIN':    ch,
            'PREV' :  prev,
            'NEXT' :  Next,
            'START': start,
            'END'  :   end,
            'LEN'  :   Len,
            'LINKS':     0,
            'SEQ'  :   seq,
            'EXT'  :     0,
            'JMOL' :  jmol,
            'FULL' : False} #thread == chain


def Threads(model):

    threads = []

    thread_id = 1

    wings = model.wings['LU'] #select which wings we are processing (LU/OLD/FULL)

    if wings:
    # handling first wing
    
        #if wings[0]['START'][2]: # if start of first wing != start of chain
        if True:#this allows zero threads
            threads.append(ExtThread(model,wings[0],'PREV',thread_id))
            wings[0]['PREV'] = thread_id
            thread_id += 1
        
        if wings[0]['CHAIN'] != wings[1]['CHAIN']: # if there is break between first and second wing
            
            #if wings[0]['END'][2] != len(model.chains[wings[0]['END'][0][0]]['RES']) - 1: # if end of first wing != end of chain
            if True:#this allows zero threads
                threads.append(ExtThread(model,wings[0],'NEXT',thread_id))
                wings[0]['NEXT'] = thread_id
                thread_id += 1
        
        #elif wings[0]['END'][2] != wings[1]['START'][2] - 1: # elif thread between first and second wing has non-zero length 
        else: #this allows zero threads
                threads.append(Thread(model,wings[0],wings[1],thread_id))
                wings[0]['NEXT'] = thread_id
                wings[1]['PREV'] = thread_id
                thread_id += 1

    # handling internal wings
    for i in range(1,len(wings)-1):

        if wings[i]['CHAIN'] != wings[i-1]['CHAIN']: # break between this and previous

            #if wings[i]['START'][2]:
            if True:#this allows zero threads
                threads.append(ExtThread(model,wings[i],'PREV',thread_id))
                wings[i]['PREV'] = thread_id
                thread_id += 1
        
        if wings[i]['CHAIN'] != wings[i+1]['CHAIN']: # break between this and next

            #if wings[i]['END'][2] != len(model.chains[wings[i]['END'][0][0]]['RES']) - 1:
            if True:#this allows zero threads
                threads.append(ExtThread(model,wings[i],'NEXT',thread_id))
                wings[i]['NEXT'] = thread_id
                thread_id += 1
    
        #elif wings[i]['END'][2] != wings[i+1]['START'][2] - 1: # elif thread between this and next wing has non-zero length
        else: #this allows zero threads
                threads.append(Thread(model,wings[i],wings[i+1],thread_id))
                wings[i]['NEXT'] = thread_id
                wings[i+1]['PREV'] = thread_id
                thread_id += 1 

    if wings:
    # handling last wing
        
        if wings[-1]['CHAIN'] != wings[-2]['CHAIN']: # if there is break between last and penult wing
            
            #if wings[-1]['START'][2]: # if start of last wing != start of chain
            if True:#this allows zero threads
                threads.append(ExtThread(model,wings[-1],'PREV',thread_id))
                wings[-1]['PREV'] = thread_id
                thread_id += 1
        
        #if wings[-1]['END'][2] != len(model.chains[wings[-1]['END'][0][0]]['RES']) - 1: # if end of last wing != end of chain
        if True:#this allows zero threads
            threads.append(ExtThread(model,wings[-1],'NEXT',thread_id))
            wings[-1]['NEXT'] = thread_id
            thread_id += 1

    #Adding full threads (threads that == full chain)

    na_chains = []
    not_empty = []
    empty     = []

    for ch in model.chains:

        if model.chains[ch]['TYPE'] in ('RNA','DNA'): na_chains.append(ch)

    for stem in model.fullstems:

        if stem['CHAIN1'] not in not_empty: not_empty.append(stem['CHAIN1'])
        if stem['CHAIN2'] not in not_empty: not_empty.append(stem['CHAIN2'])

    for ch in na_chains:

        if ch not in not_empty: empty.append(ch)

    for ch in empty:

        if not model.chains[ch]['GARBAGE']:
            
            nucl1 = model.chains[ch]['RES'][0]
            nucl2 = model.chains[ch]['RES'][-1]

            jmol = MakeJmol([x['DSSR'] for x in model.chains[ch]['RES']])

            threads.append({'ID'   : thread_id,
                            'MODEL': 1,
                            'CHAIN': ch,
                            'PREV' : None,
                            'NEXT' : None,
                            'START': [nucl1['DSSR'], 'RES', 0],
                            'END'  : [nucl2['DSSR'], 'RES', model.chains[ch]['LENGTH'] - 1],
                            'LEN'  : model.chains[ch]['LENGTH'],
                            'LINKS': 0,
                            'SEQ'  : ','.join(model.chains[ch]['SEQ2']),
                            'EXT'  : 1,
                            'FULL' : True,
                            'JMOL' : jmol}) #thread == chain
            thread_id += 1
   
    return threads

def transform(model): #clean bad threads (example: MG,MG,MG,MG,MG,MG from ligands)

    for th in model.threads:

        ch = th['CHAIN']

        if th['EXT']:

            if model.chains[ch]['GARBAGE'] or\
               model.chains[ch]['TYPE'] not in ('RNA','DNA') or\
               th['START'][1] == 'LIGANDS':

                if th['PREV']:

                    th['START'] = model.wings['FULL'][th['PREV']-1]['END']
                    th['END']   = th['START']

                elif th['NEXT']:

                    th['START'] = model.wings['FULL'][th['NEXT']-1]['START']
                    th['END']   = th['START']

                th['SEQ']   = ''
                th['LEN']   = 0

def relation(model):

    # adding THREAD to model.chains[ch][RES/LIGANDS][i] (residue['THREAD'])
    for thread in model.threads:

        if not thread['LEN']: continue      #it's important! else several nts have both wing and thread

        ch = thread['CHAIN']
        where = thread['START'][1]
        start = thread['START'][2]
        end   = thread['END'][2]

        for i in range(start,end+1): model.chains[ch][where][i]['THREAD'] = thread['ID']

    '''
    ###################################### this is code for testing
    print(model.chain_concat)
    print(model.chain_order)
    print('STEM','\t','ID','\t','START','\t','END','\t','PREV','\t',
              'NEXT','\t','PREVW','\t','NEXTW','\t','SEQ','\t','LEN')
    for w in model.wings['LU']:
        print(w['STEM'],'\t',w['ID'],'\t',w['TYPE'],'\t',w['START'],'\t',w['END'],'\t',w['PREV'],'\t',
              w['NEXT'],'\t',w['PREVW'],'\t',w['NEXTW'],'\t',w['SEQ'],'\t',w['LEN'])
    print('___________________________________________________________________________________')
    for w in model.threads:
        print(' \t',w['ID'],'\t',w['START'],'\t',w['END'],'\t',w['PREV'],'\t',
              w['NEXT'],' \t',' \t',w['SEQ'],'\t',w['LEN'],'\t',w['FULL'],w['EXT'])
    for bp in model.bpairs: print(bp['NUCL1'],bp['NUCL2'])
    ######################################
    #'''

def add(model):

    threads = Threads(model)

    model.threads = threads

    #transform(model)
    relation(model)
