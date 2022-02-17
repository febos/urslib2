
def add(model):

    chain_order  = {}   # key - chain letter; value - int
    chain_concat = []   # pairs of chain-letters that have shared bps

    for bp in model.bpairs:

        if bp['CHAIN1'] != bp['CHAIN2'] and\
           bp['CHAIN1'] + '-' +  bp['CHAIN2'] not in chain_concat and\
           bp['CHAIN2'] + '-' +  bp['CHAIN1'] not in chain_concat:

            chain_concat.append(bp['CHAIN1'] + '-' + bp['CHAIN2'])

        if bp['CHAIN1'] not in chain_order: chain_order[bp['CHAIN1']] = 0
        if bp['CHAIN2'] not in chain_order: chain_order[bp['CHAIN2']] = 0

    if model.bpairs: 
        # if bpairs is not empty but there are no RES-RES pairs
        if model.bpairs[0]['CHAIN1'] != model.bpairs[0]['CHAIN2']:

            chain_order[model.bpairs[0]['CHAIN2']] += 1 

    ####### preliminary sorting (maybe it is not necessary)
    alph = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz '
    for ch in chain_order:
        if len(ch)==1: chain_order[ch] += alph.find(ch)
        else:          chain_order[ch] += alph.find(ch[0]) + alph.find(ch[1])/10
    #######################################################
    
    for i in range(len(model.bpairs)):

        if i:

            ch1 = model.bpairs[i-1]['CHAIN1']
            where1 = model.bpairs[i-1]['NUCL1'][1] #RES/LIGANDS
        
        ch3 = model.bpairs[i]['CHAIN1']
        ch4 = model.bpairs[i]['CHAIN2']
        where3 = model.bpairs[i]['NUCL1'][1]
        where4 = model.bpairs[i]['NUCL2'][1]

        if i and ch1 != ch3 and where1 == where3 == 'RES': chain_order[ch3] = chain_order[ch1] + 10000
        if       ch3 != ch4 and where3 == where4 == 'RES': chain_order[ch4] = chain_order[ch3] + 1
            
    for ch in chain_order:
        if len(ch)==1: chain_order[ch] = float(chain_order[ch]) + alph.find(ch)/100
        else:          chain_order[ch] = float(chain_order[ch]) + alph.find(ch[0])/100 + alph.find(ch[1])/1000

    model.chain_order  = chain_order
    model.chain_concat = chain_concat
