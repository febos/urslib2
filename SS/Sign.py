'''
Created on 29.04.2014

@author: baulin
'''

''' this script emphasizes Elementary Closed Fragments (ECF), shows Signatures and Bracket Diagrams'''

def Diagram(model,cc,ID,ecfID):

    seq = sorted(cc,key=lambda x: model.chain_order[x])

    brack,slbrack,nuclseq = [],[],[]
    lens = []

    for ch in seq:
        brack   += [nucl['BRACKETS']   for nucl in model.chains[ch]['RES']]
        slbrack += [nucl['SLBRACKETS'] for nucl in model.chains[ch]['RES']]
        nuclseq += model.chains[ch]['SEQ2']
        lens.append(len(model.chains[ch]['RES']))
        model.chains[ch]['DIAGRAM'] = ID # chains.diagram

    brackets   = Brackets(brack)
    slbrackets = Brackets(slbrack)
    sllptrank  = lptRank(slbrack)
    sldibrank  = int(len(''.join(set(slbrackets)).replace('.','').replace('-',''))/2)

    seqbrack = ''
    for i in range(len(nuclseq)):
        seqbrack += nuclseq[i]+slbrackets[i]

    for i in range(len(seq)): # chains.brackets, chains.slbrackets

        model.chains[seq[i]]['BRACKETS']   =   brackets[sum(lens[:i]):sum(lens[:i])+lens[i]]
        model.chains[seq[i]]['SLBRACKETS'] = slbrackets[sum(lens[:i]):sum(lens[:i])+lens[i]]

    scheme = []
    intra_scheme = []
    flag = {'L':1,'R':-1}

    for w in model.wings['LU']:
        another = model.wings['LU'][w['ANOTHER']-1]
        if       w['CHAIN'] in seq and       w['START'][1] == 'RES' and \
           another['CHAIN'] in seq and another['START'][1] == 'RES':

            scheme.append(flag[w['TYPE']]*w['STEM'])

            if w['CHAIN'] == another['CHAIN']:
                intra_scheme.append(flag[w['TYPE']]*w['STEM'])

    # Renumeration
    id_num = {}
    num_id = {}
    num    = 1

    intra_id_num = {}
    intra_num_id = {}
    intra_num    = 1
    
    for i in range(len(scheme)):

        if scheme[i] > 0:
            id_num[scheme[i]] = num
            num_id[num]       = scheme[i]
            model.stems[scheme[i]-1]['DIAGRAM']   = ID  # stems.diagram
            model.stems[scheme[i]-1]['NUMINDIAG'] = num # stems.numindiag
            scheme[i] = num
            num += 1
        else: scheme[i] = -id_num[-scheme[i]]

    for i in range(len(intra_scheme)):

        if intra_scheme[i] > 0:
            intra_id_num[intra_scheme[i]] = intra_num
            intra_num_id[intra_num]       = intra_scheme[i]
            intra_scheme[i] = intra_num
            intra_num += 1
        else: intra_scheme[i] = -intra_id_num[-intra_scheme[i]]

    stembrack = Brackets(scheme)
    dibrank   = int(len(set(stembrack))/2)
    ecfs     = []
    ecf_inds = Ecfs(scheme)
    depths   = Depth(ecf_inds)
    parents  = Parents(ecf_inds,ecfID)

    intra_ecfID     = ecfID
    intra_stembrack = Brackets(intra_scheme)
    intra_dibrank   = int(len(set(intra_stembrack))/2)
    intra_ecf_inds  = Ecfs(intra_scheme)
    intra_depths    = Depth(intra_ecf_inds)
    intra_parents   = Parents(intra_ecf_inds,intra_ecfID)
    
    for ecf in ecf_inds:
        if ecf not in parents: parent = '\\N'
        else: parent = parents[ecf]
        ecfs.append(Ecf(model,num_id,ecf,ecf_inds,scheme,depths[ecf],parent,ID,ecfID))
        
        ecfID += 1

    for intra_ecf in intra_ecf_inds:
        if intra_ecf not in intra_parents: intra_parent = '\\N'
        else: intra_parent = intra_parents[intra_ecf]
        temp_intra_ecf = Ecf(model,intra_num_id,intra_ecf,
                             intra_ecf_inds,intra_scheme,intra_depths[intra_ecf],
                             intra_parent,ID,intra_ecfID,local=True)

        intra_ecfID += 1

        if temp_intra_ecf['SIGNATURE'] != 'aA': #adding intrachain sub-pseudoknots of interchain pseudoknots

            for temp_ecf in ecfs:

                if temp_intra_ecf['STEMS'] & temp_ecf['STEMS'] and temp_intra_ecf['STEMS'] != temp_ecf['STEMS']:

                    temp_intra_ecf['ID'] = ecfID
                    ecfs.append(temp_intra_ecf)
                    ecfID += 1

                    break
        
    depth   = max(list(depths.values())+[0,])
    lptrank = max([i['LPTRANK'] for i in ecfs]+[1,])

    return {'ID'        :            ID,
            'NODEL'     :             1,
            'ECFS'      :          ecfs,
            'SEQ'       : ','.join(seq),
            'NUCLSEQ'   :       nuclseq,
            'BRACKETS'  :      brackets,
            'SLBRACKETS':    slbrackets,
            'SLDIBRANK' :     sldibrank,
            'SLLPTRANK' :     sllptrank,
            'SCHEME'    :        scheme,
            'STEMBRACK' :     stembrack,
            'DIBRANK'   :       dibrank,
            'LPTRANK'   :       lptrank,
            'DEPTH'     :         depth,
            'SEQBRACK'  :      seqbrack}

def Ecf(model,num_id,ecf,ecf_inds,scheme,depth,parent,ID,ecfID,local=False):

    fullscheme = scheme[ecf[0]:ecf[1]+1]
    scheme     = offChildren(scheme,ecf_inds,ecf)

    ecfstems = set()

    # Renumeration
    num_num2 = {}
    num2_num = {}
    num2     = 1
    chainseq = []

    wingseq = []

    for i in range(len(scheme)):

        if scheme[i] > 0:
            num_num2[scheme[i]] = num2
            num2_num[num2]      = scheme[i]
            stem_ind = num_id[scheme[i]]-1

            if not local:
                model.stems[stem_ind]['ECF']      = ecfID    # stems.ecf
                model.stems[stem_ind]['NUMINECF'] = num2     # stems.numinecf
                model.wings['LU'][model.stems[stem_ind]['LEFT']-1]['ECF']  = ecfID # Lwing.ecf
                model.wings['LU'][model.stems[stem_ind]['RIGHT']-1]['ECF'] = ecfID # Rwing.ecf
                model.wings['LU'][model.stems[stem_ind]['LEFT']-1]['NUMINECF']  = num2  #Lwing.numinecf
                model.wings['LU'][model.stems[stem_ind]['RIGHT']-1]['NUMINECF'] = -num2 #Rwing.numinecf

            chainseq.append(model.wings['LU'][model.stems[stem_ind]['LEFT']-1]['CHAIN']) # for chainseq
            chainseq.append(model.wings['LU'][model.stems[stem_ind]['RIGHT']-1]['CHAIN'])# for chainseq
            
            scheme[i] = num2
            wingseq.append(str(num2)+':'+str(model.wings['LU'][model.stems[stem_ind]['LEFT']-1]['LEN'])+\
                           ':'+model.wings['LU'][model.stems[stem_ind]['LEFT']-1]['SEQ'])
            num2 += 1

            ecfstems.add(model.stems[stem_ind]['ID'])
        else:
            stem_ind = num_id[-scheme[i]]-1
            scheme[i] = -num_num2[-scheme[i]]
            wingseq.append(str(scheme[i])+':'+str(model.wings['LU'][model.stems[stem_ind]['RIGHT']-1]['LEN'])+\
                           ':'+model.wings['LU'][model.stems[stem_ind]['RIGHT']-1]['SEQ'])

    wingseq = '.'.join(wingseq)

    chainseq = ','.join(sorted(list(set(chainseq))))

    signature, intsign = Sign(scheme)
    brackets = Brackets(intsign)
    dibrank  = int(len(set(brackets))/2)
    lptrank  = lptRank(intsign)

    jmolbrackets = Brackets(scheme)
    jmol = '.'.join([jmolbrackets[i]+model.stems[num_id[num2_num[scheme[i]]]-1]['LOOPJMOL'] for i in range(len(jmolbrackets)) if scheme[i]>0])

    return {'ID'        :      ecfID,
            'MODEL'     :          1,
            'DIAGRAM'   :         ID,
            'FULLSCHEME': fullscheme,
            'SCHEME'    :     scheme,
            'SIGNATURE' :  signature,
            'BRACKETS'  :   brackets,
            'DIBRANK'   :    dibrank,
            'LPTRANK'   :    lptrank,
            'DEPTH'     :      depth,
            'PARENT'    :     parent,
            'WINGSEQ'   :    wingseq,
            'JMOL'      :       jmol,
            'CHAINSEQ'  :   chainseq,
            'STEMS'     :   ecfstems,
            'INTRA'     : int(local)}

def Parents(ecfs,shift): 

    parents = {} # key = ecf_inds; value = ecfID

    for i in range(len(ecfs)):
        for j in range(i-1,-1,-1):
            if  ecfs[j][0] < ecfs[i][0] and ecfs[i][1] < ecfs[j][1]:
                parents[ecfs[i]] = shift+j
                break
    return parents

def BFS(adj, cc, v):

    seen = []
    pop1 = []
    pop2 = []

    for v2 in cc:   # filling first version of pop1        
        if v2 not in seen:
            pop1.append(v2)
            seen.append(v2)
            
    while pop1:
        for v2 in pop1:
            if v2 == v: return True
            for new_v in adj[v2]:
                if new_v not in seen:
                    pop2.append(new_v)
                    seen.append(new_v)
        pop1 = pop2
        pop2 = []
        
    return False

def CC(adj):

    ccs = []

    for v in adj:

        found = 0

        for i in range(len(ccs)):

            if BFS(adj,ccs[i],v):

                ccs[i].append(v)
                found = 1

        if not found: ccs.append([v,])

    return ccs

def list_to_dict(neg_word): # key=value; value=index

    dorw = {}
    for i in range(0,len(neg_word)):
        if neg_word[i] in ('.','-'): continue
        dorw[neg_word[i]] = i

    return dorw

def Depth(ecfs):

    depth = {}

    for ecf in ecfs: depth[ecf] = 0

    for i in range(0,len(ecfs)):
        for j in range(0,len(ecfs)):
            if ecfs[j][0] < ecfs[i][0] and ecfs[j][1] > ecfs[i][1]:
                depth[ecfs[i]] += 1
    return depth

def lptRank(intsign):
    
    rank = 1
    ngis = list_to_dict(intsign)
    
    for s in [i for i in intsign if i not in ('.','-') and i > 0]:
        end = ngis[-s]
        k   = 1
        for i in range(ngis[s]+1, ngis[-s]):
            if intsign[i] in ('.','-'): continue
            if intsign[i] > 0 and ngis[-intsign[i]] > end: # if intersects 
                end = ngis[-intsign[i]] # for checking intersecting with all previous
                k += 1
        if k > rank: rank = k

    return rank

def offChildren(scheme,ecfs,ecf):

    scheme2 = []
    masc    = [1]*len(scheme)

    for m in ecfs:
        if m[0] > ecf[0]:
            for i in range(m[0],m[1]+1): masc[i] = 0

    for i in range(ecf[0],ecf[1]+1):
        if masc[i]: scheme2.append(scheme[i])

    return scheme2

def Sign(presign):
    
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    presign2 = [str(i) for i in presign]
    presign_str = ','.join(presign2)
    nigs = {}
    for i in range(len(presign)): nigs[presign[i]] = i

    change = {}

    i   = 0
    old = 0

    while i < len(presign) - 2:

        if presign[i] > 0 and presign[i+1] > 0 and\
           nigs[-presign[i]] - nigs[-presign[i+1]] == 1:
            i += 1
        
        elif i != old:

            flag = [',',',',',',','] # flags because '1,2':'1' makes '21,22'->'212'
            if old==0:                               flag[0]=''
            if i+1==len(presign2):                   flag[1]=''
            if nigs[-presign[i]]==0:                 flag[2]=''
            if nigs[-presign[old]]+1==len(presign2): flag[3]=''
            change[flag[0]+','.join(presign2[old:i+1])+flag[1]] = flag[0]+presign2[old]+flag[1]
            change[flag[2]+','.join(presign2[nigs[-presign[i]]:nigs[-presign[old]]+1])+flag[3]] = flag[2]+str(-presign[old])+flag[3]
            old = i+1
            i   = i+1
        else:
            i   += 1
            old += 1

    for ch in change:
        presign_str = presign_str.replace(ch,change[ch])
    
    presign = presign_str.split(',')
    alter = {}
    Set = []
    for i in presign:
        if abs(int(i)) not in Set: Set.append(abs(int(i)))

    for i in range(0,len(Set)):
        alter[Set[i]]  = alphabet[i]
        alter[-Set[i]] = alphabet[i].upper()

    sign = []

    for i in presign: sign.append(alter[int(i)])

    return ''.join(sign), [int(i) for i in presign]

def Ecfs(word):

    stack = []
    ecfs  = []
    dorw  = list_to_dict(word)

    for i in range(0,len(word)):

        if word[i] > 0:
            n = dorw[-word[i]]

            if not bool(stack) or stack[-1][1] > n: stack.append([i,n])
            else:
                for f in range(0,len(stack)):
                    if stack[f][1] < n:
                        stack[f][1] = n
                        while stack[f] != stack[-1]: stack.pop()
                        break

        elif i == stack[-1][1]:
            ecf = stack[-1]
            ecfs.append(tuple(ecf))
            stack.pop()
        else:
            while stack[-1][0] > dorw[-word[i]]: stack.pop()

    return sorted(ecfs,key = lambda x: x[0])

def Brackets(seq):
    
    diagram = seq[:]
    dic = {'(':1,'[':2,'{':3,'<':4,'!':5,'C':6,'6':7}
    left  = ('(','[','{','<','!','C','6')
    right = (')',']','}','>','?','D','9')

    qes = list_to_dict(seq)
    opened  = []
   
    for i in range(len(seq)):
        
        if seq[i] in ('.','-'): continue
        
        level = 0
        
        if seq[i] > 0:

            j = qes[-seq[i]]

            for num in opened:

                j0 = qes[-num]

                if j > j0 and diagram[j0] == right[level]: level += 1

            diagram[i] = left[level]
            diagram[j] = right[level]
            opened.append(seq[i])
            # sorting in order of increasing of brackets level
            opened = sorted(opened, key = lambda x: dic[diagram[qes[x]]])

        else: opened.remove(-seq[i])            
    

    return ''.join(diagram)

def Adj_plus_nucls(model): # Get Adjacency list of chains + mark nucls for bracket diagrams

    chains = {} # adjacency list of chains 
    chainpairs = []
    seen_nucls = {}

    for bp in model.bpairs:

        if bp['TYPE'] not in ('WC','WB'): continue

        if model.chains[bp['CHAIN1']]['TYPE'] in ('RNA','DNA') and bp['NUCL1'][1] == 'RES' and\
           model.chains[bp['CHAIN2']]['TYPE'] in ('RNA','DNA') and bp['NUCL2'][1] == 'RES':

            # if two WC/WB bps conflict
            if str(bp['NUCL1']) in seen_nucls or str(bp['NUCL2']) in seen_nucls: continue
            seen_nucls[str(bp['NUCL1'])] = 1
            seen_nucls[str(bp['NUCL2'])] = 1

            if bp['CHAIN1'] not in chains: chains[bp['CHAIN1']] = []

            if bp['CHAIN1'] != bp['CHAIN2']:
                
                if bp['CHAIN2'] not in chains: chains[bp['CHAIN2']] = []

                if bp['CHAIN1']+bp['CHAIN2'] not in chainpairs and\
                   bp['CHAIN2']+bp['CHAIN1'] not in chainpairs:

                    chainpairs.append(bp['CHAIN1']+bp['CHAIN2'])
                    chains[bp['CHAIN1']].append(bp['CHAIN2'])
                    chains[bp['CHAIN2']].append(bp['CHAIN1'])

            nucl1,nucl2 = bp['NUCL1'],bp['NUCL2']

            if model.chain_order[bp['CHAIN1']]*10000 + nucl1[2] < \
               model.chain_order[bp['CHAIN2']]*10000 + nucl2[2]:
                a,b = 1,-1
            else: a,b = -1,1
            
            model.chains[bp['CHAIN1']]['RES'][nucl1[2]]['SLBRACKETS'] = a*bp['ID']
            model.chains[bp['CHAIN2']]['RES'][nucl2[2]]['SLBRACKETS'] = b*bp['ID']

            if bp['STEM']:
                model.chains[bp['CHAIN1']]['RES'][nucl1[2]]['BRACKETS'] = a*bp['ID']
                model.chains[bp['CHAIN2']]['RES'][nucl2[2]]['BRACKETS'] = b*bp['ID']

    return chains              

def Main(model):

    ID = 1
    ecfID = 1
    diagrams = []
    
    ch_adj = Adj_plus_nucls(model)
    ccs    = CC(ch_adj)

    for cc in sorted(ccs,key = lambda x:x[0]):

        diagrams.append(Diagram(model,cc,ID,ecfID))
        ecfID += len(diagrams[-1]['ECFS'])
        ID += 1      

    return diagrams
    
def add(model):

    diagrams = Main(model)
    
    model.diagrams = diagrams


