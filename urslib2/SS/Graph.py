
def BFS(edges, adj_list, cc, new_e, Type):

    seen = []
    pop1 = []
    pop2 = []

    if   Type == 'NUCL' : vert1, vert2 = new_e['NUCL1'][0], new_e['NUCL2'][0]
    elif Type == 'FSTEM': vert1, vert2 = new_e['FSTEM1'],   new_e['FSTEM2']

    for old in cc:   # filling first version of pop1

        if   Type == 'NUCL' : vert3, vert4 = edges[old-1]['NUCL1'][0], edges[old-1]['NUCL2'][0]
        elif Type == 'FSTEM': vert3, vert4 = edges[old-1]['FSTEM1'],   edges[old-1]['FSTEM2']

        for vert in (vert3, vert4):

            if vert not in seen:

                pop1.append(vert)
                seen.append(vert)
    while pop1:

        for vert in pop1:

            if vert == vert1: return True
            if vert == vert2: return True

            for new_vert in adj_list[vert]:

                if new_vert not in seen:

                    pop2.append(new_vert)
                    seen.append(new_vert)
        pop1 = pop2
        pop2 = []
        
    return False

def Adjacency_List(edges, Type):

    adj_list = {}

    for e in edges:

        if   Type == 'NUCL'  : vert1, vert2 = e['NUCL1'][0], e['NUCL2'][0]
        elif Type == 'FSTEM' : vert1, vert2 = e['FSTEM1'],   e['FSTEM2']

        if vert1 == vert2: continue

        if vert1 not in adj_list: adj_list[vert1] = [vert2,]
        else:                     adj_list[vert1].append(vert2)

        if vert2 not in adj_list: adj_list[vert2] = [vert1,]
        else:                     adj_list[vert2].append(vert1)

    return adj_list

def ConnectedComponents(edges, Type): 

    ccs = []
    
    adj_list = Adjacency_List(edges, Type)

    for e in edges:

        found = 0

        for i in range(len(ccs)):

            if BFS(edges,adj_list,ccs[i],e,Type):

                ccs[i].append(e['ID'])
                found = 1

        if not found: ccs.append([e['ID'],])

    return ccs

def normalized(Dict):

    rename = {}

    new_dict = {}

    verts = sorted(list(Dict.keys()))

    for i in range(len(verts)): rename[verts[i]] = i

    for i in range(len(verts)):
        new_dict[i] = Dict[verts[i]]
        for j in range(len(new_dict[i])):
            new_dict[i][j] = rename[new_dict[i][j]]
        new_dict[i].sort()

    return new_dict

def Build_Type(cc_edges, Type):

    adj_list = Adjacency_List(cc_edges,Type)

    return invariant(adj_list)

def invariant(adj):

    
    result = []
    for i in adj: result.append(invariant_str(adj,i))
    return ','.join(sorted(result))

def invariant_str(adj, x):

    seen  = []
    global res
    res   = str(len(adj[x]))+'('
    if res[0] == '0': res='0'

    def search(seen,x):
        global res
        for i in adj[x]:
            if i not in seen:
                br = ''
                for ch in adj[i]:
                    if ch not in seen+adj[x]: br='('
                res += str(len(adj[i]))+br
                search(seen+adj[x],i)
                if br:res += ')'
                res += '.'

    def Sorting(res):
        br = res.find('(')
        if br == -1: return res
        res2 = res[:br+1]
        k = 0
        last = br+1
        pieces = []
        for i in range(len(res)):
            if i < br+1: continue
            if   res[i] == '(': k += 1
            elif res[i] == ')': k -= 1
            elif res[i]=='.' and not k:
                pieces.append(res[last:i])
                last = i+1
        #res2 += '.'.join([Sorting(i) for i in sorted(pieces)]) + ')' - uniqueness is not respected
        res2 += '.'.join(sorted([Sorting(i) for i in pieces])) + ')'

        return res2

    search([x],x)
    if res[0] != '0':res += ')'
    res = Sorting(res)
    
    return res
