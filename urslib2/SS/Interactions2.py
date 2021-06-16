
    
def Process(model, type1, type2, restype):

    needed_atoms = [[],[]] 

    xmin, ymin, zmin =  1000000,  1000000,  1000000
    xmax, ymax, zmax = -1000000, -1000000, -1000000

    for ch in model.chains:

        for pl in ('RES','LIGANDS'):
            
            for i in range(len(model.chains[ch][pl])):

                if restype[model.chains[ch][pl][i]['NAME']][0] in type1:

                    for j in range(len(model.chains[ch][pl][i]['ATOMS'])):

                        atom = model.chains[ch][pl][i]['ATOMS'][j]

                        if atom['X'] > xmax: xmax = atom['X']
                        if atom['X'] < xmin: xmin = atom['X']
                        if atom['Y'] > ymax: ymax = atom['Y']
                        if atom['Y'] < ymin: ymin = atom['Y']
                        if atom['Z'] > zmax: zmax = atom['Z']
                        if atom['Z'] < zmin: zmin = atom['Z']

                        needed_atoms[0].append([ch,pl,i,j])

                if restype[model.chains[ch][pl][i]['NAME']][0] in type2:

                    for j in range(len(model.chains[ch][pl][i]['ATOMS'])):

                        atom = model.chains[ch][pl][i]['ATOMS'][j]

                        if atom['X'] > xmax: xmax = atom['X']
                        if atom['X'] < xmin: xmin = atom['X']
                        if atom['Y'] > ymax: ymax = atom['Y']
                        if atom['Y'] < ymin: ymin = atom['Y']
                        if atom['Z'] > zmax: zmax = atom['Z']
                        if atom['Z'] < zmin: zmin = atom['Z']

                        needed_atoms[1].append([ch,pl,i,j])

    return needed_atoms, [xmin, xmax, ymin, ymax, zmin, zmax]



def Cubics(model, needed_atoms, minmax, side):

    cubics = [{},{}]
    
    xmin, xmax, ymin, ymax, zmin, zmax = minmax

    length = xmax - xmin
    width  = ymax - ymin
    height = zmax - zmin

    n, m, k = int(length//side), int(width//side), int(height//side)

    bitlen = length/n
    bitwid = width/m
    bithei = height/k

    for i in range(2):

        for atom in needed_atoms[i]:

            at = model.chains[atom[0]][atom[1]][atom[2]]['ATOMS'][atom[3]]

            xnum = int((at['X']-xmin)/bitlen-0.001)+1
            ynum = int((at['Y']-ymin)/bitwid-0.001)+1
            znum = int((at['Z']-zmin)/bithei-0.001)+1

            if (xnum,ynum,znum) not in cubics[i]: cubics[i][(xnum,ynum,znum)] = []

            cubics[i][(xnum,ynum,znum)].append(atom)

    return cubics





def Search(model, cubics, MaxDist):

    def AdjCubes(a,b,c):

        return [(a-1,b-1,c-1),(a-1,b-1,c),(a-1,b-1,c+1),(a-1,b,c-1),(a-1,b,c),(a-1,b,c+1),(a-1,b+1,c-1),(a-1,b+1,c),(a-1,b+1,c+1),
                (a,  b-1,c-1),(a,  b-1,c),(a,  b-1,c+1),(a,  b,c-1),(a,  b,c),(a,  b,c+1),(a,  b+1,c-1),(a,  b+1,c),(a,  b+1,c+1),
                (a+1,b-1,c-1),(a+1,b-1,c),(a+1,b-1,c+1),(a+1,b,c-1),(a+1,b,c),(a+1,b,c+1),(a+1,b+1,c-1),(a+1,b+1,c),(a+1,b+1,c+1)]

    def DistSquare(atom,atom2):

        at = model.chains[atom[0]][atom[1]][atom[2]]['ATOMS'][atom[3]]
        at2 = model.chains[atom2[0]][atom2[1]][atom2[2]]['ATOMS'][atom2[3]]

        return (at['X']-at2['X'])**2 + (at['Y']-at2['Y'])**2 + (at['Z']-at2['Z'])**2

    potentials = []

    for (a,b,c) in cubics[0]:

        for atom in cubics[0][(a,b,c)]:

            for (d,e,f) in AdjCubes(a,b,c):

                if (d,e,f) in cubics[1]:

                    for atom2 in cubics[1][(d,e,f)]:

                        if atom[:3] != atom2[:3] and not (atom[:2] == atom2[:2] and atom[2] > atom2[2]):
                            
                            distsquare = DistSquare(atom,atom2) 

                            if distsquare <= MaxDist**2:

                                dist = pow(distsquare,0.5)

                                potentials.append([atom,atom2,dist])

    return potentials



def Atompairs(model, type1 = '', type2 = '', dist = 4, restype = {}):

    side = dist + 1

    needed_atoms, minmax = Process(model, type1, type2, restype)
    
    cubics = Cubics(model, needed_atoms, minmax, side)

    del needed_atoms, minmax

    pairs = Search(model, cubics, dist)

    del cubics

    atompairs = []

    for pair in pairs:

        ch1,pl1,i1,j1 = pair[0]
        ch2,pl2,i2,j2 = pair[1]

        atompair = {'DSSR1': model.chains[ch1][pl1][i1]['DSSR'],
                    'DSSR2': model.chains[ch2][pl2][i2]['DSSR'],
                    'atom1': model.chains[ch1][pl1][i1]['ATOMS'][j1]['NAME'],
                    'atom2': model.chains[ch2][pl2][i2]['ATOMS'][j2]['NAME'],
                    'dist' : pair[2],
                    'type' : restype[model.chains[ch1][pl1][i1]['NAME']][0]+\
                             restype[model.chains[ch2][pl2][i2]['NAME']][0]}

        if not (atompair['type'][0]==atompair['type'][1] and i1 > i2):

            atompairs.append(atompair)

    del pairs

    return atompairs







    
