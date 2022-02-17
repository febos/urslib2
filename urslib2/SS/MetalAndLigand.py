import math

try:    import urslib2.SS.Maps as Maps
except:
    try:    import SS.Maps as Maps
    except: import Maps

def Atompair(model,atom1,pair,bypass):

    ch1,res1,at1,where1 = atom1[0],  atom1[1],  atom1[2], atom1[3]
    ch2,res2,at2,where2 = pair[0][0],pair[0][1],pair[0][2], pair[0][3]

    ratom = model.chains[ch1][where1][res1]['ATOMS'][at1]
    patom = model.chains[ch2][where2][res2]['ATOMS'][at2]

    rna  = []
    prot = []

    # chain ID
    rna.append(model.chains[ch1]['ID'])
    prot.append(model.chains[ch2]['ID'])
    # residue name
    rna.append(ratom['RESNAME'])
    prot.append(patom['RESNAME'])
    # residue dssr code
    rna.append(model.chains[ch1][where1][res1]['DSSR'])
    prot.append(model.chains[ch2][where2][res2]['DSSR'])
    # atom ID
    rna.append(ratom['ID'])
    prot.append(patom['ID'])
    # atom name
    rna.append(ratom['NAME'])
    prot.append(patom['NAME'])
    # atom element
    rna.append(ratom['ELEM'])
    prot.append(patom['ELEM'])
    # atom jmol code
    rna.append(ratom['RESNUM'][:-1]+':'+ch1+'.'+ratom['NAME'])
    rightchain = model.chains[ch2][where2][res2]['DSSR'][:model.chains[ch2][where2][res2]['DSSR'].find('.')]
    prot.append(patom['RESNUM'][:-1]+':'+rightchain+'.'+patom['NAME']+',HOH')

    dist  = round(pair[1],2)

    if bypass:
        rna,prot      = prot,rna
        atom1,pair[0] = pair[0],atom1

    return {'ID'      :               0,
            'MODEL'   :               1,
            'RCHAIN'  :          rna[0],
            'LCHAIN'  :         prot[0],
            'NUCL'    :          rna[1],
            'LIGAND'  :         prot[1],
            'RDSSR'   :          rna[2],
            'LDSSR'   :         prot[2],
            'RATOMID' :          rna[3],
            'LATOMID' :         prot[3],
            'RNAME'   :          rna[4],
            'LNAME'   :         prot[4],
            'RELEM'   :          rna[5],
            'LELEM'   :         prot[5],
            'RJMOL'   :          rna[6],
            'LJMOL'   :         prot[6],
            'DIST'    :            dist,
            'MONOPAIR':               0,
            'ADRESS'  : [atom1,pair[0]]}

def Monopair(atompairs,mon,ID):

    firstpair = atompairs[mon[0]-1]

    rjmol = firstpair['RJMOL'][:firstpair['RJMOL'].find('.')]
    pjmol = firstpair['LJMOL'][:firstpair['LJMOL'].find('.')]

    return {'ID'        :                  ID,
            'MODEL'     :                   1,
            'RCHAIN'    : firstpair['RCHAIN'],
            'LCHAIN'    : firstpair['LCHAIN'],
            'NUCL'      :   firstpair['NUCL'],
            'LIGAND'    : firstpair['LIGAND'],
            'RDSSR'     :  firstpair['RDSSR'],
            'LDSSR'     :  firstpair['LDSSR'],
            'RJMOL'     :               rjmol,
            'LJMOL'     :               pjmol,
            'APIDS'     :                 mon,
            'ATOMPAIRS' :            len(mon)}
    
def Process(model,Class):

    needed_atoms = [] # [RNA/Metal or RNA/Unknown, donor/acceptor, ch, res_index, atom_index, RES/LIGANDS]

    xmin, ymin, zmin =  1000000,  1000000,  1000000
    xmax, ymax, zmax = -1000000, -1000000, -1000000

    for ch in model.chains:

        for where in ('RES','LIGANDS'):

            for i in range(len(model.chains[ch][where])):

                if model.chains[ch][where][i]['TYPE'] in ('RNA',Class): 

                    res  = model.chains[ch][where][i]['NAME']
                    Type = model.chains[ch][where][i]['TYPE']

                    for j in range(len(model.chains[ch][where][i]['ATOMS'])):

                        atom = model.chains[ch][where][i]['ATOMS'][j]

                        if atom['X'] > xmax: xmax = atom['X']
                        if atom['X'] < xmin: xmin = atom['X']
                        if atom['Y'] > ymax: ymax = atom['Y']
                        if atom['Y'] < ymin: ymin = atom['Y']
                        if atom['Z'] > zmax: zmax = atom['Z']
                        if atom['Z'] < zmin: zmin = atom['Z']
                        #if model.chains[ch][where][i]['DSSR'] in ('0.C.2533.','0.MG.2925.'): print([Type,'',ch,i,j,where])
                        needed_atoms.append([Type,'',ch,i,j,where])

    return needed_atoms, [xmin, xmax, ymin, ymax, zmin, zmax]

def Cubics(model,needed_atoms,minmax,side):

    cubics = [{},{}] # [rna, metal/unknown]

    xmin, xmax, ymin, ymax, zmin, zmax = minmax[0], minmax[1], minmax[2], minmax[3], minmax[4], minmax[5]

    length = xmax - xmin
    width  = ymax - ymin
    height = zmax - zmin

    n, m, k = int(length//side), int(width//side), int(height//side)

    if not n: n = 1
    if not m: m = 1
    if not k: k = 1

    bitlen = length/n
    bitwid = width/m
    bithei = height/k

    for atom in needed_atoms:

        if atom[0] == 'RNA': Type = 0
        else:                Type = 1

        at = model.chains[atom[2]][atom[5]][atom[3]]['ATOMS'][atom[4]]

        if not bitlen: xnum = 1
        else:          xnum = int((at['X']-xmin)/bitlen-0.001)+1
        if not bitwid: ynum = 1
        else:          ynum = int((at['Y']-ymin)/bitwid-0.001)+1
        if not bithei: znum = 1
        else:          znum = int((at['Z']-zmin)/bithei-0.001)+1

        if (xnum,ynum,znum) not in cubics[Type]: cubics[Type][(xnum,ynum,znum)] = []
        #if [atom[2],atom[3],atom[4],atom[5]] == ['0', 2531, 0, 'RES']: print(xnum,ynum,znum)
        cubics[Type][(xnum,ynum,znum)].append([atom[2],atom[3],atom[4],atom[5]])

    return cubics

def Search(model,cubics,MaxDist):

    def AdjCubes(a,b,c):

        return [(a-1,b-1,c-1),(a-1,b-1,c),(a-1,b-1,c+1),(a-1,b,c-1),(a-1,b,c),(a-1,b,c+1),(a-1,b+1,c-1),(a-1,b+1,c),(a-1,b+1,c+1),
                (a,  b-1,c-1),(a,  b-1,c),(a,  b-1,c+1),(a,  b,c-1),(a,  b,c),(a,  b,c+1),(a,  b+1,c-1),(a,  b+1,c),(a,  b+1,c+1),
                (a+1,b-1,c-1),(a+1,b-1,c),(a+1,b-1,c+1),(a+1,b,c-1),(a+1,b,c),(a+1,b,c+1),(a+1,b+1,c-1),(a+1,b+1,c),(a+1,b+1,c+1)]

    def DistSquare(atom,atom2):

        at = model.chains[atom[0]][atom[3]][atom[1]]['ATOMS'][atom[2]]
        at2 = model.chains[atom2[0]][atom2[3]][atom2[1]]['ATOMS'][atom2[2]]

        return (at['X']-at2['X'])**2 + (at['Y']-at2['Y'])**2 + (at['Z']-at2['Z'])**2

    num = 1 # will go through Type with minimal number of atoms

    potentials = {}

    for (a,b,c) in cubics[num]:

        for atom in cubics[num][(a,b,c)]:

            for (d,e,f) in AdjCubes(a,b,c):

                if (d,e,f) in cubics[1-num]:

                    for atom2 in cubics[1-num][(d,e,f)]:
                            
                        distsquare = DistSquare(atom,atom2) 

                        if distsquare <= MaxDist**2:

                            dist = pow(distsquare,0.5)

                            if (atom[0],atom[1],atom[2],atom[3]) not in potentials: potentials[(atom[0],atom[1],atom[2],atom[3])] = []

                            potentials[(atom[0],atom[1],atom[2],atom[3])].append([atom2,dist])
    return potentials,num

def Atompairs(model,Type,dist):

    needed_atoms, minmax = Process(model,Type)
    
    cubics = Cubics(model,needed_atoms,minmax,dist)
    del needed_atoms,minmax

    potentials,bypass = Search(model,cubics,dist)
    del cubics
    #raw_atompairs = Clean(model,potentials)
    raw_atompairs = potentials 

    atompairs = []

    for atom1 in raw_atompairs:
        for pair in raw_atompairs[atom1]:
            atompairs.append(Atompair(model,atom1,pair,bypass))

    atompairs.sort(key= lambda x: (x['RJMOL'][x['RJMOL'].find(':')+1],
                                   int(x['RJMOL'][:x['RJMOL'].find(':')]),
                                   x['RJMOL'][x['RJMOL'].find('.'):]))
    ID = 1
    for i in range(len(atompairs)):
        atompairs[i]['ID'] = ID
        ID += 1

    return atompairs

def Monopairs(model,atompairs):

    monopairs = {} # (ch1,where1,res1,ch2,where2,res2): [atompair_ids]

    for ap in atompairs:

        ch1,    ch2    = ap['ADRESS'][0][0], ap['ADRESS'][1][0]
        where1, where2 = ap['ADRESS'][0][3], ap['ADRESS'][1][3]
        res1,   res2   = ap['ADRESS'][0][1], ap['ADRESS'][1][1]

        if (ch1,where1,res1,ch2,where2,res2) not in monopairs: monopairs[(ch1,where1,res1,ch2,where2,res2)] = []

        monopairs[(ch1,where1,res1,ch2,where2,res2)].append(ap['ID'])

    monopairs2 = []
    ID = 1

    for mon in sorted(monopairs.values()):
        monopairs2.append(Monopair(atompairs,mon,ID))
        ID += 1

    for mp in monopairs2:
        for i in mp['APIDS']:
            atompairs[i-1]['MONOPAIR'] = mp['ID']

    #for i in  atompairs: print(i['NUCL'],i['RJMOL'],i['DIST'],i['AMINO'],i['PJMOL'],i['ID'],i['MONOPAIR'],i['POWER'],i['RCOS'],i['PCOS'])
    #for i in monopairs2: print(i['NUCL'],i['RJMOL'],i['ATOMPAIRS'],i['AMINO'],i['PJMOL'],i['ID'])    

    return monopairs2,atompairs

def add(model):

    atompairsL, atompairsM = [], []
    monopairsL, monopairsM = [], []

    distance = 5

    if model.headers['LIGLIST']:
        atompairsL            = Atompairs(model,'Unknown',distance)
        monopairsL,atompairsL = Monopairs(model,atompairsL)

    if model.headers['METLIST']:
        atompairsM            = Atompairs(model,'Metal',distance)
        monopairsM,atompairsM = Monopairs(model,atompairsM)

    model.atompairsL = atompairsL
    model.monopairsL = monopairsL
    model.atompairsM = atompairsM
    model.monopairsM = monopairsM

    
