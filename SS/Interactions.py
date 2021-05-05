import math

try:    import urslib2.SS.Maps as Maps
except:
    try:    import SS.Maps as Maps
    except: import Maps

def Atompair(model,atom1,pair,bypass):

    ch1,res1,at1 = atom1[0],  atom1[1],  atom1[2]
    ch2,res2,at2 = pair[0][0],pair[0][1],pair[0][2]

    ratom = model.chains[ch1]['RES'][res1]['ATOMS'][at1]
    patom = model.chains[ch2]['RES'][res2]['ATOMS'][at2]

    rna  = []
    prot = []

    # chain ID
    rna.append(model.chains[ch1]['ID'])
    prot.append(model.chains[ch2]['ID'])
    # residue name
    rna.append(ratom['RESNAME'])
    prot.append(patom['RESNAME'])
    # residue dssr code
    rna.append(model.chains[ch1]['RES'][res1]['DSSR'])
    prot.append(model.chains[ch2]['RES'][res2]['DSSR'])
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
    prot.append(patom['RESNUM'][:-1]+':'+ch2+'.'+patom['NAME'])

    power = round(pair[1],2)
    dist  = round(pair[2],2)
    cos1  = round(pair[3],2)
    cos2  = round(pair[4],2)

    if bypass:
        rna,prot      = prot,rna
        cos1,cos2     = cos2,cos1
        atom1,pair[0] = pair[0],atom1

    return {'ID'      :               0,
            'MODEL'   :               1,
            'RCHAIN'  :          rna[0],
            'PCHAIN'  :         prot[0],
            'NUCL'    :          rna[1],
            'AMINO'   :         prot[1],
            'RDSSR'   :          rna[2],
            'PDSSR'   :         prot[2],
            'RATOMID' :          rna[3],
            'PATOMID' :         prot[3],
            'RNAME'   :          rna[4],
            'PNAME'   :         prot[4],
            'RELEM'   :          rna[5],
            'PELEM'   :         prot[5],
            'RJMOL'   :          rna[6],
            'PJMOL'   :         prot[6],
            'POWER'   :           power,
            'DIST'    :            dist,
            'RCOS'    :            cos1,
            'PCOS'    :            cos2,
            'MONOPAIR':               0,
            'ADRESS'  : [atom1,pair[0]]}

def Monopair(atompairs,mon,ID):

    firstpair = atompairs[mon[0]-1]

    rjmol = firstpair['RJMOL'][:firstpair['RJMOL'].find('.')]
    pjmol = firstpair['PJMOL'][:firstpair['PJMOL'].find('.')]

    return {'ID'        :                  ID,
            'MODEL'     :                   1,
            'RCHAIN'    : firstpair['RCHAIN'],
            'PCHAIN'    : firstpair['PCHAIN'],
            'NUCL'      :   firstpair['NUCL'],
            'AMINO'     :  firstpair['AMINO'],
            'RDSSR'     :  firstpair['RDSSR'],
            'PDSSR'     :  firstpair['PDSSR'],
            'RJMOL'     :               rjmol,
            'PJMOL'     :               pjmol,
            'APIDS'     :                 mon,
            'ATOMPAIRS' :            len(mon)}
    
def Process(model,allin=False):

    needed_atoms = [] # [RNA/Protein, donor/acceptor, ch, res_index, atom_index]

    xmin, ymin, zmin =  1000000,  1000000,  1000000
    xmax, ymax, zmax = -1000000, -1000000, -1000000

    for ch in model.chains:

        if model.chains[ch]['TYPE'] in ('RNA','Protein'):

            for i in range(len(model.chains[ch]['RES'])):

                if model.chains[ch]['RES'][i]['TYPE'] == model.chains[ch]['TYPE'] and\
                   (model.chains[ch]['RES'][i]['NAME'] in Maps.donors_acceptors or allin):

                    model.chains[ch]['RES'][i]['ATOMNAMEDICT'] = {}

                    res  = model.chains[ch]['RES'][i]['NAME']
                    Type = model.chains[ch]['RES'][i]['TYPE']

                    for j in range(len(model.chains[ch]['RES'][i]['ATOMS'])):

                        atom = model.chains[ch]['RES'][i]['ATOMS'][j]

                        model.chains[ch]['RES'][i]['ATOMNAMEDICT'][atom['NAME']] = j

                        if allin or atom['NAME'] in Maps.neighbors[res]:

                            if atom['X'] > xmax: xmax = atom['X']
                            if atom['X'] < xmin: xmin = atom['X']
                            if atom['Y'] > ymax: ymax = atom['Y']
                            if atom['Y'] < ymin: ymin = atom['Y']
                            if atom['Z'] > zmax: zmax = atom['Z']
                            if atom['Z'] < zmin: zmin = atom['Z']

                            needed_atoms.append([Type,'',ch,i,j])

    return needed_atoms, [xmin, xmax, ymin, ymax, zmin, zmax]

def Cubics(model,needed_atoms,minmax,side=2.5):

    cubics = [{},{}] # [rna, protein]

    xmin, xmax, ymin, ymax, zmin, zmax = minmax[0], minmax[1], minmax[2], minmax[3], minmax[4], minmax[5]

    length = xmax - xmin
    width  = ymax - ymin
    height = zmax - zmin

    #side = 2.5

    n, m, k = int(length//side), int(width//side), int(height//side)

    bitlen = length/n
    bitwid = width/m
    bithei = height/k

    for atom in needed_atoms:

        Type = {'RNA':0,'Protein':1}[atom[0]]

        at = model.chains[atom[2]]['RES'][atom[3]]['ATOMS'][atom[4]]

        xnum = int((at['X']-xmin)/bitlen-0.001)+1
        ynum = int((at['Y']-ymin)/bitwid-0.001)+1
        znum = int((at['Z']-zmin)/bithei-0.001)+1

        if (xnum,ynum,znum) not in cubics[Type]: cubics[Type][(xnum,ynum,znum)] = []

        cubics[Type][(xnum,ynum,znum)].append([atom[2],atom[3],atom[4]])

    return cubics

def Power(model,ch1,resind1,atind1,ch2,resind2,atind2,dist):

    def Teta(bor,atom1,atom2):

        v1 = [atom1['X']-bor[0],     atom1['Y']-bor[1],     atom1['Z']-bor[2]]
        v2 = [atom2['X']-atom1['X'], atom2['Y']-atom1['Y'], atom2['Z']-atom1['Z']]

        numer = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
        denom = pow(v1[0]**2+v1[1]**2+v1[2]**2,0.5)*pow(v2[0]**2+v2[1]**2+v2[2]**2,0.5)
        return numer/denom

    def Pd(d,Dopt):

        if d <= Dopt: return 1
        else        : return math.exp(-6.25*(d-Dopt)**2)

    def PaI(costeta):

        cosTopt  = 0.4 #cos(1.15radian) - optimal teta for  I-type atoms
        return math.exp(-6.25*(costeta-cosTopt)**2)
        
    def PaII(costeta):

        cosTopt = 1.0 #cos(0.0radian)  - optimal teta for II-type atoms
        return math.exp(-9*(costeta-cosTopt)**2)
    
    atom1 = model.chains[ch1]['RES'][resind1]['ATOMS'][atind1]
    atom2 = model.chains[ch2]['RES'][resind2]['ATOMS'][atind2]

    try:
        Dopt = {'OO':2.7,'ON':2.9,'NO':2.9,'NN':3.0}[atom1['NAME'][0]+atom2['NAME'][0]] # optimal distances
        neighbors1temp = Maps.neighbors[atom1['RESNAME']][atom1['NAME']][:]
        neighbors2temp = Maps.neighbors[atom2['RESNAME']][atom2['NAME']][:]

    except:
        return 0, 0, 0

    neighbors1 = []
    neighbors2 = []

    for i in range(len(neighbors1temp)):
        if neighbors1temp[i] in model.chains[ch1]['RES'][resind1]['ATOMNAMEDICT']:
            neighbors1.append(model.chains[ch1]['RES'][resind1]['ATOMNAMEDICT'][neighbors1temp[i]])
    for i in range(len(neighbors2temp)):
        if neighbors2temp[i] in model.chains[ch2]['RES'][resind2]['ATOMNAMEDICT']:
            neighbors2.append(model.chains[ch2]['RES'][resind2]['ATOMNAMEDICT'][neighbors2temp[i]])

    del neighbors1temp,neighbors2temp

    if len(neighbors1) == 0 or len(neighbors2) == 0: return 0.0, 0.0, 0.0

    atomtypes = ['', '']

    if len(neighbors1) == 1:
        bor1 = model.chains[ch1]['RES'][resind1]['ATOMS'][neighbors1[0]]
        bor1 = [bor1['X'],bor1['Y'],bor1['Z']]
        atomtypes[0] = 'I'
    else:
        bor11 = model.chains[ch1]['RES'][resind1]['ATOMS'][neighbors1[0]]
        bor12 = model.chains[ch1]['RES'][resind1]['ATOMS'][neighbors1[1]]
        bor1  = [(bor11['X']+bor12['X'])/2, (bor11['Y']+bor12['Y'])/2, (bor11['Z']+bor12['Z'])/2]
        atomtypes[0] = 'II'
    if len(neighbors2) == 1:
        bor2 = model.chains[ch2]['RES'][resind2]['ATOMS'][neighbors2[0]]
        bor2 = [bor2['X'],bor2['Y'],bor2['Z']]
        atomtypes[1] = 'I'
    else:
        bor21 = model.chains[ch2]['RES'][resind2]['ATOMS'][neighbors2[0]]
        bor22 = model.chains[ch2]['RES'][resind2]['ATOMS'][neighbors2[1]]
        bor2  = [(bor21['X']+bor22['X'])/2, (bor21['Y']+bor22['Y'])/2, (bor21['Z']+bor22['Z'])/2]
        atomtypes[1] = 'II'
        
    costeta1 = Teta(bor1,atom1,atom2)
    costeta2 = Teta(bor2,atom2,atom1)

    pa1 = {'I':PaI,'II':PaII}[atomtypes[0]](costeta1)
    pa2 = {'I':PaI,'II':PaII}[atomtypes[1]](costeta2)
    pd  = Pd(dist,Dopt)

    power = pd*pa1*pa2

    return power,costeta1,costeta2

def Search(model,cubics,MaxDist=3.5,MinPower=0.1):

    #MaxDist  = 3.5
    #MinPower = 0.1 

    def AdjCubes(a,b,c):

        return [(a-1,b-1,c-1),(a-1,b-1,c),(a-1,b-1,c+1),(a-1,b,c-1),(a-1,b,c),(a-1,b,c+1),(a-1,b+1,c-1),(a-1,b+1,c),(a-1,b+1,c+1),
                (a,  b-1,c-1),(a,  b-1,c),(a,  b-1,c+1),(a,  b,c-1),(a,  b,c),(a,  b,c+1),(a,  b+1,c-1),(a,  b+1,c),(a,  b+1,c+1),
                (a+1,b-1,c-1),(a+1,b-1,c),(a+1,b-1,c+1),(a+1,b,c-1),(a+1,b,c),(a+1,b,c+1),(a+1,b+1,c-1),(a+1,b+1,c),(a+1,b+1,c+1)]

    def DistSquare(atom,atom2):

        at = model.chains[atom[0]]['RES'][atom[1]]['ATOMS'][atom[2]]
        at2 = model.chains[atom2[0]]['RES'][atom2[1]]['ATOMS'][atom2[2]]

        return (at['X']-at2['X'])**2 + (at['Y']-at2['Y'])**2 + (at['Z']-at2['Z'])**2

    if len(cubics[0]) < len(cubics[1]):
        num = 0 # by RNA atoms (if their number less than protein atoms number)
    else:
        num = 1 # by Protein atoms (else)

    potentials = {}

    for (a,b,c) in cubics[num]:

        for atom in cubics[num][(a,b,c)]:

            for (d,e,f) in AdjCubes(a,b,c):

                if (d,e,f) in cubics[1-num]:

                    for atom2 in cubics[1-num][(d,e,f)]:
                            
                        distsquare = DistSquare(atom,atom2) 

                        if distsquare <= MaxDist**2:

                            dist = pow(distsquare,0.5)

                            power,teta1,teta2 = Power(model,atom[0],atom[1],atom[2],atom2[0],atom2[1],atom2[2],dist)

                            if power >= MinPower:

                                if (atom[0],atom[1],atom[2]) not in potentials: potentials[(atom[0],atom[1],atom[2])] = []

                                potentials[(atom[0],atom[1],atom[2])].append([atom2,power,dist,teta1,teta2])
    return potentials,num

'''def Clean(model,potentials):
    
    def Run(potentials):

        for (ch1,res1,at1) in potentials:

            a = model.chains[ch1]['RES'][res1]['ATOMS'][at1]
        
            if a['NAME'] in Maps.acceptors[a['RESNAME']] and len(potentials[(ch1,res1,at1)]) > 1:

                potentials[(ch1,res1,at1)] = sorted(potentials[(ch1,res1,at1)], key = lambda x: x[1], reverse=True)[:1]
        
            elif a['NAME'] in Maps.donors[a['RESNAME']] and len(potentials[(ch1,res1,at1)]) > Maps.donors[a['RESNAME']][a['NAME']]:

                number = Maps.donors[a['RESNAME']][a['NAME']]
                potentials[(ch1,res1,at1)] = sorted(potentials[(ch1,res1,at1)], key = lambda x: x[1], reverse=True)[:number]

        return potentials

    def Reverse(potentials):

        potentials2 = {}

        for (ch1,res1,at1) in potentials:

            for pot in potentials[(ch1,res1,at1)]:

                atom2 = (pot[0][0],pot[0][1],pot[0][2])
                
                if atom2 not in potentials2: potentials2[atom2] = []

                potentials2[atom2].append([[ch1,res1,at1],] + pot[1:])

        return potentials2

    for i in (1,2):
        potentials = Run(potentials)
        potentials = Reverse(potentials)

    return potentials'''

def Atompairs(model,side=2.5,MaxDist=3.5,MinPower=0.1,allin=False):

    needed_atoms, minmax = Process(model,allin)
    cubics = Cubics(model,needed_atoms,minmax,side)
    del needed_atoms,minmax

    potentials,bypass = Search(model,cubics,MaxDist,MinPower)
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

    monopairs = {} # (ch1,res1,ch2,res2): [atompair_ids]

    for ap in atompairs:

        ch1,   ch2 = ap['ADRESS'][0][0], ap['ADRESS'][1][0]
        res1, res2 = ap['ADRESS'][0][1], ap['ADRESS'][1][1]

        if (ch1,res1,ch2,res2) not in monopairs: monopairs[(ch1,res1,ch2,res2)] = []

        monopairs[(ch1,res1,ch2,res2)].append(ap['ID'])

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

    if 'P' in model.headers['TYPE']:
        atompairs           = Atompairs(model)
        monopairs,atompairs = Monopairs(model,atompairs)
    else:
        atompairs = []
        monopairs = []

    model.atompairs = atompairs
    model.monopairs = monopairs




    
