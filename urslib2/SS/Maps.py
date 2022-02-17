
'''donors = {'A'  : {"O2'":1, 'N6':2},
          'C'  : {"O2'":1, 'N4':2},
          'G'  : {"O2'":1, 'N1':1, 'N2':2},
          'U'  : {"O2'":1, 'N3':1},
          'ALA': {  'N':1,},
          'ARG': {  'N':1, 'NE':1,'NH1':2,'NH2':2},
          'ASN': {  'N':1,'ND2':2},
          'ASP': {  'N':1,'OD2':1},
          'CYS': {  'N':1,},
          'GLN': {  'N':1,'NE2':2},
          'GLU': {  'N':1,'OE2':1},
          'GLY': {  'N':1,},
          'HIS': {  'N':1,'ND1':1,'NE2':1},
          'ILE': {  'N':1,},
          'LEU': {  'N':1,},
          'LYS': {  'N':1, 'NZ':3},
          'MET': {  'N':1,},
          'PHE': {  'N':1,},
          'PRO': {  'N':1,},
          'SER': {  'N':1, 'OG':1},
          'THR': {  'N':1,'OG1':1},
          'TRP': {  'N':1,'NE1':1},
          'TYR': {  'N':1, 'OH':1},
          'VAL': {  'N':1,}}'''

donors_acceptors = {'A'  : ['OP1','OP2','N7','N1','N3',"O2'",'N6',"O3'","O4'","O5'"],
                    'C'  : ['OP1','OP2','O2','N3',"O2'",'N4',"O3'","O4'","O5'"],
                    'G'  : ['OP1','OP2','N7','O6','N3',"O2'",'N1','N2',"O3'","O4'","O5'"],
                    'U'  : ['OP1','OP2','O2','O4',"O2'",'N3',"O3'","O4'","O5'"],
                    '5BU': ['OP1','OP2','O2','O4',"O2'",'N3',"O3'","O4'","O5'"], ##mod
                    'ALA': ['O','N'],
                    'ARG': ['O','N','NE','NH1','NH2'],
                    'ASN': ['O','OD1','N','ND2'],
                    'ASP': ['O','OD1','N','OD2'],
                    'CYS': ['O','N'],
                    'GLN': ['O','OE1','N','NE2'],
                    'GLU': ['O','OE1','N','OE2'],
                    'GLY': ['O','N'],
                    'HIS': ['O','N','ND1','NE2'],
                    'ILE': ['O','N'],
                    'LEU': ['O','N'],
                    'LYS': ['O','N','NZ'],
                    'MET': ['O','N'],
                    'PHE': ['O','N'],
                    'PRO': ['O','N'],
                    'SER': ['O','N','OG'],
                    'THR': ['O','N','OG1'],
                    'TRP': ['O','N','NE1'],
                    'TYR': ['O','N','OH'],
                    'VAL': ['O','N']}

neighbors = {'A'  :{'OP1': ['P'],
                    'OP2': ['P'],
                    "O2'": ["C2'"],
                    "O3'": ["C3'"],
                    "O4'": ["C4'","C1'"],
                    "O5'": ["C5'", 'P'],
                    'N1' : ['C6', 'C2'],
                    'N3' : ['C2', 'C4'],
                    'N6' : ['C6'],
                    'N7' : ['C8', 'C5']},
             'C'  :{'OP1': ['P'],
                    'OP2': ['P'],
                    "O2'": ["C2'"],
                    "O3'": ["C3'"],
                    "O4'": ["C4'","C1'"],
                    "O5'": ["C5'", 'P'],
                    'O2' : ['C2'],
                    'N3' : ['C2', 'C4'],                  
                    'N4' : ['C4']},
             'G'  :{'OP1': ['P'],
                    'OP2': ['P'],
                    "O2'": ["C2'"],
                    "O3'": ["C3'"],
                    "O4'": ["C4'","C1'"],
                    "O5'": ["C5'", 'P'],
                    'N1' : ['C6', 'C2'],
                    'N2' : ['C2'],
                    'N3' : ['C2', 'C4'],
                    'O6' : ['C6'],
                    'N7' : ['C8', 'C5']},
             'U'  :{'OP1': ['P'],
                    'OP2': ['P'],
                    "O2'": ["C2'"],
                    "O3'": ["C3'"],
                    "O4'": ["C4'","C1'"],
                    "O5'": ["C5'", 'P'],
                    'O2' : ['C2'],
                    'N3' : ['C2', 'C4'],
                    'O4' : ['C4']},
             '5BU':{'OP1': ['P'],  ##mod
                    'OP2': ['P'],
                    "O2'": ["C2'"],
                    "O3'": ["C3'"],
                    "O4'": ["C4'","C1'"],
                    "O5'": ["C5'", 'P'],
                    'O2' : ['C2'],
                    'N3' : ['C2', 'C4'],
                    'O4' : ['C4']},
             'ALA':{'O'  : ['C'],
                    'N'  : ['CA']},
             'ARG':{'O'  : ['C'],
                    'N'  : ['CA'],
                    'NE' : ['CD', 'CZ'],
                    'NH1': ['CZ'],
                    'NH2': ['CZ']},
             'ASN':{'O'  : ['C'],
                    'N'  : ['CA'],
                    'OD1': ['CG'],
                    'ND2': ['CG']},
             'ASP':{'O'  : ['C'],
                    'N'  : ['CA'],
                    'OD1': ['CG'],
                    'OD2': ['CG']},
             'CYS':{'O'  : ['C'],
                    'N'  : ['CA']},
             'GLN':{'O'  : ['C'],
                    'N'  : ['CA'],
                    'OE1': ['CD'],
                    'NE2': ['CD']},
             'GLU':{'O'  : ['C'],
                    'N'  : ['CA'],
                    'OE1': ['CD'],
                    'OE2': ['CD']},
             'GLY':{'O'  : ['C'],
                    'N'  : ['CA']},
             'HIS':{'O'  : ['C'],
                    'N'  : ['CA'],
                    'ND1': ['CG', 'CE1'],
                    'NE2': ['CD2', 'CE1']},
             'ILE':{'O'  : ['C'],
                    'N'  : ['CA']},
             'LEU':{'O'  : ['C'],
                    'N'  : ['CA']},
             'LYS':{'O'  : ['C'],
                    'N'  : ['CA'],
                    'NZ' : ['CE']},
             'MET':{'O'  : ['C'],
                    'N'  : ['CA']},
             'PHE':{'O'  : ['C'],
                    'N'  : ['CA']},
             'PRO':{'O'  : ['C'],
                    'N'  : ['CA', 'CD']},
             'SER':{'O'  : ['C'],
                    'N'  : ['CA'],
                    'OG' : ['CB']},
             'THR':{'O'  : ['C'],
                    'N'  : ['CA'],
                    'OG1': ['CB']},
             'TRP':{'O'  : ['C'],
                    'N'  : ['CA'],
                    'NE1': ['CD1', 'CE2']},
             'TYR':{'O'  : ['C'],
                    'N'  : ['CA'],
                    'OH' : ['CZ']},
             'VAL':{'O'  : ['C'],
                    'N'  : ['CA']}}

'''
def temp_parse():

    files = ['A.txt','C.txt','G.txt',
             'U.txt','ALA.txt','ARG.txt',
             'ASN.txt','ASP.txt','CYS.txt',
             'GLN.txt','GLU.txt','GLY.txt',
             'HIS.txt','ILE.txt','LEU.txt',
             'LYS.txt','MET.txt','PHE.txt',
             'PRO.txt','SER.txt','THR.txt',
             'TRP.txt','TYR.txt','VAL.txt',]

    x2_donors = {}
    donors = {}
    acceptors = {}
    neighbors = {}

    for f in files:
        with open('/home/baulin/eugene/Работа/RSSDB/programs/USSR/SS/residues/'+f) as file:
            for line in file:

                while '  ' in line: line = line.replace('  ',' ')

                if line[:6] == 'RESIDU':

                    res = line[:-1].split(' ')[1]
                    x2_donors[res] = []
                    donors[res] = []
                    acceptors[res] = []
                    neighbors[res] = {}

                elif line[:6] == 'CONECT':

                    line = line[:-1].split(' ')
                    atom = line[1]

                    if atom[0] in ('O','N'):

                        bors  = []
                        hydro = 0
                        
                        for n in line[3:]:
                            if not n: break
                            if n[0] == 'H': hydro += 1
                            else: bors.append(n)

                        if hydro > 0:
                            donors[res].append(atom)
                            if hydro > 1:
                                x2_donors[res].append(atom)
                        elif len(bors)<=2:
                            acceptors[res].append(atom)

                        neighbors[res][atom] = bors

    return x2_donors,donors,acceptors,neighbors
'''  
