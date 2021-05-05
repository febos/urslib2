try: from urslib2 import Tools
except ImportError: import Tools

class Model():

    def __init__(self,modelname):

        self.headers   = {'ID':           1,    # model id
                          'DATE':        '',    # date of deposition
                          'HEADER':      '',    # header of model
                          'PDBFILE':     '',    # 4 letter id of pdbfile
                          'TITLE':       '',    # title of model
                          'KEYWDS':      [],    # key words
                          'EXPDTA':      [],    # experiment data
                          'NUMMDL':       1,    # number of models
                          'MDLTYP':      [],    # details about model
                          'AUTHOR':      '',    # authors of experiment
                          'RESOL':    '\\N',    # resolution
                          'CRYSYM':      [],    # crystal symmetry info
                          'AREMODELS':False,    # file has several models
                          'MDLNO':        1,    # model No
                          'MODRES':      {},    # modified residues
                          'TYPE':        ''}    # type of model (R, RD, RP, RPD) RNA/Protein/DNA
        self.chains    = {}
        self.molecules = {}
        self.ligands   = []         # will be deleted
        self.ids       = {}
        self.parse(modelname)
        self.mark_id_to_all()

    def parse(self, modelname):

        PDBdict = {'HEADER':'', 'TITLE':  [], 'COMPND':   [], 'SOURCE':   [],
                   'KEYWDS':[], 'EXPDTA': [], 'NUMMDL':    0, 'MDLTYP':   [],
                   'AUTHOR':[], 'REMARK2':'', 'REMARK290':[], 'REMARK465':[],
                   'DBREF': [], 'SEQRES': [], 'MODRES':   [], 'MODEL':     0,
                   'ATOM':  []} 

        with open(modelname) as model:
            for line in model:

                """parsing segments of PDB file:
                   HEADER - header of file
                   TITLE  - title of file
                   COMPND, SOURCE - information on molecules
                   KEYWDS - key words
                   EXPDTA - information on experiment
                   NUMMDL - number of models
                   MDLTYP - details about model
                   AUTHOR - authors of experiment
                   REMARK2 - resolution info
                   REMARK290 - crystal symmetry info
                   REMARK465 - missing residues
                   DBREF - chain info (start, end, id-code)
                   SEQRES - sequence of residues
                   MODRES - modified residues
                   MODEL - model No
                   ATOM - coordinate section (ATOM, HETATM, TER)
                """

                if   line[:6] == 'HEADER'      : PDBdict['HEADER'] = line

                elif line[:6] == 'TITLE '      : PDBdict['TITLE'].append(line)

                elif line[:6] == 'COMPND'      : PDBdict['COMPND'].append(line)

                elif line[:6] == 'SOURCE'      : PDBdict['SOURCE'].append(line)

                elif line[:6] == 'KEYWDS'      : PDBdict['KEYWDS'].append(line)

                elif line[:6] == 'EXPDTA'      : PDBdict['EXPDTA'].append(line)

                elif line[:6] == 'NUMMDL'      : PDBdict['NUMMDL'] = int(line[10:14])

                elif line[:6] == 'MDLTYP'      : PDBdict['MDLTYP'].append(line)

                elif line[:6] == 'AUTHOR'      : PDBdict['AUTHOR'].append(line)

                elif line[:6] == 'REMARK':

                    if int(line[7:10]) == 2    : PDBdict['REMARK2'] = line

                    elif int(line[7:10]) == 290: PDBdict['REMARK290'].append(line)

                    elif int(line[7:10]) == 465: PDBdict['REMARK465'].append(line)

                elif line[:6] in ('DBREF ',
                                  'DBREF1')    : PDBdict['DBREF'].append(line)

                elif line[:6] == 'SEQRES'      : PDBdict['SEQRES'].append(line)

                elif line[:6] == 'MODRES'      : PDBdict['MODRES'].append(line)

                elif line[:6] == 'MODEL '      : PDBdict['MODEL'] = int(line[10:14])

                elif line[:6] in ('ATOM  ','HETATM','TER   ') and line[76:78] != ' H':
                    ''' we do not consider Hydrogen atoms'''
                    PDBdict['ATOM'].append(line)

        """ The end of parsing.
        The start of filling self.headers"""

        # DATE

        day,month,year = PDBdict['HEADER'][50:52],PDBdict['HEADER'][53:56],PDBdict['HEADER'][57:59]

        '''99 -> 1999; 01 -> 2001'''
        if int(year[0]) < 2: year = str(2000 + int(year))
        else               : year = str(1900 + int(year))

        '''MAR -> 03'''
        month = Tools.months[month] 

        '''assembling of date of mysql format'''
        self.headers['DATE'] = year + '-' + month + '-' + day

        ''' delete garbage '''
        del day,month,year 

        # HEADER

        """ clean_seq: ' TRNA    FRAGMENT    ' -> 'TRNA FRAGMENT'"""
        self.headers['HEADER']   = Tools.clean_seq(PDBdict['HEADER'][10:50])

        # PDBFILE
        
        self.headers['PDBFILE']  = PDBdict['HEADER'][62:66].upper()

        del PDBdict['HEADER']

        # TITLE

        Title = ''

        for string in PDBdict['TITLE']:
            Title += string[10:80]

        self.headers['TITLE'] = Tools.clean_seq(Title)

        del Title,PDBdict['TITLE']

        # COMPND,SOURCE (self.molecules)
        
        Mol_id = 0
        icolon = 0
        token  = ''
        value  = ''

        for string in PDBdict['COMPND']:

            icolon = string.find(':') # index of char ':' in string

            if icolon != -1: # if there is ':' in string

                if string[9] == ' ': token = string[10:icolon]
                else               : token = string[11:icolon]
                '''format of COMPND = "TOKEN: VALUE;"'''
                value = string[icolon+1:80]

                if token == 'MOL_ID':

                    Mol_id = int(Tools.cut_for_int(value))
                    self.molecules[Mol_id] = Tools.Molecule(Mol_id)

                else:
                    token = token.replace('EXSPRESSION','EXPRESSION')
                    value = Tools.clean_source(value)
                    self.molecules[Mol_id][token] = value

            else: # if our string is a continuance of previous

                value = Tools.clean_source(string[11:80])
                self.molecules[Mol_id][token] += value

        ''' SOURCE has the same format as COMPND'''

        for string in PDBdict['SOURCE']:

            icolon = string.find(':')

            if icolon != -1:

                if string[9] == ' ': token = string[10:icolon]
                else               : token = string[11:icolon]

                if token == 'FRAGMENT': continue

                value = string[icolon+1:80]

                if token == 'MOL_ID': Mol_id = int(Tools.cut_for_int(value))

                else:
                    token = token.replace('EXSPRESSION','EXPRESSION')
                    value = Tools.clean_source(value)
                    self.molecules[Mol_id][token] = value

            else:

                value = Tools.clean_source(string[11:80])
                self.molecules[Mol_id][token] += value
        
        del Mol_id,icolon,token,value,PDBdict['COMPND'],PDBdict['SOURCE']

        # KEYWDS

        keywords = ''

        for keyword in PDBdict['KEYWDS']: keywords += keyword[10:80]

        keywords = Tools.clean_seq(keywords).split(', ') # 'one; two;three' -> ['one','two','three']

        self.headers['KEYWDS'] = keywords

        del keywords,PDBdict['KEYWDS']

        # EXPDTA

        expdata = ''

        for datum in PDBdict['EXPDTA']: expdata += datum[10:80]

        expdata = Tools.clean_seq(expdata).split(';')

        self.headers['EXPDTA'] = expdata

        del expdata,PDBdict['EXPDTA']        

        # NUMMDL

        if PDBdict['NUMMDL']:

            self.headers['NUMMDL'] = PDBdict['NUMMDL']

        del PDBdict['NUMMDL']

        # MDLTYP

        lines = ''

        for line in PDBdict['MDLTYP']: lines += line[10:80]

        if lines: self.headers['MDLTYP'] = Tools.clean_seq(lines).split(';')

        del lines,PDBdict['MDLTYP']

        # AUTHOR

        lines = ''

        for line in PDBdict['AUTHOR']: lines += line[10:80]

        if lines: self.headers['AUTHOR'] = Tools.clean_seq(lines)

        del lines,PDBdict['AUTHOR']

        # RESOL (resolution in angstroms)

        resstring = PDBdict['REMARK2']

        if resstring[11:38] != 'RESOLUTION. NOT APPLICABLE.':

            self.headers['RESOL'] = float(Tools.cut_for_float(resstring[23:30]))

        del resstring,PDBdict['REMARK2']

        # CRYSYM (Crystallographic Symmetry, REMARK290)

        self.headers['CRYSYM'] = PDBdict['REMARK290']

        del PDBdict['REMARK290']

        # MDLNO (Model No)

        if PDBdict['MODEL']:

            self.headers['AREMODELS'] = True
            self.headers['MDLNO'] = PDBdict['MODEL']

        del PDBdict['MODEL']    

        # MODRES

        mod, std = 0,0

        for modres in PDBdict['MODRES']:

            mod = Tools.cut_spaces(modres[12:15]) # real name of residue
            std = Tools.cut_spaces(modres[24:27]) # standard name of residue
            self.headers['MODRES'][mod] = std

        del mod,std,PDBdict['MODRES']

        """ The end of filling self.headers.
        The start of filling self.chains"""

        # Chain Start and End  (DBREF)

        chain, start, end = '',0,0
         
        for line in PDBdict['DBREF']:

            chain = line[12]
            start = int(line[14:18])
            end   = int(line[20:24])

            ''' next IF is for situation like:
            DBREF A 1 1
            DBREF A 2 20'''
            if chain not in self.chains:

                self.chains[chain] = Tools.Chain(chain)
                self.chains[chain]['START']  = start

            self.chains[chain]['END']    = end

        del chain,start,end,PDBdict['DBREF']

        # SEQRES (Sequence of Residues)

        chain = ''
        
        for line  in PDBdict['SEQRES']:

            chain = line[11]
            if chain not in self.chains: self.chains[chain] = Tools.Chain(chain) # if there is not DBREF for this chain
            self.chains[chain]['LENGTH'] = int(line[13:17])
            self.chains[chain]['SEQ']   += line[17:80]

        for ch in sorted(list(self.chains.keys())):

            """ '  A, G, C, CCC  ' -> ['A','G','C','CCC']"""
            self.chains[ch]['SEQ'] = Tools.clean_seq(self.chains[ch]['SEQ']).split(' ')

        del chain,PDBdict['SEQRES']

        # filling self.chains[chain]['MOL_ID']

        for ch in sorted(list(self.chains.keys())):

            for mol in self.molecules:

                if ch in self.molecules[mol]['CHAIN']:

                    self.chains[ch]['MOL_ID'] = mol
                    break

        # Marking TYPE of chain

        check, stat, Max, result = 0, {}, 0, ''

        for ch in sorted(list(self.chains.keys())):

            check  = min(len(self.chains[ch]['SEQ']),11) # we look at first CHECK residues (if len(SEQ) > CHECK)
            stat   = {'RNA':0,'DNA':0,'Protein':0,'Unknown':0}

            for residue in self.chains[ch]['SEQ'][:check]:

                stat[Tools.restype[residue]] += 1

            Max = max(stat.values())

            result = str(int(stat['RNA']     == Max))+\
                     str(int(stat['DNA']     == Max))+\
                     str(int(stat['Protein'] == Max))+\
                     str(int(stat['Unknown'] == Max))
            
            if   result[0] == '1': self.chains[ch]['TYPE'] = 'RNA'
            elif result[1] == '1': self.chains[ch]['TYPE'] = 'DNA'
            elif result[2] == '1': self.chains[ch]['TYPE'] = 'Protein'
            elif result[3] == '1': self.chains[ch]['TYPE'] = 'Unknown'

            else:
                print('Error. Can\'t mark a type of chain %s. Code == %s'%(ch,result))
                exit()

        del check,stat,Max,result

        # TYPE (type of model)

        r,d,p = '','',''

        for ch in sorted(list(self.chains.keys())):

            if self.chains[ch]['TYPE']   == 'RNA'    : r = 'R'
            elif self.chains[ch]['TYPE'] == 'Protein': p = 'P'
            elif self.chains[ch]['TYPE'] == 'DNA'    : d = 'D'

        self.headers['TYPE'] = r+p+d

        del r,p,d

        # MISRES (Missing Residues, REMARK 465)

        chain, residue, number = '', '', ''
        firstm, lastm, start   = 0, 0, None

        if PDBdict['REMARK465']:

            if 'SSEQI' in PDBdict['REMARK465'][5]   : start = 6
            elif 'SSEQI' in PDBdict['REMARK465'][6] : start = 7
        
            ''' next three conditions are for three different formats'''
            if self.headers['AREMODELS'] and 'MODELS' in PDBdict['REMARK465'][5]:

                firstm = int(PDBdict['REMARK465'][5][19:21])
                lastm  = int(PDBdict['REMARK465'][5][22:26])
            
                for res in PDBdict['REMARK465'][start:]:

                    if firstm <= self.headers['MDLNO'] <= lastm:

                        chain, residue = res[19], Tools.cut_spaces(res[15:18])
                        number = Tools.pdbnum(res[21:27])
                        self.chains[chain]['MISRES'].append([residue, number])

            elif self.headers['AREMODELS']:

                for res in PDBdict['REMARK465'][start:]:

                    if  self.headers['MDLNO'] == res[12:15]:

                        chain, residue = res[19], Tools.cut_spaces(res[15:18])
                        number = Tools.pdbnum(res[21:27])
                        self.chains[chain]['MISRES'].append([residue, number])
            else:
                
                for res in PDBdict['REMARK465'][start:]:

                    chain, residue = res[19], Tools.cut_spaces(res[15:18])
                    number = Tools.pdbnum(res[21:27])
                    self.chains[chain]['MISRES'].append([residue, number])

        del chain,residue,number,firstm,lastm,start,PDBdict['REMARK465']

        # ATOM

        atom,chain = None, None
        ter  = {}
        seen = {}
        for ch in sorted(list(self.chains.keys())): ter[ch] = False

        for line in PDBdict['ATOM']:

            if line[:6] in ('ATOM  ','HETATM'):

                atom = Tools.Atom(line)

                chain = atom['CHAIN']

                if chain not in self.chains: # if atom don't belong to any chain

                    if chain + atom['RESNUM'] not in seen:

                        seen[chain + atom['RESNUM']] = 1
                        res = Tools.Residue(atom)

                        if Tools.restype[res['NAME']] in ('RNA','DNA','Protein'):

                            atom['TYPE'] = 'Unknown'
                            res['TYPE']  = 'Unknown'

                        else:

                            atom['TYPE'] = Tools.restype[res['NAME']]
                            res['TYPE']  = Tools.restype[res['NAME']]

                        res['ATOMS'].append(atom)
                        self.ligands.append(res)

                    else:
                        
                        atom['TYPE'] = self.ligands[-1]['TYPE']
                        self.ligands[-1]['ATOMS'].append(atom)

                elif not ter[chain]: # if we are still in sequence of chain

                    if chain + atom['RESNUM'] not in seen:

                        seen[chain + atom['RESNUM']] = 1
                        res = Tools.Residue(atom)

                        # if .NA = .NA or .rotein = .rotein (for equality of RNA and DNA residues)
                        if self.chains[chain]['TYPE'][1:] == Tools.restype[res['NAME']][1:]:

                            atom['TYPE'] = self.chains[chain]['TYPE']
                            res['TYPE']  = self.chains[chain]['TYPE']

                        else:

                            atom['TYPE'] = Tools.restype[res['NAME']]
                            res['TYPE']  = Tools.restype[res['NAME']]

                        res['ATOMS'].append(atom)
                        self.chains[chain]['RES'].append(res)

                    else:
                        
                        atom['TYPE'] = self.chains[chain]['RES'][-1]['TYPE']
                        self.chains[chain]['RES'][-1]['ATOMS'].append(atom)

                else: # if we are after sequence of chain (various ligands)

                    if chain + atom['RESNUM'] not in seen:

                        seen[chain + atom['RESNUM']] = 1
                        res = Tools.Residue(atom)

                        atom['TYPE'] = Tools.restype[res['NAME']]
                        res['TYPE']  = Tools.restype[res['NAME']]

                        res['ATOMS'].append(atom)
                        self.chains[chain]['LIGANDS'].append(res)

                    else:
                        
                        atom['TYPE'] = self.chains[chain]['LIGANDS'][-1]['TYPE']
                        self.chains[chain]['LIGANDS'][-1]['ATOMS'].append(atom)

            else: ter[chain] = True

        del atom,chain,ter,seen,PDBdict
                        
        ### Garbage chains (that are not in SEQRES nor COMPND) ###

        added, ch = [], None

        # Looking for garbage chains

        for lig in self.ligands:

            ch = lig['CHAIN']

            if ch not in self.chains: self.chains[ch] = Tools.Chain(ch)

        for i in range(len(self.ligands)):

            ch = self.ligands[i]['CHAIN']

            if ch in self.chains: self.chains[ch]['RES'].append(self.ligands[i])

        # START, END, LENGTH, SEQ and SEQ2 for garbage chains
        for ch in sorted(list(self.chains.keys())):

            if self.chains[ch]['TYPE'] != '\\N': continue

            self.chains[ch]['GARBAGE'] = True

            for i in range(0,len(self.chains[ch]['RES'])):

                self.chains[ch]['NUMS'][self.chains[ch]['RES'][i]['PDBNUM']] = i

            self.chains[ch]['START']   = int(self.chains[ch]['RES'][0]['FLOAT'])
            self.chains[ch]['END']     = int(self.chains[ch]['RES'][-1]['FLOAT'])
            self.chains[ch]['LENGTH']  = len(self.chains[ch]['RES'])
            self.chains[ch]['SEQ2']    = [r['NAME'] for r in self.chains[ch]['RES']]
            self.chains[ch]['SEQ']     = self.chains[ch]['SEQ2']
            self.chains[ch]['LIGSEQ']  = [r['NAME'] for r in self.chains[ch]['LIGANDS']]

        check, stat, Max, result = 0, {}, 0, ''

        # Type of garbage chains
        for ch in sorted(list(self.chains.keys())):

            if self.chains[ch]['TYPE'] != '\\N': continue

            check  = min(len(self.chains[ch]['SEQ']),11) # we look at first CHECK residues (if len(SEQ) > CHECK)
            stat   = {'RNA':0,'DNA':0,'Protein':0,'Unknown':0, 'Water': 0, 'Metal': 0}

            for residue in self.chains[ch]['SEQ'][:check]: stat[Tools.restype[residue]] += 1

            stat['Unknown'] = stat['Unknown'] + stat['Metal'] + stat['Water']

            Max = max(stat.values())

            result = str(int(stat['RNA']     == Max))+\
                     str(int(stat['DNA']     == Max))+\
                     str(int(stat['Protein'] == Max))+\
                     str(int(stat['Unknown'] == Max))
            
            if   result[0] == '1': self.chains[ch]['TYPE'] = 'RNA'
            elif result[1] == '1': self.chains[ch]['TYPE'] = 'DNA'
            elif result[2] == '1': self.chains[ch]['TYPE'] = 'Protein'
            elif result[3] == '1': self.chains[ch]['TYPE'] = 'Unknown'

            else:
                print('Error. Can\'t mark a type of chain %s. Code == %s'%(ch,result))
                exit()

        # Type of residues and atoms in garbage chains
        for ch in sorted(list(self.chains.keys())):

            if self.chains[ch]['GARBAGE']:

                for res in self.chains[ch]['RES']:

                    if self.chains[ch]['TYPE'][1:] == Tools.restype[res['NAME']][1:]:

                        res['TYPE'] = self.chains[ch]['TYPE']
                        for atom in res['ATOMS']: atom['TYPE'] = self.chains[ch]['TYPE']

                    else:

                        res['TYPE']  = Tools.restype[res['NAME']]
                        for atom in res['ATOMS']: atom['TYPE'] = Tools.restype[res['NAME']]

        del added,ch,check,stat,Max,result, self.ligands

        # Inserting MISRES

        misseq, misres, number, name, pnum, icode = None, None, None, None, None, None

        for ch in sorted(list(self.chains.keys())):

            misseq = self.chains[ch]['MISRES']

            while misseq:

                misres = misseq.pop(0)
                name   = misres[0]
                number = misres[1]
                icode  = number[-1]
                if icode == ' ': icode = ''
                pnum   = int(number[:-1])
                misres = Tools.MisResidue()
                misres['NAME']   = name
                misres['PDBNUM'] = number
                misres['FLOAT']  = Tools.pdbnum_to_float(number)
                misres['CHAIN']  = ch
                misres['DSSR']   = ch + '.' + name + '.' + str(pnum) + '.' + icode

                if Tools.restype[misres['NAME']] == self.chains[ch]['TYPE']:

                    misres['TYPE'] = self.chains[ch]['TYPE']

                else: misres['TYPE'] = 'Unknown'

                # The order below is important

                if not self.chains[ch]['RES'] or misres['FLOAT'] > self.chains[ch]['RES'][-1]['FLOAT']:

                    self.chains[ch]['RES'].append(misres)
                    continue

                if misres['FLOAT'] < self.chains[ch]['RES'][0]['FLOAT']:

                    self.chains[ch]['RES'].insert(0,misres)
                    continue

                for i in range(1,len(self.chains[ch]['RES'])):

                    if self.chains[ch]['RES'][i-1]['FLOAT'] <\
                       misres['FLOAT'] <\
                       self.chains[ch]['RES'][i]['FLOAT']:

                        ''' next if is for situation like
                        1503 -> 1500 2003 1501 1502 1504'''
                        if i+1 == len(self.chains[ch]['RES']) or\
                           self.chains[ch]['RES'][i-1]['FLOAT'] <\
                           misres['FLOAT'] <\
                           self.chains[ch]['RES'][i+1]['FLOAT']:

                            self.chains[ch]['RES'].insert(i,misres)
                            break

        del misres,number,name,pnum,icode

        # Filling self.chains[ch]['NUMS'] + adding new START,END,LENGTH and SEQ2

        for ch in sorted(list(self.chains.keys())):

            for i in range(0,len(self.chains[ch]['RES'])):

                self.chains[ch]['NUMS'][self.chains[ch]['RES'][i]['PDBNUM']] = i

            self.chains[ch]['START']   = int(self.chains[ch]['RES'][0]['FLOAT'])
            self.chains[ch]['END']     = int(self.chains[ch]['RES'][-1]['FLOAT'])
            self.chains[ch]['LENGTH']  = len(self.chains[ch]['RES'])
            self.chains[ch]['SEQ2']    = [r['NAME'] for r in self.chains[ch]['RES']]
            self.chains[ch]['LIGSEQ']  = [r['NAME'] for r in self.chains[ch]['LIGANDS']]
        
    def mark_id_to_all(self):

        self.ids = {'CHAIN' : 1, 'NUCL': 1, 'AMINO':  1, 'LIGAND': 1, 'ATOM': 1}

        chaintype, restype, atomtype = None, None, None 

        for ch in sorted(list(self.chains.keys())):

            chaintype = self.chains[ch]['TYPE'] 

            self.chains[ch]['ID'] = self.ids['CHAIN']
            self.ids['CHAIN'] += 1

            for i in range(0,len(self.chains[ch]['RES'])):

                restype = self.chains[ch]['RES'][i]['TYPE']

                if restype in ('RNA','DNA'):

                    self.chains[ch]['RES'][i]['ID'] = self.ids['NUCL']
                    self.ids['NUCL'] += 1

                elif restype == 'Protein':

                    self.chains[ch]['RES'][i]['ID'] = self.ids['AMINO']
                    self.ids['AMINO'] += 1

                elif restype == 'Unknown':

                    self.chains[ch]['RES'][i]['ID'] = self.ids['LIGAND']
                    self.ids['LIGAND'] += 1

                for j in range(0,len(self.chains[ch]['RES'][i]['ATOMS'])):

                    self.chains[ch]['RES'][i]['ATOMS'][j]['ID'] = self.ids['ATOM']
                    self.ids['ATOM'] += 1

            for i in range(0,len(self.chains[ch]['LIGANDS'])):

                ligtype = self.chains[ch]['LIGANDS'][i]['TYPE']

                if ligtype in ('RNA','DNA'):

                    self.chains[ch]['LIGANDS'][i]['ID'] = self.ids['NUCL']
                    self.ids['NUCL'] += 1

                elif ligtype == 'Protein':

                    self.chains[ch]['LIGANDS'][i]['ID'] = self.ids['AMINO']
                    self.ids['AMINO'] += 1

                elif ligtype == 'Unknown':

                    self.chains[ch]['LIGANDS'][i]['ID'] = self.ids['LIGAND']
                    self.ids['LIGAND'] += 1

                for j in range(0,len(self.chains[ch]['LIGANDS'][i]['ATOMS'])):

                    self.chains[ch]['LIGANDS'][i]['ATOMS'][j]['ID'] = self.ids['ATOM']
                    self.ids['ATOM'] += 1

        for Id in self.ids: self.ids[Id] -= 1

        del chaintype,restype,atomtype

