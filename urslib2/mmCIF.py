#import os
#import PDB
import glob

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
                          'TYPE':        '',    # type of model (R, RD, RP, RPD) RNA/Protein/DNA
                          'LIGLIST':     '',    # ','.join(ligands)
                          'METLIST':     '',    # ','.join(metals)
                          'MASKEDCHS':   {},    # {fake_chain: designated_chain}
                          'CHAINBIO':    {},}    
        self.chains     = {}
        self.molecules  = {}
        self.ids        = {}
        self.restype    = {}         # residue: type
        self.mmcif_dict = {}
        self.allwords   = {}         # for simplesearch
        self.CIFparse(modelname)
        self.CIFprocess(modelname)
        self.mark_id_to_all()
        del self.mmcif_dict

    def CIFparse(self, modelname):

        #parser     = MMCIFParser()
        #self.mmcif_dict = MMCIF2Dict.MMCIF2Dict(modelname)
        #self.structure  = parser.get_structure('PHA-L', modelname)

        def FightQuotes(line):

            line2 = ''
            count = ''

            for ch in line:
                if ch in ("'",'"'):
                    line2 += ch
                    if not count: count = ch
                    elif count==ch: count = ''
                elif ch in (' ', '\t'):
                    if count: line2 += '&&&&&'
                    else    : line2 += ' '
                else: line2 += ch
            return line2

        mmcif_dict = {} # key:value or table:[[headers],[values]*n]

        intable = False # after loop_ until #
        intext  = False # we are between ; and ;
        inrow   = False # after loop_ and headers or just before ; for simple-value
        current = None  # last non-finished entity (simple-value or table)

        with open(modelname) as model:
            for line in model:

                if ("'" in line and line.count("'")>1) or \
                   ('"' in line and line.count('"')>1): line = FightQuotes(line)

                if   line[0]  == '#'    : intable = False
                elif line[:5] == 'loop_': intable = True

                elif intext:

                    text = line[:-1]

                    if not text: continue

                    if text[0] == ';': text = text[1:]

                    if not text:
                        intext = False
                        if intable and len(mmcif_dict[current][-1]) == len(mmcif_dict[current][0]):
                            inrow = False
                    else:
                        if intable:
                            mmcif_dict[current][-1][-1] += text
                        else      : mmcif_dict[current] += text
                elif inrow:

                    if line[0] == ';':

                        if intable:
                            mmcif_dict[current][-1] += [line[1:-1]]
                        else:
                            mmcif_dict[current] += line[1:-1]
                            inrow = False
                        intext = True

                    elif not intable:
                        mmcif_dict[current] = line[:-1]
                        inrow = False

                    else:
                        row = line.split()
                        mmcif_dict[current][-1] += row
                        if len(mmcif_dict[current][-1]) == len(mmcif_dict[current][0]):
                            inrow = False

                elif line[0]  == '_':

                    string = line[:-1].split()

                    if intable:

                        table,row = string[0].split('.')
                        if table not in mmcif_dict:
                            mmcif_dict[table] = [[]]
                            current = table
                        mmcif_dict[table][0].append(row)

                    else:
                        if len(string) == 2:
                            mmcif_dict[string[0]] = string[1]
                        else:
                            mmcif_dict[string[0]] = ''
                            current = string[0]
                            inrow = True

                elif intable:

                    if line[0] == ';':
                        mmcif_dict[current].append([line[1:-1]])
                        intext,inrow = True,True
                    else:
                        row = line.split()

                        if len(row) == len(mmcif_dict[current][0]):
                            mmcif_dict[current].append(row)
                        else:
                            mmcif_dict[current].append(row)
                            inrow = True

        del intable,inrow,intext,current,table,row,string
        self.mmcif_dict = mmcif_dict

    def CIFprocess(self,modelname):

        def CleanText(text):
            text = text.replace('&&&&&',' ')
            if text[0]==text[-1]=="'": text = text[1:-1]
            if text[0]==text[-1]=='"': text = text[1:-1]
            while text[0]  in (' ','_',';'): text = text[1:]
            while text[-1] in (' ','_',';'): text = text[:-1]
            if text[0]==text[-1]=="'": text = text[1:-1]
            if text[0]==text[-1]=='"': text = text[1:-1]
            return text

        #headers.DATE
        try:
            for i in range(len(self.mmcif_dict['_database_PDB_rev'][0])):
                if self.mmcif_dict['_database_PDB_rev'][0][i] == 'date_original':
                    self.headers['DATE'] = self.mmcif_dict['_database_PDB_rev'][1][i]
                    break
        except:
            if '_database_PDB_rev.date_original' in self.mmcif_dict:
                self.headers['DATE'] = self.mmcif_dict['_database_PDB_rev.date_original']
                if self.headers['DATE'] == '?':
                    self.headers['DATE'] = self.mmcif_dict['_database_PDB_rev.date']
            else:
                try:
                    for i in range(len(self.mmcif_dict['_pdbx_database_status'][0])):
                        if self.mmcif_dict['_pdbx_database_status'][0][i] == 'recvd_initial_deposition_date':
                            self.headers['DATE'] = self.mmcif_dict['_pdbx_database_status'][1][i]
                            break
                except:
                    if '_pdbx_database_status.recvd_initial_deposition_date' in self.mmcif_dict:
                        self.headers['DATE'] = self.mmcif_dict['_pdbx_database_status.recvd_initial_deposition_date']
                    if self.headers['DATE'] == '?':
                        print('DATE is unknown')

        self.allwords[self.headers['DATE']] = 1

        #headers.HEADER
        self.headers['HEADER'] = self.mmcif_dict['_struct_keywords.pdbx_keywords'].replace('&&&&&',' ').replace("'","").replace('  ',' ') \
                                 if '_struct_keywords.pdbx_keywords' in self.mmcif_dict else '?'

        for i in self.headers['HEADER'].split():
            if i not in self.allwords: self.allwords[i] = 1

        #headers.PDBFILE
        self.headers['PDBFILE'] = self.mmcif_dict['_entry.id'] if '_entry.id' in self.mmcif_dict else '?'

        self.allwords[self.headers['PDBFILE']] = 1

        #headers.TITLE
        self.headers['TITLE'] = CleanText(self.mmcif_dict['_struct.title']) if '_struct.title' in self.mmcif_dict else '?'

        #headers.CHAINBIO
        if '_pdbx_struct_assembly_gen' in self.mmcif_dict:
            #print(self.mmcif_dict['_pdbx_struct_assembly_gen'])
            for row in self.mmcif_dict['_pdbx_struct_assembly_gen'][1:]:
                for chh in row[2].split(','):
                    if chh not in self.headers['CHAINBIO']:
                        self.headers['CHAINBIO'][chh] = row[0]

        for i in self.headers['TITLE'].split():
            if i not in self.allwords: self.allwords[i] = 1

        #headers.KEYWDS
        self.headers['KEYWDS'] = [CleanText(x) for x in CleanText(self.mmcif_dict['_struct_keywords.text']).split(',') if x not in ('',"'",'"',"' ",'" ')] \
                                 if '_struct_keywords.text' in self.mmcif_dict else []

        for i in self.headers['KEYWDS']:
            for j in i.split():
                if j not in self.allwords: self.allwords[i] = 1

        #headers.EXPDTA
        if '_exptl.method' in self.mmcif_dict:
            self.headers['EXPDTA'] = self.mmcif_dict['_exptl.method'][1:-1].replace('&&&&&',' ').split(';')
        elif '_exptl' in self.mmcif_dict:
            for i in range(len(self.mmcif_dict['_exptl'][0])):
                if self.mmcif_dict['_exptl'][0][i] == 'method':
                    self.headers['EXPDTA'] = [CleanText(x[i]) for x in self.mmcif_dict['_exptl'][1:]]
                    break

        for i in self.headers['EXPDTA']:
            for j in i.split():
                if j not in self.allwords: self.allwords[i] = 1

        #headers.NUMMDL
        dot = modelname.rfind('.')
        pattern = modelname[:dot]+'.cif*'
        self.headers['NUMMDL'] = len(glob.glob(pattern))

        #headers.MDLTYP
        if '_pdbx_coordinate_model' in self.mmcif_dict:
            temp_dict = {}
            for ent in self.mmcif_dict['_pdbx_coordinate_model'][1:]:
                if ent[1][1:-1] not in temp_dict: temp_dict[ent[1][1:-1]] = []
                temp_dict[ent[1][1:-1]].append(ent[0])
            for x in sorted(list(temp_dict.keys())):
                self.headers['MDLTYP'].append(x.replace('&&&&&',' ')+', CHAIN '+', '.join(temp_dict[x]))
            del temp_dict

        elif '_pdbx_coordinate_model.type' in self.mmcif_dict:
            self.headers['MDLTYP'].append(CleanText(self.mmcif_dict['_pdbx_coordinate_model.type'])+', CHAIN '\
                                          +CleanText(self.mmcif_dict['_pdbx_coordinate_model.asym_id']))

        elif '_struct.pdbx_model_type_details' in self.mmcif_dict and \
             self.mmcif_dict['_struct.pdbx_model_type_details']!= '?':

            self.headers['MDLTYP'] += self.mmcif_dict['_struct.pdbx_model_type_details'].replace('&&&&&',' ').split(' ; ')

        if self.headers['MDLTYP'] == ['.']: self.headers['MDLTYP'] = []

        for i in self.headers['MDLTYP']:
            for j in i.split():
                if j not in self.allwords: self.allwords[i] = 1

        #headers.AUTHOR
        if '_audit_author.name' in self.mmcif_dict:
            self.headers['AUTHOR'] = ''.join(self.mmcif_dict['_audit_author.name'][1:-1].split(',&&&&&')[::-1]).replace('&&&&&',' ')

        elif '_audit_author' in self.mmcif_dict:
            if   self.mmcif_dict['_audit_author'][0][0] == 'name': flag = 0
            elif self.mmcif_dict['_audit_author'][0][1] == 'name': flag = 1

            for k in range(len(self.mmcif_dict['_audit_author'][1:])):

                xx = self.mmcif_dict['_audit_author'][k+1][flag] # remove terminal quotes
                if xx[0]==xx[-1] and xx[0] in ("'",'"'):
                    self.mmcif_dict['_audit_author'][k+1][flag] = self.mmcif_dict['_audit_author'][k+1][flag][1:-1]

                if self.mmcif_dict['_audit_author'][k+1][flag][-1]==',': # fix last symbol if it's comma
                    self.mmcif_dict['_audit_author'][k+1][flag] = self.mmcif_dict['_audit_author'][k+1][flag][:-1]+'.'

                self.mmcif_dict['_audit_author'][k+1][flag] = self.mmcif_dict['_audit_author'][k+1][flag].replace(', ',',&&&&&').replace(',',',&&&&&').replace('&&&&&&&&&&','&&&&&')

                # fix for missing dot
                if ',&&&&&' in self.mmcif_dict['_audit_author'][k+1][flag] and self.mmcif_dict['_audit_author'][k+1][flag][-1]!='.':
                    self.mmcif_dict['_audit_author'][k+1][flag] += '.'

            self.headers['AUTHOR'] = ','.join([''.join(x[flag].split(',&&&&&')[::-1])
                                               for x in self.mmcif_dict['_audit_author'][1:]]).replace('&&&&&',' ')

        self.headers['AUTHOR'] = self.headers['AUTHOR'].replace(' ,',',')

        for i in self.headers['AUTHOR'].split(','):
                if i not in self.allwords: self.allwords[i] = 1

        #headers.RESOL
        if '_refine.ls_d_res_high' in self.mmcif_dict and self.mmcif_dict['_refine.ls_d_res_high'] not in ('.','?'):
            self.headers['RESOL'] = float(self.mmcif_dict['_refine.ls_d_res_high'])
        elif '_em_3d_reconstruction.resolution' in self.mmcif_dict and self.mmcif_dict['_em_3d_reconstruction.resolution'] not in ('.','?'):
            self.headers['RESOL'] = float(self.mmcif_dict['_em_3d_reconstruction.resolution'])
        elif '_reflns_shell.d_res_high' in self.mmcif_dict and self.mmcif_dict['_reflns_shell.d_res_high'] not in ('.','?'):
            self.headers['RESOL'] = float(self.mmcif_dict['_reflns_shell.d_res_high'])

        if type(self.headers['RESOL'])==float and str(self.headers['RESOL']) not in self.allwords: self.allwords[str(self.headers['RESOL'])] = 1

        #headers.CRYSYM - ommited

        #headers.AREMODELS
        if self.headers['NUMMDL'] > 1:
            self.headers['AREMODELS'] = True

        #headers.MDLNO
        if modelname[-1] in '0123456789':
            self.headers['MDLNO'] = int(modelname[dot+4:])
        del dot,pattern

        #headers.MODRES
        if '_pdbx_struct_mod_residue.id' in self.mmcif_dict:
            name1 = self.mmcif_dict['_pdbx_struct_mod_residue.label_comp_id']
            name2 = self.mmcif_dict['_pdbx_struct_mod_residue.parent_comp_id']
            self.headers['MODRES'][name1] = name2
            if name1 not in self.allwords: self.allwords[name1] = 1
            if CleanText(self.mmcif_dict['_pdbx_struct_mod_residue.details']) not in self.allwords:
                self.allwords[CleanText(self.mmcif_dict['_pdbx_struct_mod_residue.details'])] = 1
            del name1,name2
        elif '_pdbx_struct_mod_residue' in self.mmcif_dict:
            ind1,ind2,ind3 = 3,8,9
            for i in range(len(self.mmcif_dict['_pdbx_struct_mod_residue'][0])):
                if   self.mmcif_dict['_pdbx_struct_mod_residue'][0][i] == 'label_comp_id' : ind1 = i
                elif self.mmcif_dict['_pdbx_struct_mod_residue'][0][i] == 'parent_comp_id': ind2 = i
                elif self.mmcif_dict['_pdbx_struct_mod_residue'][0][i] == 'details'       : ind3 = i
            for ent in self.mmcif_dict['_pdbx_struct_mod_residue'][1:]:
                self.headers['MODRES'][ent[ind1]] = ent[ind2]
                if ent[ind1] not in self.allwords: self.allwords[ent[ind1]] = 1
                if CleanText(ent[ind3]) not in self.allwords: self.allwords[CleanText(ent[ind3])] = 1

        for ent in self.headers['MODRES']:
            if self.headers['MODRES'][ent] in ('.','?'): self.headers['MODRES'][ent] = ''

        #headers.TYPE
        Types = {'polyribonucleotide'     :'RNA', "'polypeptide(L)'":'Protein', '"polypeptide(L)"':'Protein',
                 'polydeoxyribonucleotide':'DNA', "'polypeptide(D)'":'Protein', '"polypeptide(D)"':'Protein',
                 'other' : 'RNA', "'peptide&&&&&nucleic&&&&&acid'": 'PNA',
                 "'polydeoxyribonucleotide/polyribonucleotide&&&&&hybrid'":'RNA'} # Are you sure?

        r,p,d,x,y = '','','',1,4
        
        if '_entity_poly' in self.mmcif_dict:
            for i in range(len(self.mmcif_dict['_entity_poly'][0])):
                if   self.mmcif_dict['_entity_poly'][0][i]=='type': x = i
                elif self.mmcif_dict['_entity_poly'][0][i]=='pdbx_seq_one_letter_code': y = i

            for ent in self.mmcif_dict['_entity_poly'][1:]:
                if ent[x] not in self.allwords: self.allwords[ent[x]] = 1
                if ent[y] not in self.allwords: self.allwords[ent[y]] = 1
                if   Types[ent[x]]=='RNA'    : r = 'R'
                elif Types[ent[x]] in ('Protein','PNA'): p = 'P'
                elif Types[ent[x]]=='DNA'    : d = 'D'

        elif '_entity_poly.type' in self.mmcif_dict:
            typ = Types[self.mmcif_dict['_entity_poly.type']]
            if self.mmcif_dict['_entity_poly.type'] not in self.allwords: self.allwords[self.mmcif_dict['_entity_poly.type']] = 1
            if self.mmcif_dict['_entity_poly.pdbx_seq_one_letter_code'] not in self.allwords:
                self.allwords[self.mmcif_dict['_entity_poly.pdbx_seq_one_letter_code']] = 1
            if   typ=='RNA'    : r = 'R'
            elif typ in ('Protein','PNA'): p = 'P'
            elif typ=='DNA'    : d = 'D'
            del typ

        self.headers['TYPE'] = r+p+d
        del r,p,d,x

        # initialize MOLECULES and CHAINS
        if '_entity_poly' not in self.mmcif_dict and '_entity_poly.entity_id' in self.mmcif_dict:

            mol = int(self.mmcif_dict['_entity_poly.entity_id'])
            self.molecules[mol] = Tools.Molecule(mol)
            chains = self.mmcif_dict['_entity_poly.pdbx_strand_id']
            self.molecules[mol]['CHAIN'] = chains
            for ch in chains.split(','):
                self.chains[ch] = Tools.Chain(ch)
                self.chains[ch]['MOL_ID'] = mol
                self.chains[ch]['SEQ'] = self.mmcif_dict['_entity_poly.pdbx_seq_one_letter_code']
                self.chains[ch]['TYPE'] = Types[self.mmcif_dict['_entity_poly.type']]
            del mol,chains

        elif '_entity_poly' in self.mmcif_dict:

            tt = {'pdbx_strand_id':6,'type':1,'pdbx_seq_one_letter_code':4,'entity_id':0} # temp_token_indexes

            for i in range(len(self.mmcif_dict['_entity_poly'][0])):

                if self.mmcif_dict['_entity_poly'][0][i] in tt: tt[self.mmcif_dict['_entity_poly'][0][i]] = i

            for mol in self.mmcif_dict['_entity_poly'][1:]:

                self.molecules[int(mol[tt['entity_id']])] = Tools.Molecule(int(mol[tt['entity_id']]))
                chains = mol[tt['pdbx_strand_id']]
                self.molecules[int(mol[tt['entity_id']])]['CHAIN'] = chains
                for ch in chains.split(','):
                    self.chains[ch] = Tools.Chain(ch)
                    self.chains[ch]['MOL_ID'] = int(mol[tt['entity_id']])
                    self.chains[ch]['SEQ'] = mol[tt['pdbx_seq_one_letter_code']]
                    self.chains[ch]['TYPE'] = Types[mol[tt['type']]]
            del chains

        #fill analogs of COMPND and SOURCE
        # in case of |molecules|>1
        for loop in ('_entity','_entity_name_com','_pdbx_entity_src_syn','_entity_src_nat','_entity_src_gen'):
            if loop not in self.mmcif_dict: continue

            tokens = {}

            for i in range(len(self.mmcif_dict[loop][0])):
                if self.mmcif_dict[loop][0][i] in ('id','entity_id'): indOfID = i
                if self.mmcif_dict[loop][0][i] in Tools.CIFTOPDB_Mol[loop]:
                    tokens[i] = Tools.CIFTOPDB_Mol[loop][self.mmcif_dict[loop][0][i]]

            for ent in self.mmcif_dict[loop][1:]:

                IDent = int(ent[indOfID])
                if IDent not in self.molecules: continue

                for j in tokens:
                    if ent[j] != '?':
                        if tokens[j] == 'SYNTHETIC':
                            if ent[j]=='syn':
                                self.molecules[IDent]['SYNTHETIC']  = 'YES'
                                self.molecules[IDent]['ENGINEERED'] = 'YES'
                            elif ent[j]=='man':
                                self.molecules[IDent]['ENGINEERED'] = 'YES'
                        else:
                            self.molecules[IDent][tokens[j]] = CleanText(ent[j])

        # in case of |molecules|==1
        if '_pdbx_entity_src_syn.entity_id' in self.mmcif_dict:
            ent = int(self.mmcif_dict['_pdbx_entity_src_syn.entity_id'])
            if ent in self.molecules:
                if self.mmcif_dict['_pdbx_entity_src_syn.organism_scientific'] != '?':
                    self.molecules[ent]['ORGANISM_SCIENTIFIC'] = CleanText(self.mmcif_dict['_pdbx_entity_src_syn.organism_scientific'])
                if self.mmcif_dict['_pdbx_entity_src_syn.organism_common_name'] != '?':
                    self.molecules[ent]['ORGANISM_COMMON'] = CleanText(self.mmcif_dict['_pdbx_entity_src_syn.organism_common_name'])
                if self.mmcif_dict['_pdbx_entity_src_syn.ncbi_taxonomy_id'] != '?':
                    self.molecules[ent]['ORGANISM_TAXID'] = CleanText(self.mmcif_dict['_pdbx_entity_src_syn.ncbi_taxonomy_id'])
                if self.mmcif_dict['_pdbx_entity_src_syn.details'] != '?':
                    self.molecules[ent]['OTHER_DETAILS'] = CleanText(self.mmcif_dict['_pdbx_entity_src_syn.details'])

        if '_entity_name_com.entity_id' in self.mmcif_dict:
            ent = int(self.mmcif_dict['_entity_name_com.entity_id'])
            if ent in self.molecules:
                self.molecules[ent]['SYNONYM'] = CleanText(self.mmcif_dict['_entity_name_com.name'])

        if '_entity.id' in self.mmcif_dict:
            ent = int(self.mmcif_dict['_entity.id'])
            if ent in self.molecules:
                if self.mmcif_dict['_entity.pdbx_description'] != '?':
                    self.molecules[ent]['MOLECULE'] = CleanText(self.mmcif_dict['_entity.pdbx_description'])
                if '_entity.pdbx_fragment' in self.mmcif_dict and self.mmcif_dict['_entity.pdbx_fragment'] != '?':
                    self.molecules[ent]['FRAGMENT'] = CleanText(self.mmcif_dict['_entity.pdbx_fragment'])
                if '_entity.pdbx_ec' in self.mmcif_dict and self.mmcif_dict['_entity.pdbx_ec'] != '?':
                    self.molecules[ent]['EC'] = CleanText(self.mmcif_dict['_entity.pdbx_ec'])
                if '_entity.pdbx_mutation' in self.mmcif_dict and self.mmcif_dict['_entity.pdbx_mutation'] != '?':
                    self.molecules[ent]['MUTATION'] = CleanText(self.mmcif_dict['_entity.pdbx_mutation'])
                if self.mmcif_dict['_entity.details'] != '?':
                    self.molecules[ent]['OTHER_DETAILS'] = CleanText(self.mmcif_dict['_entity.details'])
                if self.mmcif_dict['_entity.src_method'] != 'nat':
                    if self.mmcif_dict['_entity.src_method'] == 'syn':
                        self.molecules[ent]['ENGINEERED'] = 'YES'
                        self.molecules[ent]['SYNTHETIC']  = 'YES'
                    elif self.mmcif_dict['_entity.src_method'] == 'man':
                        self.molecules[ent]['ENGINEERED'] = 'YES'

        if '_entity_src_nat.entity_id' in self.mmcif_dict:
            ent = int(self.mmcif_dict['_entity_src_nat.entity_id'])
            if ent in self.molecules:
                if self.mmcif_dict['_entity_src_nat.pdbx_fragment'] != '?':
                    self.molecules[ent]['FRAGMENT'] = CleanText(self.mmcif_dict['_entity_src_nat.pdbx_fragment'])
                if self.mmcif_dict['_entity_src_nat.pdbx_organism_scientific'] != '?':
                    self.molecules[ent]['ORGANISM_SCIENTIFIC'] = CleanText(self.mmcif_dict['_entity_src_nat.pdbx_organism_scientific'])
                if self.mmcif_dict['_entity_src_nat.common_name'] != '?':
                    self.molecules[ent]['ORGANISM_COMMON'] = CleanText(self.mmcif_dict['_entity_src_nat.common_name'])
                if self.mmcif_dict['_entity_src_nat.pdbx_ncbi_taxonomy_id'] != '?':
                    self.molecules[ent]['ORGANISM_TAXID'] = CleanText(self.mmcif_dict['_entity_src_nat.pdbx_ncbi_taxonomy_id'])
                if self.mmcif_dict['_entity_src_nat.strain'] != '?':
                    self.molecules[ent]['STRAIN'] = CleanText(self.mmcif_dict['_entity_src_nat.strain'])
                if self.mmcif_dict['_entity_src_nat.pdbx_variant'] != '?':
                    self.molecules[ent]['VARIANT'] = CleanText(self.mmcif_dict['_entity_src_nat.pdbx_variant'])
                if self.mmcif_dict['_entity_src_nat.pdbx_cell_line'] != '?':
                    self.molecules[ent]['CELL_LINE'] = CleanText(self.mmcif_dict['_entity_src_nat.pdbx_cell_line'])
                if self.mmcif_dict['_entity_src_nat.pdbx_atcc'] != '?':
                    self.molecules[ent]['ATCC'] = CleanText(self.mmcif_dict['_entity_src_nat.pdbx_atcc'])
                if self.mmcif_dict['_entity_src_nat.pdbx_organ'] != '?':
                    self.molecules[ent]['ORGAN'] = CleanText(self.mmcif_dict['_entity_src_nat.pdbx_organ'])
                if self.mmcif_dict['_entity_src_nat.tissue'] != '?':
                    self.molecules[ent]['TISSUE'] = CleanText(self.mmcif_dict['_entity_src_nat.tissue'])
                if self.mmcif_dict['_entity_src_nat.pdbx_cell'] != '?':
                    self.molecules[ent]['CELL'] = CleanText(self.mmcif_dict['_entity_src_nat.pdbx_cell'])
                if self.mmcif_dict['_entity_src_nat.pdbx_organelle'] != '?':
                    self.molecules[ent]['ORGANELLE'] = CleanText(self.mmcif_dict['_entity_src_nat.pdbx_organelle'])
                if self.mmcif_dict['_entity_src_nat.pdbx_secretion'] != '?':
                    self.molecules[ent]['SECRETION'] = CleanText(self.mmcif_dict['_entity_src_nat.pdbx_secretion'])
                if self.mmcif_dict['_entity_src_nat.pdbx_cellular_location'] != '?':
                    self.molecules[ent]['CELLULAR_LOCATION'] = CleanText(self.mmcif_dict['_entity_src_nat.pdbx_cellular_location'])
                if self.mmcif_dict['_entity_src_nat.pdbx_plasmid_name'] != '?':
                    self.molecules[ent]['PLASMID'] = CleanText(self.mmcif_dict['_entity_src_nat.pdbx_plasmid_name'])
                if self.mmcif_dict['_entity_src_nat.details'] != '?':
                    self.molecules[ent]['OTHER_DETAILS'] = CleanText(self.mmcif_dict['_entity_src_nat.details'])

        if '_entity_src_gen.entity_id' in self.mmcif_dict:
            ent = int(self.mmcif_dict['_entity_src_gen.entity_id'])
            if ent in self.molecules:
                if self.mmcif_dict['_entity_src_gen.pdbx_gene_src_fragment'] != '?':
                    self.molecules[ent]['FRAGMENT'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_gene_src_fragment'])
                if self.mmcif_dict['_entity_src_gen.pdbx_gene_src_scientific_name'] != '?':
                    self.molecules[ent]['ORGANISM_SCIENTIFIC'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_gene_src_scientific_name'])
                if self.mmcif_dict['_entity_src_gen.gene_src_common_name'] != '?':
                    self.molecules[ent]['ORGANISM_COMMON'] = CleanText(self.mmcif_dict['_entity_src_gen.gene_src_common_name'])
                if self.mmcif_dict['_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id'] != '?':
                    self.molecules[ent]['ORGANISM_TAXID'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id'])
                if self.mmcif_dict['_entity_src_gen.gene_src_strain'] != '?':
                    self.molecules[ent]['STRAIN'] = CleanText(self.mmcif_dict['_entity_src_gen.gene_src_strain'])
                if self.mmcif_dict['_entity_src_gen.pdbx_gene_src_variant'] != '?':
                    self.molecules[ent]['VARIANT'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_gene_src_variant'])
                if self.mmcif_dict['_entity_src_gen.pdbx_gene_src_cell_line'] != '?':
                    self.molecules[ent]['CELL_LINE'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_gene_src_cell_line'])
                if self.mmcif_dict['_entity_src_gen.pdbx_gene_src_atcc'] != '?':
                    self.molecules[ent]['ATCC'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_gene_src_atcc'])
                if self.mmcif_dict['_entity_src_gen.pdbx_gene_src_organ'] != '?':
                    self.molecules[ent]['ORGAN'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_gene_src_organ'])
                if self.mmcif_dict['_entity_src_gen.gene_src_tissue'] != '?':
                    self.molecules[ent]['TISSUE'] = CleanText(self.mmcif_dict['_entity_src_gen.gene_src_tissue'])
                if self.mmcif_dict['_entity_src_gen.pdbx_gene_src_cell'] != '?':
                    self.molecules[ent]['CELL'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_gene_src_cell'])
                if self.mmcif_dict['_entity_src_gen.pdbx_gene_src_organelle'] != '?':
                    self.molecules[ent]['ORGANELLE'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_gene_src_organelle'])
                if self.mmcif_dict['_entity_src_gen.pdbx_gene_src_cellular_location'] != '?':
                    self.molecules[ent]['CELLULAR_LOCATION'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_gene_src_cellular_location'])
                if self.mmcif_dict['_entity_src_gen.pdbx_gene_src_gene'] != '?':
                    self.molecules[ent]['GENE'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_gene_src_gene'])
                if self.mmcif_dict['_entity_src_gen.pdbx_host_org_scientific_name'] != '?':
                    self.molecules[ent]['EXPRESSION_SYSTEM'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_host_org_scientific_name'])
                if self.mmcif_dict['_entity_src_gen.host_org_common_name'] != '?':
                    self.molecules[ent]['EXPRESSION_SYSTEM_COMMON'] = CleanText(self.mmcif_dict['_entity_src_gen.host_org_common_name'])
                if self.mmcif_dict['_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id'] != '?':
                    self.molecules[ent]['EXPRESSION_SYSTEM_TAXID'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id'])
                if self.mmcif_dict['_entity_src_gen.pdbx_host_org_strain'] != '?':
                    self.molecules[ent]['EXPRESSION_SYSTEM_STRAIN'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_host_org_strain'])
                if self.mmcif_dict['_entity_src_gen.pdbx_host_org_variant'] != '?':
                    self.molecules[ent]['EXPRESSION_SYSTEM_VARIANT'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_host_org_variant'])
                if self.mmcif_dict['_entity_src_gen.pdbx_host_org_cell_line'] != '?':
                    self.molecules[ent]['EXPRESSION_SYSTEM_CELL_LINE'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_host_org_cell_line'])
                if self.mmcif_dict['_entity_src_gen.pdbx_host_org_atcc'] != '?':
                    self.molecules[ent]['EXPRESSION_SYSTEM_ATCC_NUMBER'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_host_org_atcc'])
                if self.mmcif_dict['_entity_src_gen.pdbx_host_org_organ'] != '?':
                    self.molecules[ent]['EXPRESSION_SYSTEM_ORGAN'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_host_org_organ'])
                if self.mmcif_dict['_entity_src_gen.pdbx_host_org_tissue'] != '?':
                    self.molecules[ent]['EXPRESSION_SYSTEM_TISSUE'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_host_org_tissue'])
                if self.mmcif_dict['_entity_src_gen.pdbx_host_org_cell'] != '?':
                    self.molecules[ent]['EXPRESSION_SYSTEM_CELL'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_host_org_cell'])
                if self.mmcif_dict['_entity_src_gen.pdbx_host_org_organelle'] != '?':
                    self.molecules[ent]['EXPRESSION_SYSTEM_ORGANELLE'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_host_org_organelle'])
                if self.mmcif_dict['_entity_src_gen.pdbx_host_org_cellular_location'] != '?':
                    self.molecules[ent]['EXPRESSION_SYSTEM_CELLULAR_LOCATION'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_host_org_cellular_location'])
                if self.mmcif_dict['_entity_src_gen.pdbx_host_org_vector_type'] != '?':
                    self.molecules[ent]['EXPRESSION_SYSTEM_VECTOR_TYPE'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_host_org_vector_type'])
                if self.mmcif_dict['_entity_src_gen.pdbx_host_org_vector'] != '?':
                    self.molecules[ent]['EXPRESSION_SYSTEM_VECTOR'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_host_org_vector'])
                if self.mmcif_dict['_entity_src_gen.plasmid_name'] != '?':
                    self.molecules[ent]['EXPRESSION_SYSTEM_PLASMID'] = CleanText(self.mmcif_dict['_entity_src_gen.plasmid_name'])
                if self.mmcif_dict['_entity_src_gen.pdbx_host_org_gene'] != '?':
                    self.molecules[ent]['EXPRESSION_SYSTEM_GENE'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_host_org_gene'])
                if self.mmcif_dict['_entity_src_gen.pdbx_description'] != '?':
                    self.molecules[ent]['OTHER_DETAILS'] = CleanText(self.mmcif_dict['_entity_src_gen.pdbx_description'])

        # resolve dependencies
        for i in self.molecules:
            self.molecules[i]['OTHER_DETAILS'] = self.molecules[i]['OTHER_DETAILS'].replace('linkerC-terminal','linker; C-terminal')
            if self.molecules[i]['EXPRESSION_SYSTEM_VECTOR_TYPE']:
                self.molecules[i]['ENGINEERED'] = 'YES'
            if self.molecules[i]['ORGANISM_SCIENTIFIC'].upper() == 'SYNTHETIC CONSTRUCT' or\
               self.molecules[i]['OTHER_DETAILS'].upper() == 'RNA SYNTHESIS':
                self.molecules[i]['ENGINEERED'] = 'YES'
                self.molecules[i]['SYNTHETIC']  = 'YES'

        for i in self.molecules:
            for k,v in self.molecules[i].items():
                if str(v).upper()=='YES' and k not in self.allwords: self.allwords[k] = 1
                else:
                    for w in str(v).split():
                        if w not in self.allwords: self.allwords[w] = 1

        # Residue Types (self.residues) - mmcif._chem_comp  ; liglist + metlist
        liglist,metlist = [],[]

        if '_chem_comp' in self.mmcif_dict:

            ind0 = self.mmcif_dict['_chem_comp'][0].index('id')
            ind1 = self.mmcif_dict['_chem_comp'][0].index('type')
            try:
                ind3 = self.mmcif_dict['_chem_comp'][0].index('name')
            except:
                ind3 = ind0
            

            for ent in self.mmcif_dict['_chem_comp'][1:]:

                for j in ent:
                    for jj in CleanText(j).split():
                        if jj not in self.allwords: self.allwords[jj] = 1

                TType = Tools.CIFrestype([CleanText(ent[ind1]),CleanText(ent[ind3])])
                self.restype[ent[ind0]] = TType
                if   TType=='Unknown': liglist.append(ent[ind0])
                elif TType=='Metal'  : metlist.append(ent[ind0])

            if liglist: self.headers['LIGLIST'] = ','+','.join(liglist)+','
            if metlist: self.headers['METLIST'] = ','+','.join(metlist)+','
            del liglist,metlist,TType,ind0,ind1,ind3

        self.allwords = ' '.join([str(x) for x in list(self.allwords.keys())])

        # Residues to Chains (chains.RES) ; _pdbx_poly_seq_scheme ; Type=Unknown if Type(chain)!=Type(residue) or CIFID == '.'  !!
        cif_resind = {}

        if '_pdbx_poly_seq_scheme' in self.mmcif_dict:

            inds = [0]*len(self.mmcif_dict['_pdbx_poly_seq_scheme'][0])
            inds[2]  = self.mmcif_dict['_pdbx_poly_seq_scheme'][0].index('seq_id')
            inds[3]  = self.mmcif_dict['_pdbx_poly_seq_scheme'][0].index('mon_id')
            inds[5]  = self.mmcif_dict['_pdbx_poly_seq_scheme'][0].index('pdb_seq_num')
            inds[9]  = self.mmcif_dict['_pdbx_poly_seq_scheme'][0].index('pdb_strand_id')
            inds[10] = self.mmcif_dict['_pdbx_poly_seq_scheme'][0].index('pdb_ins_code')

            for ent in self.mmcif_dict['_pdbx_poly_seq_scheme'][1:]:
                res = Tools.ResidueCIF(ent,inds)
                res['TYPE'] = self.restype[res['NAME']]
                self.chains[res['CHAIN']]['RES'].append(res)
                cif_resind[(res['CHAIN'],res['CIFID'])] = len(self.chains[res['CHAIN']]['RES'])-1
            del inds

        # ATOM+HETATM (_atom_site)
        pdb_ligind = {}

        longestchain = ''
        maxlen = 0

        if self.chains:

            for ch in self.chains:
                if len(self.chains[ch]['SEQ'])>maxlen: maxlen,longestchain = len(self.chains[ch]['SEQ']), ch

            del ch,maxlen

        if '_atom_site' in self.mmcif_dict:

            for ent in self.mmcif_dict['_atom_site'][1:]:

                if ent[2] == 'H': continue # We do not consider Hydrogen atoms

                atom = Tools.AtomCIF(ent,self.mmcif_dict['_atom_site'][0])
                atom['TYPE'] = self.restype[atom['RESNAME']] if atom['RESNAME'] in self.restype \
                               else Tools.restype[atom['RESNAME']] if atom['RESNAME'] in Tools.restype else 'Unknown' 

                if type(atom['CIFID']) == int:

                    if atom['CHAIN'] not in self.chains:
                        self.chains[atom['CHAIN']] = Tools.Chain(atom['CHAIN'])

                    if (atom['CHAIN'],atom['CIFID']) not in cif_resind:
                        res = Tools.ResidueCIFAtom(atom)
                        self.chains[atom['CHAIN']]['RES'].append(res)
                        cif_resind[(atom['CHAIN'],atom['CIFID'])] = len(self.chains[atom['CHAIN']]['RES'])-1
                    
                    self.chains[atom['CHAIN']]['RES'][cif_resind[(atom['CHAIN'],atom['CIFID'])]]['ATOMS'].append(atom)
                    self.chains[atom['CHAIN']]['RES'][cif_resind[(atom['CHAIN'],atom['CIFID'])]]['CIFCHAIN'] = atom['CIFCHAIN']
                else:
                    if (atom['CHAIN'],atom['RESNUM']) not in pdb_ligind:
                        if atom['CHAIN'] not in self.chains:
                            self.headers['MASKEDCHS'][atom['CHAIN']] = longestchain
                            atom['CHAIN'] = longestchain
                        res = Tools.ResidueCIFAtom(atom)
                        res['ATOMS'].append(atom)
                        self.chains[atom['CHAIN']]['LIGANDS'].append(res)
                        pdb_ligind[(atom['CHAIN'],atom['RESNUM'])] = len(self.chains[atom['CHAIN']]['LIGANDS'])-1
                    else:
                        self.chains[atom['CHAIN']]['LIGANDS'][pdb_ligind[(atom['CHAIN'],atom['RESNUM'])]]['ATOMS'].append(atom)
                        self.chains[atom['CHAIN']]['LIGANDS'][pdb_ligind[(atom['CHAIN'],atom['RESNUM'])]]['CIFCHAIN'] = atom['CIFCHAIN']
            del cif_resind,pdb_ligind,atom,res

        # START,END,SEQ,LENGTH,LIGSEQ,SEQ2 for chains + MISS=True for residues with no atoms
        for ch in sorted(list(self.chains.keys())):

            for i in range(len(self.chains[ch]['RES'])):
                self.chains[ch]['NUMS'][self.chains[ch]['RES'][i]['PDBNUM']] = i
                if not self.chains[ch]['RES'][i]['ATOMS']:
                    self.chains[ch]['RES'][i]['MISS'] = True
                    self.chains[ch]['RES'][i]['BRACKETS'] = '-'
                    self.chains[ch]['RES'][i]['SLBRACKETS'] = '-'

            for i in range(len(self.chains[ch]['LIGANDS'])):
                if not self.chains[ch]['LIGANDS'][i]['ATOMS']:
                    self.chains[ch]['LIGANDS'][i]['MISS'] = True
                    self.chains[ch]['LIGANDS'][i]['BRACKETS'] = '-'
                    self.chains[ch]['LIGANDS'][i]['SLBRACKETS'] = '-'

            self.chains[ch]['START']   = int(self.chains[ch]['RES'][0]['FLOAT'])
            self.chains[ch]['END']     = int(self.chains[ch]['RES'][-1]['FLOAT'])
            self.chains[ch]['LENGTH']  = len(self.chains[ch]['RES'])
            self.chains[ch]['SEQ2']    = [r['NAME'] for r in self.chains[ch]['RES']]
            self.chains[ch]['SEQ']     = self.chains[ch]['SEQ2']
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

'''
def Check(pdb,cif):

    if sorted(list(pdb.headers.keys())) != sorted(list(cif.headers.keys())):
        print(pdb.headers['PDBFILE'],pdb.headers['MDLNO'],'Wrong set of headers.keys()')
        return
    for key in ('HEADER','MODRES','AUTHOR'):
        if key=='CRYSYM': continue
        if pdb.headers[key] != cif.headers[key]:
            if type(pdb.headers[key]) == str and pdb.headers[key] == cif.headers[key].upper(): continue
            if type(pdb.headers[key]) == str and pdb.headers[key].replace(', ',',') == cif.headers[key].upper(): continue
            if key=='TITLE' and pdb.headers[key].replace(' ','') == cif.headers[key].upper().replace(' ',''): continue
            if key=='KEYWDS':
                c = 0
                for word in pdb.headers[key]:
                    if word.replace('/','-') not in [x.upper() for x in cif.headers[key]]:
                        c+=1
                if c <3: continue
                if any([',' in x for x in pdb.headers[key]]): continue
            if key=='RESOL':
                if type(pdb.headers[key])==float and type(cif.headers[key])==float:
                    if abs(pdb.headers[key]-cif.headers[key])<0.01: continue
            if key=='AUTHOR':
                if pdb.headers[key].replace(' ','.') == cif.headers[key].upper(): continue
                print(pdb.headers['PDBFILE'],pdb.headers['MDLNO'],'Wrong headers.key="%s"'%key)
                print('CIF:',cif.headers[key])
                print('PDB:',pdb.headers[key])
                continue
            if key=='MDLTYP' and len(pdb.headers[key]) == len(cif.headers[key]): continue
            print(pdb.headers['PDBFILE'],pdb.headers['MDLNO'],'Wrong headers.key="%s"'%key)
            print('CIF:',cif.headers[key])
            print('PDB:',pdb.headers[key])
            continue
    if sorted(list(pdb.molecules.keys())) != sorted(list(cif.molecules.keys())):
        print(pdb.headers['PDBFILE'],pdb.headers['MDLNO'],'Wrong set of molecules.keys()')
        return
    for key in pdb.molecules:
        if pdb.molecules[key] != cif.molecules[key]:
            for key2 in pdb.molecules[key]:
                if key2 == 'EXTENSION': continue
                try:
                    if pdb.molecules[key][key2] != cif.molecules[key][key2]:
                        if type(pdb.molecules[key][key2]) == str and pdb.molecules[key][key2].upper().replace(' ','') == cif.molecules[key][key2].upper().replace(' ',''): continue
                        if type(pdb.molecules[key][key2]) == str and pdb.molecules[key][key2][:-1].upper().replace(' ','') == cif.molecules[key][key2].upper().replace(' ',''): continue
                        if key2=='MOLECULE' and pdb.molecules[key][key2].replace('&&&&&',' ') == cif.molecules[key][key2].replace('GP *GP','GP*GP'): continue
                        if key2=='MOLECULE' and pdb.molecules[key][key2] == cif.molecules[key][key2][:-1]: continue
                        if key2=='CHAIN' and pdb.molecules[key][key2].replace(' ','') == cif.molecules[key][key2]: continue
                        if key2=='CHAIN' and pdb.molecules[key][key2].replace(',','') == cif.molecules[key][key2]: continue
                        if key2=='MUTATION' and pdb.molecules[key][key2]=='YES' and cif.molecules[key][key2]: continue
                        if key2=='SYNONYM' and pdb.molecules[key][key2].replace('INGPROT','ING PROT').replace('TICTRAN','TIC TRAN') == cif.molecules[key][key2].upper(): continue
                        if key2=='ENGINEERED' or key2=='SYNTHETIC':
                            print(pdb.headers['PDBFILE'],pdb.headers['MDLNO'],'Wrong molecule="%s", key=%s'%(key,key2))
                            print('CIF:',cif.molecules[key][key2])
                            print('PDB:',pdb.molecules[key][key2])
                            continue
                        if key2=='OTHER_DETAILS' and pdb.molecules[key][key2].replace('ON:LIN','ON: LIN').replace('TAG C','TAGC').replace('KER C','KERC').replace('ON:STR','ON: STR') == cif.molecules[key][key2].upper(): continue
                        if key2 in ('OTHER_DETAILS','EXPRESSION_SYSTEM','SYNONYM','FRAGMENT','GENE'):
                            print(pdb.headers['PDBFILE'],pdb.headers['MDLNO'],'Wrong molecule="%s", key=%s'%(key,key2))
                            print('CIF:',cif.molecules[key][key2])
                            print('PDB:',pdb.molecules[key][key2])
                            continue
                        print(pdb.headers['PDBFILE'],pdb.headers['MDLNO'],'Wrong molecule="%s", key=%s'%(key,key2))
                        print('CIF:',cif.molecules[key][key2])
                        print('PDB:',pdb.molecules[key][key2])
                        continue
                except:
                    print(pdb.headers['PDBFILE'],pdb.headers['MDLNO'], 'wtf???',key,"'",key2,"'")
    if sorted(list(pdb.chains.keys())) != sorted(list(cif.chains.keys())):
        print(pdb.headers['PDBFILE'],pdb.headers['MDLNO'],'Wrong set of chains.keys()')
        print('PDB:',list(pdb.chains.keys()))
        print('CIF:',list(cif.chains.keys()))
    for key in pdb.chains:
        if key not in cif.chains: continue
        for key2 in pdb.chains[key]:
            if key2=='RES':
                for i in range(min(len(pdb.chains[key]['RES']),len(cif.chains[key]['RES']))):
                    if {x:pdb.chains[key]['RES'][i][x] for x in pdb.chains[key]['RES'][i] if x not in ('CIFID','ATOMS','TYPE','ID')}  != {x:cif.chains[key]['RES'][i][x] for x in cif.chains[key]['RES'][i] if x not in ('CIFID','ATOMS','TYPE','ID')} and\
                       len(pdb.chains[key]['RES'][i]['ATOMS'])==len(cif.chains[key]['RES'][i]['ATOMS']):
                        diff = [x for x in pdb.chains[key]['RES'][i] if pdb.chains[key]['RES'][i][x]!=cif.chains[key]['RES'][i][x]]
                        if not diff: continue
                        if diff == ['ATOMS']:
                            flag = False
                            for j in range(len(pdb.chains[key]['RES'][i]['ATOMS'])):
                                diff2 = [x for x in pdb.chains[key]['RES'][i]['ATOMS'][j]\
                                         if pdb.chains[key]['RES'][i]['ATOMS'][j][x]!=cif.chains[key]['RES'][i]['ATOMS'][j][x]]
                                if diff2 and diff2 != ['NUM']:
                                    flag = True
                            if not flag: continue
                        if sorted([x['NAME'] for x in pdb.chains[key]['RES'][i]['ATOMS']]) == sorted([x['NAME'] for x in cif.chains[key]['RES'][i]['ATOMS']]):
                            continue
                        if diff in (['ATOMS','TYPE'],['TYPE','ATOMS']):
                            print('PDB:',pdb.chains[key]['RES'][i]['TYPE'],'CIF:',cif.chains[key]['RES'][i]['TYPE'])
                            continue
                        print(pdb.headers['PDBFILE'],pdb.headers['MDLNO'],'Wrong residue %s in chain %s (keys %s)'%(i,key,','.join(diff)))
                        continue
            elif key2=='LIGANDS':
                for i in range(min(len(pdb.chains[key]['LIGANDS']),len(cif.chains[key]['LIGANDS']))):
                    if {x:pdb.chains[key]['LIGANDS'][i][x] for x in pdb.chains[key]['LIGANDS'][i] if x not in ('CIFID','ATOMS','TYPE','ID')} != {x:cif.chains[key]['LIGANDS'][i][x] for x in cif.chains[key]['LIGANDS'][i] if x not in ('CIFID','ATOMS','TYPE','ID')} and\
                       len(pdb.chains[key]['LIGANDS'][i]['ATOMS'])==len(cif.chains[key]['LIGANDS'][i]['ATOMS']):
                        diff = [x for x in pdb.chains[key]['LIGANDS'][i] if pdb.chains[key]['LIGANDS'][i][x]!=cif.chains[key]['LIGANDS'][i][x]]
                        if not diff: continue
                        if diff == ['ATOMS']:
                            flag = False
                            for j in range(len(pdb.chains[key]['LIGANDS'][i]['ATOMS'])):
                                diff2 = [x for x in pdb.chains[key]['LIGANDS'][i]['ATOMS'][j]\
                                         if pdb.chains[key]['LIGANDS'][i]['ATOMS'][j][x]!=cif.chains[key]['LIGANDS'][i]['ATOMS'][j][x]]
                                if diff2 and diff2 != ['NUM']:
                                    flag = True
                            if not flag: continue
                        print(pdb.headers['PDBFILE'],pdb.headers['MDLNO'],'Wrong ligand %s in chain %s (keys %s)'%(i,key,','.join(diff)))
                        continue
            elif pdb.chains[key][key2] != cif.chains[key][key2]:
                if key2=='MOL_ID': continue
                print(pdb.headers['PDBFILE'],pdb.headers['MDLNO'],'Wrong key %s in chain %s'%(key2,key))
                print('PDB:',pdb.chains[key][key2])
                print('CIF:',cif.chains[key][key2])
                continue
    if sorted(list(pdb.ids.keys())) != sorted(list(cif.ids.keys())):
        print(pdb.headers['PDBFILE'],pdb.headers['MDLNO'],'Wrong set of ids.keys()')
        return
    for key in pdb.ids:
        if pdb.ids[key] != cif.ids[key]:
            if key in ('NUCL','LIGAND') and pdb.ids['NUCL']+pdb.ids['LIGAND']==cif.ids['NUCL']+cif.ids['LIGAND']: continue
            print(pdb.headers['PDBFILE'],pdb.headers['MDLNO'],'Wrong ids.key="%s"'%key)
            print('PDB:',pdb.ids)
            print('CIF:',cif.ids)
            break

    #print('Models are the same!!')
    return 15

def Check_all():

    files = glob.glob('/home/baulin/eugene/Работа/RSSDB/pdb/models/*.pdb*')

    k = 0
    for file in sorted(files):
        k += 1

        slash = file.rfind('/')
        dot = file.rfind('.')
        ciffile = file[:slash].replace('pdb','mmCIF')+file[slash:dot].lower()+file[dot:].replace('pdb','cif')
        if not os.path.exists(ciffile) or not ciffile.endswith('cif1'): continue
        print(ciffile)

        pdb = PDB.Model(file)
        cif = Model(ciffile)

        if Check(pdb,cif) != 15:
            #continue
            return pdb,cif

    return 'FUCK','YEAH!!'

pdb,cif = Check_all()
#cif = Model('/home/baulin/eugene/Работа/RSSDB/mmCIF/models/4jxz.cif1')
'''
