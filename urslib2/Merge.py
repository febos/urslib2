
import os

try:
    from urslib2 import mmCIF
    from urslib2 import PDB
    from urslib2 import DSSR
    from urslib2 import Tools
except ImportError:
    import mmCIF
    import PDB
    import DSSR
    import Tools


class Model():

    def __init__(self,pdbmodel,outmodel,form='CIF'):

        if os.path.basename(pdbmodel)[:4] != os.path.basename(outmodel)[:4]: raise TypeError('Inappropriate file names!')

        if   form == 'CIF': pdb = mmCIF.Model(pdbmodel)
        elif form == 'PDB': pdb =   PDB.Model(pdbmodel)

        dssr = DSSR.Model(outmodel)

        self.headers = {}

        for key,value in pdb.headers.items():  self.headers[key] = value
        for key,value in dssr.headers.items(): self.headers[key] = value

        self.chains     = pdb.chains
        self.molecules  = pdb.molecules
        self.ids        = pdb.ids
        self.allwords   = pdb.allwords
        self.restype    = pdb.restype
        self.bpairs     = dssr.bpairs
        self.lumults    = dssr.multiplets       
        self.helices    = dssr.helices
        self.lustems    = dssr.stems
        self.lonepairs  = dssr.lonepairs        # will be merged with self.lustems and then deleted
        self.non_pairs  = dssr.non_pairs
        self.lu_loops   = dssr.lu_loops
        self.non_loops  = dssr.non_loops
        self.kissing    = dssr.kissing
        self.a_minors   = dssr.a_minors
        self.u_turns    = dssr.u_turns
        self.ribzips    = dssr.ribzips
        self.k_turns    = dssr.k_turns
        self.phosphates = dssr.phosphates
        self.lubrackets = dssr.brackets
        self.stacks     = dssr.stacks
        self.abcaps     = dssr.abcaps

        self.add_nucl_id()
        self.merge_lustems_lones()
        self.add_lubrackets()
        self.add_bound()
        self.add_nucl_bps()

    def add_nucl_id(self):

        ### Filling self.nucl_id with keys ###
        
        self.nucl_id    = {}

        for a in self.a_minors:   self.nucl_id[a['NUCL']] = 0
        
        for a in self.bpairs:

            self.nucl_id[a['NUCL1']] = 0
            self.nucl_id[a['NUCL2']] = 0
        
        for a in self.lu_loops['BULGE']:

            for n in a['NUCLS']:  self.nucl_id[n] = 0
        
        for a in self.lu_loops['HAIRPIN']:

            for n in a['NUCLS']:  self.nucl_id[n] = 0
        
        for a in self.lu_loops['INTERNAL']:

            for n in a['LNUCLS']: self.nucl_id[n] = 0
            for n in a['RNUCLS']: self.nucl_id[n] = 0
        
        for a in self.lu_loops['JUNCTION']:

            for nn in a['NUCLS']:

                for n in nn:      self.nucl_id[n] = 0
        
        for a in self.lumults:

            for n in a['NUCLS']:  self.nucl_id[n] = 0
        
        for a in self.non_loops:

            self.nucl_id[a['START']] = 0
            self.nucl_id[a['END']]   = 0
        
        for a in self.non_pairs:

            self.nucl_id[a['NUCL1']] = 0
            self.nucl_id[a['NUCL2']] = 0
        
        for a in self.u_turns:

            self.nucl_id[a['NUCL1']] = 0
            self.nucl_id[a['NUCL2']] = 0
        #print(self.nucl_id)
        for a in self.ribzips:

            for n in a['NUCLS']: self.nucl_id[n] = 0

        ### Filling self.nucl_id with values ###

        List, ch, pdbnum = None, None, None

        for nucl in sorted(self.nucl_id.keys()):
            
            List   = nucl.split('.')
            ch     = List[0]
            pdbnum = Tools.pdbnum(List[2]+List[3])

            if ch not in self.chains: ch = self.headers['MASKEDCHS'][ch] # For nucls from fake chains

            if pdbnum in self.chains[ch]['NUMS']: self.nucl_id[nucl] = ['RES', self.chains[ch]['NUMS'][pdbnum]]
            else: # if not in RES

                for i in range(len(self.chains[ch]['LIGANDS'])):

                    if self.chains[ch]['LIGANDS'][i]['PDBNUM'] == pdbnum: self.nucl_id[nucl] = ['LIGANDS', i]

        del List,ch,pdbnum

        ### Adding values from self.nucl_id to self.bpairs   (C.N.P.I -> [C.N.P.I, RES/LIGAND, index])
            
        for bp in self.bpairs:

            bp['NUCL1'] = [bp['NUCL1'],] + self.nucl_id[bp['NUCL1']]
            bp['NUCL2'] = [bp['NUCL2'],] + self.nucl_id[bp['NUCL2']]

        ### Adding values from self.nucl_id to self.ribzips   (C.N.P.I -> [C.N.P.I, RES/LIGAND, index])
        
        for Zip in self.ribzips:

            for i in range(len(Zip['NUCLS'])):

                n = [Zip['NUCLS'][i],] + self.nucl_id[Zip['NUCLS'][i]]
                Zip['NUCLS'][i] = n
                self.chains[n[0].split('.')[0]][n[1]][n[2]]['ZIP'] = Zip['ID']

        ### Adding values from self.nucl_id to self.lumults   (C.N.P.I -> [C.N.P.I, RES/LIGAND, index])

        for Mult in self.lumults:

            for i in range(len(Mult['NUCLS'])):

                n = [Mult['NUCLS'][i],] + self.nucl_id[Mult['NUCLS'][i]]
                Mult['NUCLS'][i] = n
                self.chains[n[0].split('.')[0]][n[1]][n[2]]['LUMULT'] = Mult['ID']

    def merge_lustems_lones(self):

        self.lustems = self.lustems + self.lonepairs        # Merging
        
        del self.lonepairs

        self.lustems.sort(key = lambda x: x['PAIRS'][0])    # Sorting by id of first bp

        stemDSSRid_stemID = {None:None,} # {old_id : new_id,}

        for i in range(len(self.lustems)):

            self.lustems[i]['ID'] = i+1  # adding new id

            stemDSSRid_stemID[self.lustems[i]['DSSRID']] = self.lustems[i]['ID']    # DSSRID == old; ID == new

        for bp       in self.bpairs: bp['LUSTEM'] = stemDSSRid_stemID[bp['LUSTEM']]

        for bulge    in self.lu_loops['BULGE']:

            bulge['PREV'] = stemDSSRid_stemID[bulge['PREV']]
            bulge['NEXT'] = stemDSSRid_stemID[bulge['NEXT']]

        for hairpin  in self.lu_loops['HAIRPIN']: hairpin['CLOSING'] = stemDSSRid_stemID[hairpin['CLOSING']]

        for internal in self.lu_loops['INTERNAL']:

            internal['PREV'] = stemDSSRid_stemID[internal['PREV']]
            internal['NEXT'] = stemDSSRid_stemID[internal['NEXT']]

        for junction in self.lu_loops['JUNCTION']:

            for i in range(len(junction['CLOSING'])):

                junction['CLOSING'][i] = stemDSSRid_stemID[junction['CLOSING'][i]]

        for kiss     in self.kissing: kiss['KISS'] = stemDSSRid_stemID[kiss['KISS']]

        for k_turn   in self.k_turns:

            k_turn['STEM1'] = stemDSSRid_stemID[k_turn['STEM1']]
            k_turn['STEM2'] = stemDSSRid_stemID[k_turn['STEM2']]

        del stemDSSRid_stemID

    def add_lubrackets(self): # adding brackets from dssr to pdb-chains

        for brack in self.lubrackets:

            if brack['CHAIN'] in self.chains :

                self.chains[brack['CHAIN']]['LUBRACKETS'] = brack['DIAGRAM']

    def add_bound(self): # adding number of bound nucleotides to pdb-chains

        seen = {}

        for bp in self.bpairs:

            if bp['NUCL1'][0] not in seen:

                seen[bp['NUCL1'][0]] = 0

                if bp['NUCL1'][1] == 'RES':

                    self.chains[bp['CHAIN1']]['BOUND'] += 1
                    if bp['TYPE'] in ('WC','WB'): self.chains[bp['CHAIN1']]['BOUNDWCWB']+= 1

                elif bp['NUCL1'][1] == 'LIGANDS':

                    self.chains[bp['CHAIN1']]['LIGBOUND'] += 1


            if bp['NUCL2'][0] not in seen:

                seen[bp['NUCL2'][0]] = 0

                if bp['NUCL2'][1] == 'RES':

                    self.chains[bp['CHAIN2']]['BOUND'] += 1
                    if bp['TYPE'] in ('WC','WB'): self.chains[bp['CHAIN2']]['BOUNDWCWB']+= 1

                elif bp['NUCL2'][1] == 'LIGANDS':

                    self.chains[bp['CHAIN2']]['LIGBOUND'] += 1

    def add_nucl_bps(self): # increasing nucl['BPS']

        for bp in self.bpairs:

            for nucl in (bp['NUCL1'], bp['NUCL2']):

                self.chains[nucl[0].split('.')[0]][nucl[1]][nucl[2]]['BPS'] += 1










        
