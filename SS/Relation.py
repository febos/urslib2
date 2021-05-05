 




def DSSRdict(model): # add dict {dssr: (chain,index)}

    dssrnucls = {}
    for ch in model.chains:
        for place in ('RES','LIGANDS'):
            for i in range(len(model.chains[ch][place])):
                dssrnucls[model.chains[ch][place][i]['DSSR']] = [ch,place,i]

    model.dssrnucls = dssrnucls

def ModifyStemsLoops(model): #add neighbors to stems and loops + add list 'LOOPS' to threads

    for lt in model.loops:
        for i in range(len(model.loops[lt])):
            loop = model.loops[lt][i]
            for tloop in loop['TLOOP']:
                if 'LOOPS' not in model.threads[tloop['THREAD']-1]:
                    model.threads[tloop['THREAD']-1]['LOOPS'] = []
                model.threads[tloop['THREAD']-1]['LOOPS'].append([lt,loop['ID'],loop['PTYPE']])
            stems = []
            stems.append(loop['STEM'])
            for wloop in loop['WLOOP']:
                stems.append(model.wings['LU'][wloop['WING']-1]['STEM'])
            for floop in loop['FLOOP']:
                stems.append(floop['STEM1'])
                stems.append(floop['STEM2'])
                
            model.loops[lt][i]['NEIGHBORS'] = [x for x in list(set(stems)) if x!='\\N']
            for s in model.loops[lt][i]['NEIGHBORS']:
                if 'NEIGHBORS' not in model.stems[s-1]:
                    model.stems[s-1]['NEIGHBORS'] = []
                model.stems[s-1]['NEIGHBORS'].append([lt,loop['ID'],loop['PTYPE']])

    for i in range(len(model.threads)):

        model.threads[i]['NEIGHBORS'] = []

        if 'LOOPS' not in model.threads[i]: model.threads[i]['LOOPS'] = []
        
        for loop in model.threads[i]['LOOPS']:
            model.threads[i]['NEIGHBORS'] += model.loops[loop[0]][loop[1]-1]['NEIGHBORS']
        model.threads[i]['NEIGHBORS'] = list(set(model.threads[i]['NEIGHBORS']))

    for i in range(len(model.wings['LU'])):

        stem_i = model.wings['LU'][i]['STEM']
        model.wings['LU'][i]['NEIGHBORS'] = model.stems[stem_i - 1]['NEIGHBORS']


def NuclSS(self,dssr):

    ch,pl,i = self.dssrnucls[dssr]
    if self.chains[ch][pl][i]['WING']: return 'S'
    elif self.chains[ch][pl][i]['THREAD']:
        t = self.chains[ch][pl][i]['THREAD']
        res = []
        if not self.threads[t-1]['LOOPS']: return 'NA'
        for loop in self.threads[t-1]['LOOPS']:
            res.append(loop[0][0]+loop[2])
        return ''.join(sorted(list(set(res))))
    else:
        return 'NA'
    

def NuclRelation(self,dssr1,dssr2): # dssr1,dssr2 -> SM/LC/NR/LR

    ch1,pl1,i1 = self.dssrnucls[dssr1]
    ch2,pl2,i2 = self.dssrnucls[dssr2]

    w1 = self.chains[ch1][pl1][i1]['WING']
    t1 = self.chains[ch1][pl1][i1]['THREAD']
    w2 = self.chains[ch2][pl2][i2]['WING']
    t2 = self.chains[ch2][pl2][i2]['THREAD']

    if (not w1 and not t1) or (not w2 and not t2):

        return 'NA'

    if w1 and w2:


        if self.wings['LU'][w1-1]['STEM'] == self.wings['LU'][w2-1]['STEM']:

            return 'SM'
        
        else:

            neibs1 = set([tuple(nei) for nei in self.wings['LU'][w1-1]['NEIGHBORS']])
            neibs2 = set([tuple(nei) for nei in self.wings['LU'][w2-1]['NEIGHBORS']]) 

            if neibs1 & neibs2:

                return 'NR'

            return 'LR'

    if w1 and t2:

        if self.wings['LU'][w1-1]['STEM'] in self.threads[t2-1]['NEIGHBORS']:

            return 'LC'
        
        else:

            return 'LR'

    if t1 and w2:

        if self.wings['LU'][w2-1]['STEM'] in self.threads[t1-1]['NEIGHBORS']:

            return 'LC'

        else:

            return 'LR'

    if t1 and t2:

        if t1 == t2:

            return 'SM'

        loops1 = [x[0]+str(x[1]) for x in self.threads[t1-1]['LOOPS']]
        loops2 = [x[0]+str(x[1]) for x in self.threads[t2-1]['LOOPS']]

        if set(loops1) & set(loops2):

            return 'SM'

        else:

            neibs1 = set(self.threads[t1-1]['NEIGHBORS'])
            neibs2 = set(self.threads[t2-1]['NEIGHBORS'])

            if neibs1 & neibs2:

                return 'NR'

            return 'LR'

def AnnotateLinks(model):

    for link in model.links:

        dssr1 = link['NUCL1'][0]
        dssr2 = link['NUCL2'][0]

        link['SS1'] = model.NuclSS(dssr1)
        link['SS2'] = model.NuclSS(dssr2)
        link['REL'] = model.NuclRelation(dssr1,dssr2)


def add(model):

    import types
    DSSRdict(model)
    ModifyStemsLoops(model)
    model.NuclRelation = types.MethodType(NuclRelation,model)
    model.NuclSS = types.MethodType(NuclSS,model)

    AnnotateLinks(model)

    
