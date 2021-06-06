import os
import time
import glob
import subprocess

try:
    try:
        from urslib2.config_private import dssr_path
    except:
        from config_private import dssr_path
except:
    try:
        from urslib2.config import dssr_path
    except:
        from config import dssr_path


try: from urslib2 import Tools
except ImportError: import Tools

def run(modelpath,outfolder):

    modelname  = os.path.basename(modelpath)
    outname    = modelname[:5] + 'out' + modelname[8:]
    torname    = modelname[:5] + 'tor' + modelname[8:]
    multname   = modelname[:5] + 'mult' + modelname[8:]
    aminorname = modelname[:5] + 'aminor' + modelname[8:]

    if outfolder[-1] != '/':
        outfolder += '/'    

    args=[dssr_path,
          '-more',
          '-non-pair',
          '-u-turn',
          '-po4',
          '-idstr=long',
          '-i='+modelpath,
          '-o='+outfolder+outname]

    problematic_args = ('-non-pair','-u-turn','-po4') # deleting of this tokens can help to avoid error of dssr

    dssr = subprocess.Popen(args)   # run dssr
    dssr.wait()                     # wait for it
    returnvalue = dssr.returncode   # get return code
    message = ''
    error = 0


    if returnvalue or not glob.glob(outfolder+outname): # if some error

        error = 1

        print ('Some problem with '+ modelname + '\nTrying to add or delete some tokens of dssr:')
        message += 'Some problem with '+ modelname + '. Trying to add or delete some tokens of dssr:'

        for i in ('0','1','2','01','02','12','012'): #different sets of deleted tokens (indexes of problematic_args)

            args2 = []
            for j in args: args2.append(j)

            deleted = []

            for index in i:

                args2.remove(problematic_args[int(index)])
                deleted.append(problematic_args[int(index)])

            dssr = subprocess.Popen(args2)
            dssr.wait()
            returnvalue = dssr.returncode

            if not returnvalue and glob.glob(outfolder+outname): break

        if returnvalue or not glob.glob(outfolder+outname): # if deleting didn't help - trying to add -altloc=(B/C/D/E/F/G)

            for i in 'BCDEFG':

                args2 = []
                for j in args: args2.append(j)
                args2.append('-altloc='+i)

                dssr = subprocess.Popen(args2)
                dssr.wait()
                returnvalue = dssr.returncode

                if not returnvalue and glob.glob(outfolder+outname):

                    added = args2[-1]
                    break

            if returnvalue or not glob.glob(outfolder+outname):
            
                print("It didn't solve the problem. " + outname + ' hasn\'t been created.')
                message += "\nIt didn't solve the problem. " + outname + ' hasn\'t been created.'

            else:

                print('Added: '+added)
                print("It solved the problem! " + outname + ' has been created!')
                message += '\nAdded: '+added
                message += "\nIt solved the problem! " + outname + ' has been created!'

        else:

            print('Deleted: ' + ','.join(deleted))
            print("It solved the problem! " + outname + ' has been created!')
            message += '\nDeleted: ' + ','.join(deleted)
            message += "\nIt solved the problem! " + outname + ' has been created!'

    else: print(outname+'\tis successfully created',end='\n')

    tor = subprocess.Popen(['mv','dssr-torsions.txt',
                            outfolder.replace('out','tor')+torname])   # move torsions file
    tor.wait()

    mult = subprocess.Popen(['mv','dssr-multiplets.pdb',
                            outfolder.replace('out','mult')+multname])   # move multiplets file

    mult.wait()

    aminor = subprocess.Popen(['mv','dssr-Aminors.pdb',
                            outfolder.replace('out','aminor')+aminorname])   # move multiplets file

    aminor.wait()

    return error,message

def run_all(infolder='',outfolder='',form='cif',temp=False):

    good_input = 0

    if temp: good_input = 1

    while not good_input:
        
        if os.path.exists(infolder) and os.path.exists(outfolder): good_input = 1

        elif infolder or outfolder:

            print('Bad paths. Try again...')
            infolder  = input('Please enter the path to folder with %s models: '%form)
            outfolder = input('Please enter the path to output folder: ')

        else:

            infolder  = input('Please enter the path to folder with %s models: '%form)
            outfolder = input('Please enter the path to output folder: ')

    Time = time.time()

    if infolder[-1] != '/':
        infolder += '/'  

    models = glob.glob(infolder+'*.%s*'%form)

    models.sort()

    num     = len(models)
    counter = 1

    errors = []
    messages = []

    for model in models:

        print(str(counter)+'/'+str(num),end='\t')
        error,code = run(model,outfolder)
        counter += 1

        if error:

            errors.append(model)
            messages.append(code)

    if errors:

        print('\n\nProblematic and unhandled models: \n')

        for i in range(len(errors)):

            print(errors[i])
            print(messages[i])

    else: print ('\n\nAll models are processed well\n')

    Time = time.time() - Time    
    print('Time: %s h %s min %s sec'%(int(Time//3600),int((Time%3600)//60),int(Time%60)))


class Model():

    def __init__(self,outpath):

        self.headers    = {'FILE':     '', 'CHAINSNUM': 0,
                           'CHAINS':   [], 'NUCLS':     0,
                           'ATOMS':     0, 'WATERS':    0,
                           'METALSNUM': 0, 'METALS':    {},
                           'BPAIRS':    0, 'LUMULTS':   0,
                           'HELICES':   0, 'STEMS':     0,
                           'LONES':     0, 'NON-PAIRS': 0,
                           'HAIRPINS':  0, 'BULGES':    0,
                           'INTERNALS': 0, 'JUNCTIONS': 0,
                           'PSEUDO':False, 'NON-LOOPS': 0,
                           'KISSINGS':  0, 'A-MINORS':  0,
                           'U-TURNS':   0, 'RIBZIPS':   0,
                           'PHOSPHATES':0, 'K-TURNS':   0}
        self.bpairs     = []
        self.multiplets = []
        self.helices    = []
        self.stems      = []
        self.lonepairs  = []
        self.non_pairs  = []
        self.lu_loops   = {'HAIRPIN': [],
                           'BULGE':   [],
                           'INTERNAL':[],
                           'JUNCTION':[]}
        self.non_loops  = []
        self.kissing    = []
        self.a_minors   = []
        self.u_turns    = []
        self.ribzips    = []
        self.k_turns    = []
        self.phosphates = []
        self.brackets   = []
        self.stacks     = []
        self.abcaps     = []


        self.parse(outpath)

    def parse(self, outpath):

        DSSRdict = {'HEADER':  [], 'BP':      [], 'MULTIPLET':[], 'HELIX':   [],
                    'STEM':    [], 'LONE':    [], 'NON-PAIR': [], 'HAIRPIN': [],
                    'BULGE':   [], 'INTERNAL':[], 'JUNCTION': [], 'PSEUDO':  [],
                    'NON-LOOP':[], 'KISSING': [], 'A-MINOR':  [], 'U-TURN':  [],
                    'ZIPPER':  [], 'K-TURN':  [], 'PO4':      [], 'BRACKETS':[],
                    '0':       [], 'STACKING':[], 'ATOM-BASE':[],}

        with open(outpath) as out:

            current = ''
            Next    =  0

            for line in out:

                if line[:6] == '******':

                    Next = 1
                    continue

                else:

                    if Next:

                        if   'Note: By defa' in line: current = 'HEADER'
                        elif 'base pair'     in line: current = 'BP'
                        elif 'multiplet'     in line: current = 'MULTIPLET'
                        elif ' heli'         in line: current = 'HELIX'
                        elif ' stem'         in line: current = 'STEM'
                        elif ' isolated WC/' in line: current = 'LONE'
                        elif 'non-pairing'   in line: current = 'NON-PAIR'
                        elif 'hairpin'       in line: current = 'HAIRPIN'
                        elif 'bulge'         in line: current = 'BULGE'
                        elif 'internal'      in line: current = 'INTERNAL'
                        elif 'junction'      in line: current = 'JUNCTION'
                        elif 'pseudoknot'    in line: current = 'PSEUDO'
                        elif 'non-loop'      in line: current = 'NON-LOOP'
                        elif 'kissing'       in line: current = 'KISSING'
                        elif 'A-minor'       in line: current = 'A-MINOR'
                        elif 'U-turn'        in line: current = 'U-TURN'
                        elif 'zipper'        in line: current = 'ZIPPER'
                        elif 'kink turn'     in line: current = 'K-TURN'
                        elif 'phosphate'     in line: current = 'PO4'
                        elif 'dot-bracket'   in line: current = 'BRACKETS'
                        elif 'stacks' in line and 'coaxial' not in line: current = 'STACKING'
                        elif 'atom-base cap' in line: current = 'ATOM-BASE'
                        else                        : current = '0'

                        Next = 0

                    DSSRdict[current].append(line)

            

        del current,Next,DSSRdict['0']

        # HEADERS (self.headers)

        headers = DSSRdict['HEADER']
        #print(headers)
        while headers[0][:4] != 'Date': headers = headers[1:]
        
        self.headers['FILE'] = Tools.cut_spaces(headers[1][10:-1])

        opn = headers[2].find('[')

        self.headers['CHAINSNUM'] = int(headers[2][26:opn])
        chains = headers[2][opn+1:-2].split(',')

        if self.headers['CHAINSNUM']:

            for ch in chains: self.headers['CHAINS'].append([ch[:ch.find('=')], int(ch[ch.find('=')+1:])])

        self.headers['NUCLS']  = int(headers[3][23:])
        self.headers['ATOMS']  = int(headers[4][17:])
        self.headers['WATERS'] = int(headers[5][18:])

        opn = headers[6].find('[')

        self.headers['METALSNUM'] = int(headers[6][18:opn])

        metals = headers[6][opn+1:-2].split(',')
        eq = 0

        if self.headers['METALSNUM']:

            for met in metals:

                eq = met.find('=')
                self.headers['METALS'][met[:eq].upper()] = int(met[eq+1:])

        del headers,opn,chains,metals,eq,DSSRdict['HEADER']
        
        # BASE PAIRS (self.bpairs)

        pairs = DSSRdict['BP']

        if pairs: self.headers['BPAIRS'] = int(pairs[0][7:pairs[0].find('base')])

        for i in range(2,len(pairs)-1):

            if pairs[i][3] != ' ': self.bpairs.append(Tools.BP(pairs[i:i+6]))

        del pairs,DSSRdict['BP']

        # MULTIPLETS (self.multiplets)

        multiplets = DSSRdict['MULTIPLET']

        if multiplets: self.headers['LUMULTS'] = int(multiplets[0][7:multiplets[0].find('multi')])

        for i in range(1,len(multiplets)-1):

            self.multiplets.append(Tools.LuMult(multiplets[i]))

        del multiplets,DSSRdict['MULTIPLET']

        # STEMS (self.stems)

        steps = {}

        stems = DSSRdict['STEM']

        if stems: self.headers['STEMS'] = int(stems[0][7:stems[0].find('stem')])

        start, end = 0,0

        for i in range(6,len(stems)):

            if stems[i] == '\n' or stems[i][:6] == '  ----':  # end of stem

                end = i
                st, ststeps = Tools.LuStem(stems[start:end])
                self.stems.append(st)
                for st in ststeps:
                    if st not in steps: steps[st] = ststeps[st]

            elif stems[i][6] == '#': start = i #start of stem

        del start,end,stems,DSSRdict['STEM']

        # HELICES (self.helices)

        helices = DSSRdict['HELIX']

        if helices: self.headers['HELICES'] = int(helices[0][7:helices[0].find('heli')])

        start, end = 0,0

        for i in range(11,len(helices)):

            if helices[i] == '\n' or helices[i][:6] == '  ----':  # end of helix

                end = i
                hel,helsteps =  Tools.LuHelix(helices[start:end])
                self.helices.append(hel)
                for st in helsteps:
                    if st not in steps: steps[st] = helsteps[st]

            elif helices[i][7] == '#': start = i #start of helix

        del start,end,helices,DSSRdict['HELIX']
        
        # LONE PAIRS (self.lonepairs)

        lpairs = DSSRdict['LONE']

        if lpairs: self.headers['LONES'] = int(lpairs[0][7:lpairs[0].find('isolated ')])

        for pair in lpairs[4:-1]: self.lonepairs.append(Tools.LuLone(pair))

        del lpairs,DSSRdict['LONE']

        # PAIR-LONE, PAIR-STEM and PAIR-HELIX relations

        pair_lone  = {}
        pair_stem  = {}
        pair_id    = {'\\N':'\\N',}
        lone_helix = {}
        pair_helix = {}

        for helix in self.helices:

            for pair in helix['PAIRS']: pair_helix[pair] = helix['ID']

        for stem  in self.stems:

            for pair in stem['PAIRS']:  pair_stem[pair]  = stem['ID']

        for lone in self.lonepairs:

            pair  = lone['PAIRS'][0]
            pair2 = pair.split('.')
            pair2 = '.'.join(pair2[4:]) + '.' + '.'.join(pair2[:4])

            pair_lone[pair] = lone['ID']

        for bp in self.bpairs:

            pair  = bp['NUCL1'] + '.' + bp['NUCL2']
            pair2 = bp['NUCL2'] + '.' + bp['NUCL1']

            pair_id[pair]  = bp['ID']
            pair_id[pair2] = bp['ID']

            if   pair  in pair_helix: bp['HELIX']  = pair_helix[pair]
            elif pair2 in pair_helix: bp['HELIX']  = pair_helix[pair2]
            if   pair  in pair_stem:  bp['LUSTEM'] = pair_stem[pair]
            elif pair2 in pair_stem:  bp['LUSTEM'] = pair_stem[pair2]
            elif pair  in pair_lone:  bp['LUSTEM'] = pair_lone[pair]
            elif pair2 in pair_lone:  bp['LUSTEM'] = pair_lone[pair2]

        for pair in steps: self.bpairs[pair_id[pair]-1]['STEP'] = steps[pair]

        del steps

        for lone in self.lonepairs: lone['PAIRS'][0] = pair_id[lone['PAIRS'][0]]

        for stem in self.stems:

            for i in range(len(stem['PAIRS'])):  stem['PAIRS'][i]  = pair_id[stem['PAIRS'][i]]

        for helix in self.helices:

            for i in range(len(helix['PAIRS'])): helix['PAIRS'][i] = pair_id[helix['PAIRS'][i]]

        del pair_lone,pair_stem,lone_helix,pair_helix

        # STEMS and LONES for HELICES (number of stems and lones in a helix)

        for stem in self.stems:
            if stem['HELIX']: self.helices[stem['HELIX']-1]['STEMS'] += 1
            for hhel in stem['OTHERHELICES']: self.helices[hhel-1]['STEMS'] += 1

        for lone in self.lonepairs:

            if lone['HELIX']: self.helices[lone['HELIX']-1]['LONES'] += 1
            for hhel in lone['OTHERHELICES']: self.helices[hhel-1]['LONES'] += 1

        # NON-PAIRING INTERACTIONS (self.non_pairs)

        npairs = DSSRdict['NON-PAIR']

        if npairs: self.headers['NON-PAIRS'] = int(npairs[0][7:npairs[0].find('non-')])

        for pair in npairs[1:-1]: self.non_pairs.append(Tools.LuNon_pair(pair))

        del npairs,DSSRdict['NON-PAIR']

        """ LOOPS (self.lu_loops) """
        # HAIRPINS (self.lu_loops['HAIRPIN'])

        hairpins = DSSRdict['HAIRPIN']

        if hairpins: self.headers['HAIRPINS'] = int(hairpins[0][7:hairpins[0].find('hair')])

        for i in range(1,len(hairpins)-2):

            if hairpins[i][3] != ' ': self.lu_loops['HAIRPIN'].append(Tools.LuHairpin(hairpins[i:i+4]))

        del hairpins,DSSRdict['HAIRPIN']
                                                                      
        # BULGES (self.lu_loops['BULGE'])

        bulges = DSSRdict['BULGE']

        if bulges: self.headers['BULGES'] = int(bulges[0][7:bulges[0].find('bulg')])

        for i in range(1,len(bulges)-4):

            if bulges[i][5] == 'b': self.lu_loops['BULGE'].append(Tools.LuBulge(bulges[i:i+5]))

        del bulges,DSSRdict['BULGE']

        # INTERNAL LOOPS (self.lu_loops['INTERNAL'])

        internals = DSSRdict['INTERNAL']

        if internals: self.headers['INTERNALS'] = int(internals[0][7:internals[0].find('inter')])

        for i in range(1,len(internals)-4):

            if internals[i][3] != ' ': self.lu_loops['INTERNAL'].append(Tools.LuInternal(internals[i:i+5]))

        del internals,DSSRdict['INTERNAL']

        # JUNCTIONS (self.lu_loops['JUNCTION'])

        junctions = DSSRdict['JUNCTION']
        hyph, num = 0, 0

        if junctions: self.headers['JUNCTIONS'] = int(junctions[0][7:junctions[0].find('junct')])

        for i in range(1,len(junctions)-1):

            if junctions[i][3] != ' ':

                hyph = junctions[i].find('-')
                num  = int(junctions[i][5:hyph])
                self.lu_loops['JUNCTION'].append(Tools.LuJunction(junctions[i:i+num+3]))

        del hyph,num,junctions,DSSRdict['JUNCTION']

        # PSEUDO (self.headers['PSEUDO'])

        if  DSSRdict['PSEUDO']: self.headers['PSEUDO'] = True
        del DSSRdict['PSEUDO']

        # NON-LOOP SEGMENTS (self.non_loops)

        non_loops = DSSRdict['NON-LOOP']

        if non_loops: self.headers['NON-LOOPS'] = int(non_loops[0][7:non_loops[0].find('non-l')])

        for line in non_loops[1:-1]: self.non_loops.append(Tools.LuNon_loop(line))

        del non_loops,DSSRdict['NON-LOOP']

        # KISSING LOOPS (self.kissing)

        kissing = DSSRdict['KISSING']

        if kissing: self.headers['KISSINGS'] = int(kissing[0][7:kissing[0].find('kiss')])

        for line in kissing[1:-1]: self.kissing.append(Tools.LuKissing(line))

        del kissing,DSSRdict['KISSING']

        # A-MINOR MOTIFS (self.a_minors)

        a_minors = DSSRdict['A-MINOR']

        if a_minors: self.headers['A-MINORS'] = int(a_minors[0][7:a_minors[0].find('A-min')])

        for i in range(1,len(a_minors)-3):

            if a_minors[i][3] != ' ':

                self.a_minors.append(Tools.LuA_minor(a_minors[i:i+3]))
                self.a_minors[-1]['PAIR'] = pair_id[self.a_minors[-1]['PAIR']]

        del a_minors,DSSRdict['A-MINOR']

        # U-TURNS (self.u_turns)
        
        u_turns = DSSRdict['U-TURN']

        if u_turns: self.headers['U-TURNS'] = int(u_turns[0][7:u_turns[0].find('U-turn')])

        for line in u_turns[1:-1]: self.u_turns.append(Tools.LuU_turn(line))

        del u_turns,DSSRdict['U-TURN']

        # RIBOSE ZIPPERS (self.ribzips)

        zippers = DSSRdict['ZIPPER']

        if zippers: self.headers['RIBZIPS'] = int(zippers[0][7:zippers[0].find('ribose')])

        for line in zippers[1:-1]: self.ribzips.append(Tools.LuZipper(line))

        del zippers,DSSRdict['ZIPPER']

        # POSSIBLE KINK TURNS (self.k_turns)

        k_turns = DSSRdict['K-TURN']

        if k_turns: self.headers['K-TURNS'] = int(k_turns[0][7:k_turns[0].find('possible')])

        for i in range(1,len(k_turns)-1):

            if k_turns[i][3] != ' ':

                self.k_turns.append(Tools.LuK_turn(k_turns[i:i+4]))
                self.k_turns[-1]['PAIR'] = pair_id[self.k_turns[-1]['PAIR']]
                
        del k_turns,pair_id,DSSRdict['K-TURN']

        # STACKS (self.stacks)

        stacks = DSSRdict['STACKING']

        for s in stacks:

            if 'nts=' in s:

                li = s.split()
                self.stacks.append([int(li[1][4:]),['.'.join(x.split('.')[2:]) for x in li[3].split(',')]])

        # PHOSPHATE INTERACTIONS (self.phosphates)

        """ LATER """
        del DSSRdict['PO4']

        # ATOM-BASE CAPPING

        abcaps = DSSRdict['ATOM-BASE']

        for i in range(4,len(abcaps)-1):

            self.abcaps.append(Tools.LuAtomBaseCap(abcaps[i]))
        
        # BRACKETS (self.brackets)

        brackets = DSSRdict['BRACKETS']

        for i in range(len(brackets)-1):

            if brackets[i][:5] == '>'+self.headers['FILE'][:4]:
                self.brackets.append(Tools.LuBrackets(brackets[i:i+3]))

        del brackets,DSSRdict['BRACKETS']
        

if __name__ == '__main__':

    run_all()












