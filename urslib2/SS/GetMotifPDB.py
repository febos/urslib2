 



def GetMotifPDB(self,dssrlist,output):

    result = {}

    outp = open(output,'w')

    for dssr in dssrlist:

        ch,pl,i = self.dssrnucls[dssr]

        for atom in self.chains[ch][pl][i]['ATOMS']:

            res = 'ATOM  '
            res += ' '*(5-len(str(atom['NUM'])))
            res += str(atom['NUM'])
            res += '  '
            res += atom['NAME']
            res += ' '*(3-len(str(atom['NAME'])))
            while len(res)<16: res+=' '
            res += atom['ALTLOC']
            res += ' '*(3-len(str(atom['RESNAME'])))
            res += atom['RESNAME']
            res += ' '
            res += atom['CHAIN'][:1]
            res += ' '*(4-len(str(atom['RESNUM'][:-1])))
            res += atom['RESNUM']
            res += ' '*3
            res += ' '*(8-len(str(atom['X'])))
            res += str(atom['X'])
            res += ' '*(8-len(str(atom['Y'])))
            res += str(atom['Y'])
            res += ' '*(8-len(str(atom['Z'])))
            res += str(atom['Z'])
            res += ' '*(6-len(str(atom['OCCUP'])))
            res += str(atom['OCCUP'])
            res += ' '*(6-len(str(atom['TEMPF'])))
            res += str(atom['TEMPF'])
            while len(res)<76: res+=' '
            res += atom['ELEM']
            while len(res)<80: res+=' '

            result[atom['NUM']] = res

    for an in sorted(result.keys()): outp.write(result[an]+'\n')

    outp.close()   




def add(model):

    import types
    model.GetMotifPDB = types.MethodType(GetMotifPDB,model)


'''{'CHAIN': 'R',
 'BONDS': 0,
 'RESNAME': 'G',
 'ID': 8231,
 'CIFID': 17,
 'Y': 56.232,
 'X': 119.527,
 'Z': -36.631,
 'NUM': 339,
 'OCCUP': 1.0,
 'TYPE': 'RNA',
 'RESNUM': '617 ',
 'NAME': 'P',
 'ELEM': 'P',
 'ALTLOC': ' '}
'TEMPF': '33.36'
'''

'''
ATOM    339  P     G R 617     119.527  56.232 -36.631   1.0 33.36          P   
ATOM    340  OP1   G R 617      119.99  57.138 -37.725   1.0 35.94          O   
ATOM    341  OP2   G R 617     118.097  56.451  -36.33   1.0 43.29          O   
ATOM    342  O5'   G R 617     119.737  54.726 -37.113   1.0 35.74          O   

ATOM    339  P     G R 617     119.527  56.232 -36.631  1.00 33.36           P  
ATOM    340  OP1   G R 617     119.990  57.138 -37.725  1.00 35.94           O  
ATOM    341  OP2   G R 617     118.097  56.451 -36.330  1.00 43.29           O  
ATOM    342  O5'   G R 617     119.737  54.726 -37.113  1.00 35.74           O  
ATOM    343  C5'   G R 617     118.854  54.141 -38.087  1.00 38.29           C  
ATOM    344  C4'   G R 617     118.644  52.665 -37.811  1.00 40.94           C  
ATOM    345  O4'   G R 617     117.866  52.066 -38.870  1.00 43.28           O  
ATOM    346  C3'   G R 617     119.969  51.905 -37.777  1.00 41.36           C  

'''
