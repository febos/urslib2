

def Organism(organism):

    org = organism.lower()

    org = org.replace('bean pod','bean-pod')
    org = org.replace('cytoplasmic polyhedrosis virus','cypovirus')
    org = org.replace('spinacia','spinacea')
    org = org.replace('eschericia','escherichia')
    org = org.replace('enterobacterio','enterobacteria')
    org = org.replace('type 1','1')
    org = org.replace('type 2','2')
    org = org.replace('[','')
    org = org.replace(']','')

    while '  ' in org: org = org.replace('  ',' ')

    orgsplit = org.split(' ')

    singles = {'acidaminococcus':'acidaminococcus',
               'bacterium': 'unidentified bacterium',
               'canis':'canis',
               'coxsackievirus':'coxsackievirus',
               'cypovirus': 'cypovirus',
               'deltaproteobacteria':'deltaproteobacteria',
               'drosophila':'drosophila',
               'enterovirus':'enterovirus',
               'geobacter':'geobacter',
               'lambdapapillomavirus':'lambdapapillomavirus',
               'lampyridae':'lampyridae',
               'leishmania':'leishmania',
               'metagenome':'metagenome',
               'neurospora':'neurospora',
               'none':'unspecified',
               'norovirus':'norovirus',
               '\\N':'unspecified',
               '':'unspecified',
               'proteobacteria':'proteobacteria',
               'pyrococcus':'pyrococcus',
               'rattus':'rattus',
               'reovirus':'reovirus',
               't7likevirus':'t7likevirus',
               'thermus':'thermus',
               'triticum':'triticum',
               'saccharomyces':'saccharomyces',
               'streptomyces':'streptomyces',
               'synthetic':'synthetic',
               'unclassified':'unidentified',
               'unidentified':'unidentified'}

    if len(orgsplit)==1: return singles[org]

    doubles = {'synthetic dna':'synthetic',
               'synthetic construct':'synthetic',
               'uncultured organism':'unidentified',
               'marine metagenome':'metagenome'}

    if len(orgsplit)==2:
        if org in doubles: return doubles[org]
        if orgsplit[1] in ('sp','sp.'): return orgsplit[0]
        if 'virus' in orgsplit[0]: return orgsplit[0]
        return org

    if 'influenza' in org and 'virus' in org: return 'influenza virus'

    triples = {'canis lupus familiaris': 'canis familiaris',
               'simian 11 rotavirus': 'simian rotavirus',
               'bunyavirus la crosse': 'bunyavirus la crosse'}

    if len(orgsplit)==3:

        if org in triples: return triples[org]
        if 'synthetic' in org: return 'synthetic'
        if orgsplit[1] == 'phage' or orgsplit[2] == 'phage' or 'virus' in orgsplit[2]: return org
        if orgsplit[1] in ('sp','sp.'): return orgsplit[0]
        return ' '.join(orgsplit[:2])

    quadruples = {'rna transcription vector pbrdi1': 'rna transcription vector pbrdi1',
                  'mopeia lassa reassortant 29':'mopeia lassa reassortant'}

    if len(orgsplit)==4:
        if org in quadruples: return quadruples[org]
        if 'synthetic' in org: return 'synthetic'
        if orgsplit[3]=='virus': return org
        if orgsplit[1] in ('sp','sp.','22'): return orgsplit[0]
        if orgsplit[2] in ('subsp','subsp.'): return ' '.join(orgsplit[:2])
        if 'virus' in orgsplit[2]: return ' '.join(orgsplit[:3])
        if ',' in org: return org.split(',')[0]
        return ' '.join(orgsplit[:2])

    multiples = {'insect cell expression vector ptie1':'insect cell expression vector ptie1'}

    if len(orgsplit)>4:
        if org in multiples: return multiples[org]
        if 'synthetic' in org: return 'synthetic'
        for i in range(len(orgsplit)):
            if orgsplit[i]=='virus': return ' '.join(orgsplit[:i+1])
        if orgsplit[1] in ('sp','sp.'): return orgsplit[0]
        if orgsplit[2] in ('subsp','subsp.','var','var.','str','str.'): return ' '.join(orgsplit[:2])
        if orgsplit[2][0]=='(': return ' '.join(orgsplit[:2])
        return ' '.join(orgsplit[:2])


def add(model):

    for mol in model.molecules:
        model.molecules[mol]['ORGCAT'] = Organism(model.molecules[mol]['ORGANISM_SCIENTIFIC'])

