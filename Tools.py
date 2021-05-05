# this module is a set of auxiliary functions (made by febos)

bptype_dict = {'WC'          :'WC', 'Imino'    :'IM', 'Wobble'    :'WB',
               'n'           :'NA', 'Platform' :'PL', 'Sheared'   :'SH',
               'rHoogsteen'  :'RH', 'Hoogsteen':'HG', 'rWC'       :'RW',
               'rWC (Levitt)':'LV', 'Linker'   :'LN', 'Calcutta'  :'CL',
               'CA_loop'     :'CA', '--'       :'NA', 'rWobble'   :'RB',
               'Metal'       :'MT', '~Wobble'  :'-W', '~Hoogsteen':'-H',
               '~rHoogsteen' :'-R', '~Shear'   :'-S', '~rWobble'  :'-B',
               '~Sheared'    :'-S'}

CIFTOPDB_Mol = {'_entity_src_nat':      {'pdbx_fragment'           : 'FRAGMENT',
                                         'pdbx_organism_scientific': 'ORGANISM_SCIENTIFIC',
                                         'common_name'             : 'ORGANISM_COMMON',
                                         'pdbx_ncbi_taxonomy_id'   : 'ORGANISM_TAXID',
                                         'strain'                  : 'STRAIN',
                                         'pdbx_variant'            : 'VARIANT',
                                         'pdbx_cell_line'          : 'CELL_LINE',
                                         'pdbx_atcc'               : 'ATCC',
                                         'pdbx_organ'              : 'ORGAN',
                                         'tissue'                  : 'TISSUE',
                                         'pdbx_cell'               : 'CELL',
                                         'pdbx_organelle'          : 'ORGANELLE',
                                         'pdbx_secretion'          : 'SECRETION',
                                         'pdbx_cellular_location'  : 'CELLULAR_LOCATION',
                                         'pdbx_plasmid_name'       : 'PLASMID',
                                         'details'                 : 'OTHER_DETAILS'},
                '_entity':              {'pdbx_description' : 'MOLECULE',
                                         'pdbx_fragment'    : 'FRAGMENT',
                                         'pdbx_ec'          : 'EC',
                                         'pdbx_mutation'    : 'MUTATION',
                                         'details'          : 'OTHER_DETAILS',
                                         'src_method'       : 'SYNTHETIC'},
                '_entity_name_com':     {'name': 'SYNONYM'},
                '_pdbx_entity_src_syn': {'organism_scientific' : 'ORGANISM_SCIENTIFIC',
                                         'organism_common_name': 'ORGANISM_COMMON',
                                         'ncbi_taxonomy_id'    : 'ORGANISM_TAXID',
                                         'details'             : 'OTHER_DETAILS'},
                '_entity_src_gen':      {'pdbx_gene_src_fragment'         : 'FRAGMENT',
                                         'pdbx_gene_src_scientific_name'  : 'ORGANISM_SCIENTIFIC',
                                         'gene_src_common_name'           : 'ORGANISM_COMMON',
                                         'pdbx_gene_src_ncbi_taxonomy_id' : 'ORGANISM_TAXID',
                                         'gene_src_strain'                : 'STRAIN',
                                         'pdbx_gene_src_variant'          : 'VARIANT',
                                         'pdbx_gene_src_cell_line'        : 'CELL_LINE',
                                         'pdbx_gene_src_atcc'             : 'ATCC',
                                         'pdbx_gene_src_organ'            : 'ORGAN',
                                         'gene_src_tissue'                : 'TISSUE',
                                         'pdbx_gene_src_cell'             : 'CELL',
                                         'pdbx_gene_src_organelle'        : 'ORGANELLE',
                                         'pdbx_gene_src_cellular_location': 'CELLULAR_LOCATION',
                                         'pdbx_gene_src_gene'             : 'GENE',
                                         'pdbx_host_org_scientific_name'  : 'EXPRESSION_SYSTEM',
                                         'host_org_common_name'           : 'EXPRESSION_SYSTEM_COMMON',
                                         'pdbx_host_org_ncbi_taxonomy_id' : 'EXPRESSION_SYSTEM_TAXID',
                                         'pdbx_host_org_strain'           : 'EXPRESSION_SYSTEM_STRAIN',
                                         'pdbx_host_org_variant'          : 'EXPRESSION_SYSTEM_VARIANT',
                                         'pdbx_host_org_cell_line'        : 'EXPRESSION_SYSTEM_CELL_LINE',
                                         'pdbx_host_org_atcc'             : 'EXPRESSION_SYSTEM_ATCC_NUMBER',
                                         'pdbx_host_org_organ'            : 'EXPRESSION_SYSTEM_ORGAN',
                                         'pdbx_host_org_tissue'           : 'EXPRESSION_SYSTEM_TISSUE',
                                         'pdbx_host_org_cell'             : 'EXPRESSION_SYSTEM_CELL',
                                         'pdbx_host_org_organelle'        : 'EXPRESSION_SYSTEM_ORGANELLE',
                                         'pdbx_host_org_cellular_location': 'EXPRESSION_SYSTEM_CELLULAR_LOCATION',
                                         'pdbx_host_org_vector_type'      : 'EXPRESSION_SYSTEM_VECTOR_TYPE',
                                         'pdbx_host_org_vector'           : 'EXPRESSION_SYSTEM_VECTOR',
                                         'plasmid_name'                   : 'EXPRESSION_SYSTEM_PLASMID',
                                         'pdbx_host_org_gene'             : 'EXPRESSION_SYSTEM_GENE',
                                         'pdbx_description'               : 'OTHER_DETAILS'}}


# dictionary to convert pdb number to float number
letter_weight = {'A':0.03, 'B':0.06, 'C':0.09, 'D':0.12, 'E':0.15,
                 'F':0.18, 'G':0.21, 'H':0.24, 'I':0.27, 'J':0.30,
                 'K':0.33, 'L':0.36, 'M':0.39, 'N':0.42, 'O':0.45,
                 'P':0.48, 'Q':0.51, 'R':0.54, 'S':0.57, 'T':0.60,
                 'U':0.63, 'V':0.66, 'W':0.69, 'X':0.72, 'Y':0.75,
                 'Z':0.78, ' ':0.00}

# dictionary of month to handle date header of pdb files
months = {'JAN':'01', 'FEB':'02', 'MAR':'03', 'APR':'04',
          'MAY':'05', 'JUN':'06', 'JUL':'07', 'AUG':'08',
          'SEP':'09', 'OCT':'10', 'NOV':'11', 'DEC':'12',}

# dictionary of residue's types
restype = { '+T':'DNA', '2DT':'DNA', '2PR':'DNA', '3DA':'DNA', '5AA':'DNA',
           '5AT':'DNA', '5CG':'DNA', '5IT':'DNA', '5IU':'DNA', '5PC':'DNA',
           '8OG':'DNA', 'A5L':'DNA', 'APN':'DNA',  'AS':'DNA', 'BGM':'DNA',
           'BGR':'DNA', 'BRU':'DNA', 'C6G':'DNA', 'CAR':'DNA', 'CBR':'DNA',
           'CFL':'DNA', 'CPN':'DNA', 'CSL':'DNA',  'DA':'DNA', 'DAD':'DNA',
            'DC':'DNA', 'DCP':'DNA', 'DCT':'DNA', 'DCZ':'DNA', 'DDG':'DNA',
           'DFC':'DNA', 'DFG':'DNA',  'DG':'DNA', 'DG3':'DNA', 'DGT':'DNA',
            'DI':'DNA', 'DOC':'DNA', 'DRT':'DNA',  'DT':'DNA',  'DU':'DNA',
           'DUT':'DNA', 'GFL':'DNA', 'GPN':'DNA',  'GS':'DNA', 'LCC':'DNA',
           'MMT':'DNA', 'PDU':'DNA', 'PGN':'DNA', 'PST':'DNA',  'SC':'DNA',
           'SDG':'DNA', 'T3P':'DNA', 'TAF':'DNA', 'TCP':'DNA', 'THY':'DNA',
           'TPN':'DNA',  'TT':'DNA', 'TTP':'DNA', 'U33':'DNA', 'UMS':'DNA',
           'XAR':'DNA', 'XCR':'DNA', 'XGR':'DNA', 'XTR':'DNA', 'XUG':'DNA',
            '+A':'RNA',  '+C':'RNA',  '+G':'RNA',  '+I':'RNA',  '+U':'RNA',
            '0A':'RNA',  '0C':'RNA',  '0G':'RNA',  '0U':'RNA', '10C':'RNA',
           '125':'RNA', '126':'RNA', '127':'RNA', '12A':'RNA', '1AP':'RNA',
           '1DP':'RNA', '1MA':'RNA', '1MG':'RNA', '1RN':'RNA', '1SC':'RNA',
           '23G':'RNA', '2AD':'RNA', '2AR':'RNA', '2AU':'RNA', '2BA':'RNA',
           '2BP':'RNA', '2DA':'RNA', '2IA':'RNA', '2MA':'RNA', '2MG':'RNA',
           '2MU':'RNA', '31H':'RNA', '31M':'RNA', '3AD':'RNA', '3AU':'RNA',
           '3AY':'RNA', '3TD':'RNA', '4OC':'RNA', '4SC':'RNA', '4SU':'RNA',
           '55C':'RNA', '5BU':'RNA', '5CF':'RNA', '5CM':'RNA', '5FU':'RNA',
           '5GP':'RNA', '5IC':'RNA', '5MC':'RNA', '5MU':'RNA', '5NC':'RNA',
           '6AP':'RNA', '6CT':'RNA', '6FC':'RNA', '6FU':'RNA', '6GO':'RNA',
           '6GU':'RNA', '6HA':'RNA', '6HC':'RNA', '6HG':'RNA', '6HT':'RNA',
           '6IA':'RNA', '6MP':'RNA', '6MZ':'RNA', '70U':'RNA', '7AT':'RNA',
           '7MG':'RNA', '8AN':'RNA', '9DG':'RNA',   'A':'RNA', 'A23':'RNA',
           'A2F':'RNA', 'A2L':'RNA', 'A2M':'RNA', 'A3P':'RNA', 'A44':'RNA',
           'A5M':'RNA', 'A5O':'RNA', 'A6A':'RNA', 'A6C':'RNA', 'A6G':'RNA',
           'A6U':'RNA', 'A7E':'RNA', 'A9Z':'RNA', 'ABR':'RNA', 'ABS':'RNA',
           'AD2':'RNA', 'ADE':'RNA', 'ADI':'RNA', 'ADP':'RNA', 'AET':'RNA',
           'AF2':'RNA', 'AG9':'RNA', 'AMD':'RNA', 'AMO':'RNA', 'AMP':'RNA',
           'ANP':'RNA', 'ANZ':'RNA', 'AP7':'RNA', 'APC':'RNA', 'AT7':'RNA',
           'ATD':'RNA', 'ATL':'RNA', 'ATP':'RNA', 'AVC':'RNA', 'AZT':'RNA',
           'BLS':'RNA', 'BRO':'RNA',   'C':'RNA', 'C2E':'RNA', 'C2L':'RNA',
           'C43':'RNA', 'C5L':'RNA', 'C5P':'RNA', 'CB2':'RNA', 'CBV':'RNA',
           'CCC':'RNA', 'CF2':'RNA', 'CFZ':'RNA', 'CG1':'RNA',  'CH':'RNA',
           'CH1':'RNA', 'CIR':'RNA', 'CM0':'RNA', 'CMP':'RNA', 'CMR':'RNA',
           'CP1':'RNA', 'CSG':'RNA', 'CYT':'RNA', 'DGP':'RNA', 'DHU':'RNA',
           'DNR':'RNA', 'DX4':'RNA',   'E':'RNA', 'EDA':'RNA', 'EDC':'RNA',
           'EEM':'RNA', 'EPE':'RNA', 'F3N':'RNA', 'F3O':'RNA', 'FAG':'RNA',
           'FHU':'RNA', 'FMU':'RNA', 'FYA':'RNA',   'G':'RNA', 'G0B':'RNA',
           'G25':'RNA', 'G2L':'RNA', 'G2P':'RNA', 'G46':'RNA', 'G48':'RNA',
           'G7M':'RNA', 'GAO':'RNA', 'GCK':'RNA', 'GCP':'RNA', 'GDO':'RNA',
           'GDP':'RNA', 'GF2':'RNA', 'GH3':'RNA', 'GMP':'RNA', 'GN7':'RNA',
           'GNG':'RNA', 'GNP':'RNA', 'GOM':'RNA', 'GRB':'RNA', 'GSR':'RNA',
           'GSS':'RNA', 'GSU':'RNA', 'GTP':'RNA', 'GUA':'RNA', 'GUN':'RNA',
           'H2U':'RNA', 'HEU':'RNA', 'HPA':'RNA',   'I':'RNA',  'IC':'RNA',
            'IG':'RNA', 'IGU':'RNA', 'IMC':'RNA', 'IMP':'RNA', 'INO':'RNA',
           'IPN':'RNA',  'IU':'RNA', 'KIR':'RNA', 'L94':'RNA', 'LCA':'RNA',
           'LCG':'RNA', 'LGP':'RNA', 'LKC':'RNA', 'M1G':'RNA', 'M2G':'RNA',
           'M4C':'RNA', 'M5M':'RNA', 'MA6':'RNA', 'MAD':'RNA', 'MAU':'RNA',
           'MCY':'RNA', 'MEP':'RNA', 'MGT':'RNA', 'MIA':'RNA', 'MNU':'RNA',
           'MSP':'RNA', 'MTU':'RNA',   'N':'RNA', 'N5C':'RNA', 'N5M':'RNA',
           'N6G':'RNA', 'N79':'RNA', 'NEA':'RNA', 'O2C':'RNA', 'OMA':'RNA',
           'OMC':'RNA', 'OMG':'RNA', 'OMU':'RNA', 'ONE':'RNA',   'P':'RNA',
           'P2U':'RNA', 'P5P':'RNA', 'PGP':'RNA', 'PLR':'RNA', 'PPU':'RNA',
           'PPZ':'RNA', 'PQ0':'RNA', 'PQ1':'RNA', 'PRF':'RNA', 'PRN':'RNA',
           'PSD':'RNA', 'PSU':'RNA',  'PU':'RNA', 'PUY':'RNA', 'PYO':'RNA',
           'QSI':'RNA', 'QUO':'RNA',   'R':'RNA', 'RIA':'RNA', 'RMP':'RNA',
           'RPC':'RNA', 'RSP':'RNA', 'RSQ':'RNA', 'RUS':'RNA',   'S':'RNA',
           'S4C':'RNA', 'S4U':'RNA', 'S6G':'RNA', 'SAH':'RNA', 'SAM':'RNA',
           'SFG':'RNA', 'SLD':'RNA', 'SMP':'RNA', 'SRA':'RNA', 'SSU':'RNA',
           'SUR':'RNA',   'T':'RNA', 'T23':'RNA', 'T2T':'RNA', 'T39':'RNA',
           'T6A':'RNA', 'TEP':'RNA', 'THM':'RNA', 'THX':'RNA', 'TLC':'RNA',
           'TLN':'RNA', 'TM2':'RNA', 'TP1':'RNA',  'TS':'RNA', 'TSB':'RNA',
           'TSP':'RNA', 'TTE':'RNA',   'U':'RNA', 'U2L':'RNA', 'U2N':'RNA',
           'U34':'RNA', 'U36':'RNA', 'U37':'RNA', 'U3H':'RNA', 'U5P':'RNA',
           'U8U':'RNA', 'UAR':'RNA', 'UCL':'RNA', 'UCP':'RNA', 'UD5':'RNA',
           'UF2':'RNA', 'UFP':'RNA', 'UFT':'RNA', 'UMP':'RNA', 'UPV':'RNA',
           'UR3':'RNA', 'URA':'RNA', 'URI':'RNA', 'URU':'RNA', 'US5':'RNA',
           'UTP':'RNA', 'UZR':'RNA', 'VAA':'RNA',   'X':'RNA', 'XAN':'RNA',
           'XBB':'RNA',   'Y':'RNA', 'Y5P':'RNA',  'YG':'RNA', 'YMP':'RNA',
           'YYG':'RNA',
           '04X':'Protein', '0TD':'Protein', '2R1':'Protein', '2R3':'Protein',
           '2QY':'Protein', '2QZ':'Protein', '3GL':'Protein', '5OH':'Protein',
           'ALA':'Protein', 'ARF':'Protein', 'ARG':'Protein', 'ASN':'Protein',
           'ASP':'Protein', 'BB9':'Protein', 'CMT':'Protein', 'CSX':'Protein',
           'CYS':'Protein', 'DBU':'Protein', 'DCY':'Protein', 'DHA':'Protein',
           'DPP':'Protein', 'GLN':'Protein', 'GLU':'Protein', 'GLY':'Protein',
           'HIS':'Protein', 'HSO':'Protein', 'HYP':'Protein', 'ILE':'Protein',
           'ILX':'Protein', 'KBE':'Protein', 'LEU':'Protein', 'LLP':'Protein',
           'LYS':'Protein', 'MET':'Protein', 'MH6':'Protein', 'MHW':'Protein',
           'MLZ':'Protein', 'MSE':'Protein', 'MVA':'Protein', 'NH2':'Protein',
           'OCS':'Protein', 'PHA':'Protein', 'PHE':'Protein', 'PRO':'Protein',
           'SAR':'Protein', 'SER':'Protein', 'THR':'Protein', 'TPO':'Protein',
           'TRP':'Protein', 'TRX':'Protein', 'TS9':'Protein', 'TYR':'Protein',
           'VAL':'Protein',
           '004':'Unknown', '0EC':'Unknown', '0G6':'Unknown', '13T':'Unknown',
           '16D':'Unknown', '1F2':'Unknown', '1F3':'Unknown', '1F4':'Unknown',
           '1PE':'Unknown', '1TU':'Unknown', '20V':'Unknown', '218':'Unknown',
           '29G':'Unknown', '29H':'Unknown', '2HP':'Unknown', '2OP':'Unknown',
           '2PE':'Unknown', '2QB':'Unknown', '2QC':'Unknown', '2TB':'Unknown',
           '2ZE':'Unknown', '2ZY':'Unknown', '2ZZ':'Unknown', '34G':'Unknown',
           '365':'Unknown', '38E':'Unknown', '3AT':'Unknown', '3AW':'Unknown',
           '3DR':'Unknown', '3H3':'Unknown', '3HE':'Unknown', '3J2':'Unknown',
           '3J6':'Unknown', '3K5':'Unknown', '3K8':'Unknown', '3KD':'Unknown',
           '3KF':'Unknown', '3L2':'Unknown', '3TS':'Unknown', '42B':'Unknown',
           '574':'Unknown', '5AZ':'Unknown', '6HS':'Unknown', '6MN':'Unknown',
           '773':'Unknown', '7DG':'Unknown', '84T':'Unknown', 'A2P':'Unknown',
           'A5A':'Unknown', 'AB6':'Unknown', 'AB9':'Unknown', 'ACA':'Unknown',
           'ACE':'Unknown', 'ACP':'Unknown', 'ACT':'Unknown', 'ACY':'Unknown',
           'ADN':'Unknown', 'AF3':'Unknown', 'AG2':'Unknown', 'AGS':'Unknown',
           'AKN':'Unknown', 'ALF':'Unknown', 'AM2':'Unknown', 'ANM':'Unknown',
           'B12':'Unknown', 'B1Z':'Unknown', 'BCM':'Unknown', 'BDG':'Unknown',
           'BDR':'Unknown', 'BEF':'Unknown', 'BFT':'Unknown', 'BGC':'Unknown',
           'BMA':'Unknown', 'BME':'Unknown', 'BO4':'Unknown',  'BR':'Unknown',
           'BTN':'Unknown', 'BU1':'Unknown', 'C7P':'Unknown', 'CAC':'Unknown',
           'CAI':'Unknown', 'CIT':'Unknown',  'CL':'Unknown', 'CLM':'Unknown',
           'CLY':'Unknown', 'CMY':'Unknown', 'CNC':'Unknown', 'CNY':'Unknown',
           'CPT':'Unknown', 'CTC':'Unknown', 'CTP':'Unknown', 'CTY':'Unknown',
           'CYY':'Unknown', 'D2C':'Unknown', 'D2X':'Unknown', 'DAI':'Unknown',
           'DAL':'Unknown', 'DAR':'Unknown', 'DBB':'Unknown', 'DJF':'Unknown',
           'DIO':'Unknown', 'DMA':'Unknown', 'DOL':'Unknown', 'DPO':'Unknown',
           'DPR':'Unknown', 'DST':'Unknown', 'DTP':'Unknown', 'DUR':'Unknown',
           'EDE':'Unknown', 'EDO':'Unknown', 'EFZ':'Unknown', 'EM1':'Unknown',
           'EMK':'Unknown', 'ENX':'Unknown', 'EOH':'Unknown', 'ERN':'Unknown',
           'ERY':'Unknown',   'F':'Unknown',  'FB':'Unknown', 'FFO':'Unknown',
           'FGA':'Unknown', 'FLC':'Unknown', 'FME':'Unknown', 'FMN':'Unknown',
           'FMT':'Unknown', 'FOU':'Unknown', 'FOZ':'Unknown', 'FPD':'Unknown',
           'FUA':'Unknown', 'FUC':'Unknown', 'FUL':'Unknown', 'G19':'Unknown',
           'G34':'Unknown', 'G6P':'Unknown', 'G80':'Unknown', 'GAL':'Unknown',
           'GAU':'Unknown', 'GE1':'Unknown', 'GE2':'Unknown', 'GE3':'Unknown',
           'GET':'Unknown', 'GIR':'Unknown', 'GLA':'Unknown', 'GLP':'Unknown',
           'GND':'Unknown', 'GOL':'Unknown', 'H4B':'Unknown', 'HFA':'Unknown',
           'HMT':'Unknown', 'HRG':'Unknown', 'HYG':'Unknown', 'I2A':'Unknown',
           'IDG':'Unknown', 'IEL':'Unknown', 'ILA':'Unknown', 'IOD':'Unknown',
           'IPA':'Unknown', 'IPH':'Unknown', 'IRI':'Unknown', 'ISH':'Unknown',
           'ISI':'Unknown', 'IUM':'Unknown', 'JOS':'Unknown', 'JS4':'Unknown',
           'JS5':'Unknown', 'JS6':'Unknown', 'KAN':'Unknown', 'KSG':'Unknown',
           'L8H':'Unknown', 'LC2':'Unknown', 'LHA':'Unknown', 'LIV':'Unknown',
           'LLL':'Unknown', 'LMA':'Unknown', 'LMS':'Unknown', 'LYA':'Unknown',
           'M2M':'Unknown', 'M5Z':'Unknown', 'MAN':'Unknown', 'MEA':'Unknown',
           'MES':'Unknown', 'MGR':'Unknown', 'MHT':'Unknown', 'MHU':'Unknown',
           'MHV':'Unknown', 'MIX':'Unknown', 'MLI':'Unknown', 'MMC':'Unknown',
           'MPD':'Unknown', 'MRC':'Unknown', 'MRD':'Unknown', 'MT9':'Unknown',
           'MTT':'Unknown', 'MUB':'Unknown', 'MUL':'Unknown', 'MYL':'Unknown',
           'MYN':'Unknown', 'N30':'Unknown', 'N33':'Unknown', 'NAG':'Unknown',
           'NCO':'Unknown', 'NDG':'Unknown', 'NEB':'Unknown', 'NEG':'Unknown',
           'NF2':'Unknown', 'NH4':'Unknown', 'NHE':'Unknown', 'NME':'Unknown',
           'NMY':'Unknown', 'NO1':'Unknown', 'NO3':'Unknown', 'NPM':'Unknown',
           'NVA':'Unknown', 'NVP':'Unknown', 'OHX':'Unknown', 'OLZ':'Unknown',
           'ON0':'Unknown', 'P12':'Unknown', 'P13':'Unknown', 'P14':'Unknown',
           'P1P':'Unknown', 'P24':'Unknown', 'P6G':'Unknown', 'PA1':'Unknown',
           'PA2':'Unknown', 'PA3':'Unknown', 'PAE':'Unknown', 'PAR':'Unknown',
           'PCY':'Unknown', 'PDI':'Unknown', 'PEG':'Unknown', 'PEV':'Unknown',
           'PG4':'Unknown', 'PGE':'Unknown', 'PGV':'Unknown', 'PMZ':'Unknown',
           'PO2':'Unknown', 'PO4':'Unknown', 'POP':'Unknown', 'PPV':'Unknown',
           'PRI':'Unknown', 'PRL':'Unknown', 'PUP':'Unknown', 'PYI':'Unknown',
           'PYY':'Unknown', 'QUA':'Unknown', 'R14':'Unknown', 'RAP':'Unknown',
           'RBF':'Unknown', 'RHD':'Unknown', 'RIB':'Unknown', 'RIO':'Unknown',
           'ROS':'Unknown', 'ROX':'Unknown', 'RPO':'Unknown', 'RS3':'Unknown',
           'RTP':'Unknown', 'RU6':'Unknown', 'S9L':'Unknown', 'SCM':'Unknown',
           'SCY':'Unknown', 'SE4':'Unknown', 'SEP':'Unknown', 'SF4':'Unknown',
           'SIN':'Unknown', 'SIS':'Unknown', 'SJP':'Unknown', 'SLZ':'Unknown',
           'SO4':'Unknown', 'SPD':'Unknown', 'SPK':'Unknown', 'SPM':'Unknown',
           'SPR':'Unknown', 'SPS':'Unknown', 'SRY':'Unknown', 'SS0':'Unknown',
           'SSA':'Unknown', 'STD':'Unknown', 'SUC':'Unknown', 'SVN':'Unknown',
           'T17':'Unknown', 'T1C':'Unknown', 'T8B':'Unknown', 'TAC':'Unknown',
           'TAO':'Unknown', 'TEL':'Unknown', 'THF':'Unknown', 'TOA':'Unknown',
           'TOB':'Unknown', 'TOC':'Unknown', 'TOY':'Unknown', 'TPP':'Unknown',
           'TPS':'Unknown', 'TRS':'Unknown', 'TS6':'Unknown', 'TSE':'Unknown',
           'TYE':'Unknown', 'TYK':'Unknown', 'UAL':'Unknown', 'UAM':'Unknown',
           'UBD':'Unknown', 'UDP':'Unknown', 'UNK':'Unknown', 'UNL':'Unknown',
           'UNX':'Unknown', 'VIB':'Unknown', 'VIF':'Unknown', 'VIR':'Unknown',
           'VO4':'Unknown', 'WIN':'Unknown', 'WO2':'Unknown', 'XXX':'Unknown',
           'ZBA':'Unknown', 'ZIT':'Unknown', 'ZLD':'Unknown', 'ZUK':'Unknown',
           '3CO':'Metal', 'AU3':'Metal',  'BA':'Metal',  'CA':'Metal',
            'CD':'Metal',  'CO':'Metal',  'CS':'Metal', 'EU3':'Metal',
           'FE2':'Metal',  'HG':'Metal',  'IR':'Metal', 'IR3':'Metal',
             'K':'Metal',  'LU':'Metal',  'MG':'Metal',  'MN':'Metal',
            'NA':'Metal',  'NI':'Metal',  'OS':'Metal',  'PB':'Metal',
            'PT':'Metal', 'PT4':'Metal',  'RB':'Metal',  'RU':'Metal',
            'SR':'Metal',  'TB':'Metal',  'TL':'Metal',  'ZN':'Metal',
           'HOH':'Water'} 

def CIFrestype(params):

    typ,name = params

    Metals = {'CADMIUM ION'       :1, 'POTASSIUM ION'     :1, 'MAGNESIUM ION'   :1, 'ZINC ION'         : 1,
              'STRONTIUM ION'     :1, 'SODIUM ION'        :1, 'CALCIUM ION'     :1, 'MANGANESE (II) ION':1,
              'LUTETIUM (III) ION':1, 'EUROPIUM (III) ION':1, 'NICKEL (II) ION' :1, 'BARIUM ION'        :1,
              'CESIUM ION'        :1, 'THALLIUM (I) ION'  :1, 'MERCURY (II) ION':1, 'OSMIUM ION'        :1,
              'GOLD 3+ ION'       :1, 'COBALT (II) ION'   :1, 'COBALT (III) ION':1, 'PLATINUM (II) ION' :1,
              'IRIDIUM (III) ION' :1, 'FE (II) ION'       :1, 'RUBIDIUM ION'    :1, 'IRIDIUM ION'       :1,
              'LEAD (II) ION'     :1, 'TERBIUM(III) ION'  :1, 'RUTHENIUM ION'   :1, 'PLATINUM (IV) ION' :1}

    Proteins = {'PEPTIDE LINKING':1, 'L-PEPTIDE LINKING':1,'D-PEPTIDE LINKING':1,
                'D-GAMMA-PEPTIDE, C-DELTA LINKING':1, 'L-PEPTIDE NH3 AMINO TERMINUS':1}

    RNAS  = {'L-RNA LINKING':1, 'RNA LINKING':1,'RNA OH 3 PRIME TERMINUS':1}

    if typ.upper() in ('NON-POLYMER','SACCHARIDE','PEPTIDE-LIKE','D-SACCHARIDE','L-SACCHARIDE'):
        if name in Metals:  return 'Metal'
        elif name=='WATER': return 'Water'
        else:               return 'Unknown'
    else:
        if   typ.upper() in  RNAS      : return 'RNA'
        elif typ.upper()=='DNA LINKING': return 'DNA'
        elif typ.upper() in Proteins   : return 'Protein'
        else:
            print('"%s" is inappropriate type of residue'%typ)
            return 'Unknown'


def check_na(dist): # check for n/a in distances and torsions

    if 'n/a' in dist: return '\\N'
    else            : return float(dist)

def clean_seq(value): # delete double spaces and then delete spaces from start and end of line

    value = cut_double_spaces(value)
    value = value.replace('- ','-')
    if value[0] == ' ' : value = value[1:]
    if value[-1]== ' ': value = value[:-1]

    return value

def clean_source(value): # delete double spaces and then delete space or ";" from end of line

    value = cut_double_spaces(value)
    while value[-1] == ' ': value = value[:-1]
    while value[0] == ' ': value = value[1:]
    if value[-1] == ';': value = value[:-1]

    return value

def cut_double_spaces(string) -> 'String with no double spaces':

    string2 = ''
    last_is_space = 0

    for i in range(len(string)):

        if not(string[i] == ' ' and last_is_space): string2 += string[i]

        if string[i] == ' ': last_is_space = 1
        else               : last_is_space = 0

    return string2

def cut_for_int(string) -> 'String is ready for int()':

    string2 = ''

    for char in string:

        if char in '1234567890': string2 += char

    return string2

def cut_for_float(string) -> 'String is ready for float()':

    string2 = ''

    for char in string:

        if char in '1234567890 .': string2 += char

    return string2

def cut_spaces(string) -> 'String with no spaces':

    string2 = ''

    for char in string:

        if char != ' ': string2 += char

    return string2

def pdbnum(number) -> "'int' + char (letter or space)":

    number = cut_spaces(number)
    if number[-1] in '1234567890': number += ' '

    return number

def pdbnum_to_float(number) -> 'like "128A" -> 128.03':

    number = pdbnum(number)

    return float(number[:-1]) + letter_weight[number[-1]]

def Atom(string):

    return {'ID'     :                         0,       # ID of atom
            'NUM'    :         int(string[6:11]),       # number from PDB file
            'NAME'   : cut_spaces(string[12:16]),       # Name of atom in residue
            'ALTLOC' :                string[16],       # Alternate location
            'RESNAME': cut_spaces(string[17:20]),       # Name of residue
            'CHAIN'  :                string[21],       # letter identifier of chain
            'TYPE'   :                        '',       # RNA, Protein, Unknown, Metal or Water
            'RESNUM' :     pdbnum(string[22:27]),       # PDBNUM of residue
            'X'      :      float(string[30:38]),       # x coordinate
            'Y'      :      float(string[38:46]),       # y coordinate
            'Z'      :      float(string[46:54]),       # z coordinate
            'OCCUP'  :      float(string[54:60]),       # Occupancy
            'ELEM'   : cut_spaces(string[76:78]),       # Chemical element
            'BONDS'  :                         0}

def AtomCIF(row):

    name = row[3]
    if name[-1] == '"': name = name[:-1]
    if name[0]  == '"': name = name[1:]

    cifid = row[8]
    if cifid == '.': cifid = '\\N'
    else: cifid = int(cifid)

    if len(row)<23:
        chch = row[18]
        resnum = row[16]+row[9]
    else:
        chch = row[23]
        resnum = row[21]+row[9]
    if resnum[-1]=='?': resnum = resnum[:-1]+' '
    

    return {'ID'     :                       0,       # ID of atom
            'NUM'    :             int(row[1]),       # number from PDB file
            'NAME'   :                    name,       # Name of atom in residue
            'ALTLOC' : row[4].replace('.',' '),       # Alternate location
            'RESNAME':                  row[5],       # Name of residue
            'CHAIN'  :                    chch,       # letter identifier of chain
            'TYPE'   :                      '',       # RNA, Protein, Unknown, Metal or Water
            'RESNUM' :                  resnum,       # PDBNUM of residue
            'X'      :          float(row[10]),       # x coordinate
            'Y'      :          float(row[11]),       # y coordinate
            'Z'      :          float(row[12]),       # z coordinate
            'OCCUP'  :          float(row[13]),       # Occupancy
            'TEMPF'  :                 row[14],
            'ELEM'   :                  row[2],       # Chemical element
            'BONDS'  :                       0,
            'CIFID'  :                   cifid}       # CIFID of residue

def BP(lines): # parsing from DSSR

    num   = int(lines[0][:5])
    first = cut_double_spaces(lines[0][5:-1]).split(' ')
    
    if first[0] == '': first = first[1:] # if bp_num >= 10000 (one char shift)

    if len(first) < 7: first.insert(3,'n')           # if there is no bondtype
    if len(first) > 7: first = first[:3] + [first[3]+' '+first[4],] + first[5:] # if rWC (Levitt)

    nucl1       = '.'.join(first[0].split('.')[2:])  # chain.nucl.pdbnum.icode
    nucl2       = '.'.join(first[1].split('.')[2:])
    bond        = first[2]
    bondtype    = bptype_dict[first[3]]
    bondclasses = first[4:]
    
    chis    = lines[1][8:-2].split('] [')

    third   = lines[2][7:-1].split(' ')
    dists   = third[:-1]
    torsion = third[-1]
    c1c1    = check_na(dists[0][11:])
    n1n9    = check_na(dists[1][9:])
    c6c8    = check_na(dists[2][9:])
    tor     = check_na(torsion[19:])

    fourth    = lines[3][7:]
    colon     = fourth.find(':')
    hbondsnum = int(fourth[8:colon-1])
    hbonds    = fourth[colon+3:-2].split(',')

    fifth  = lines[5][16:-1]
    fifth2 = cut_double_spaces(fifth[1:-1]).split(' ')

    return {'ID'       :                 num,
            'NUCL1'    :               nucl1,
            'NUCL2'    :               nucl2,
            'PAIR'     : nucl1 + '-' + nucl2,
            'BOND'     :                bond,
            'TYPE'     :            bondtype,
            'CLASS'    :         bondclasses,
            'CHAIN1'   : nucl1.split('.')[0],
            'CHAIN2'   : nucl2.split('.')[0],
            'INFO1'    :             chis[0],
            'INFO2'    :             chis[1],
            'DIST1'    :                c1c1,
            'DIST2'    :                n1n9,
            'DIST3'    :                c6c8,
            'TOR'      :                 tor,
            'HBONDSNUM':           hbondsnum,
            'HBONDS'   :              hbonds,
            'PARAMS'   :               fifth,
            'SHEAR'    :           fifth2[0],
            'STRETCH'  :           fifth2[1],
            'STAGGER'  :           fifth2[2],
            'BUCKLE'   :           fifth2[3],
            'PROPELLER':           fifth2[4],
            'OPENING'  :           fifth2[5],
            'STEM'     :                None,
            'OLDSTEM'  :                None,
            'FULLSTEM' :                None,
            'REVSTEM'  :                None,
            'LUSTEM'   :                None,
            'LINK'     :                None,
            'HELIX'    :                None,
            'NUCLMULT' :                None,
            'STEP'     :                None}
            
def Chain(letter):

    return {'ID'         :      0,   # ID of chain
            'START'      :  '\\N',   # int(pdbnum) of first residue
            'END'        :  '\\N',   # int(pdbnum) of last residue
            'SEQ'        :     '',   # sequence list from SEQRES
            'SEQ2'       :     '',   # real sequence list
            'LIGSEQ'     :     [],   # sequence of ligands
            'LENGTH'     :  '\\N',   # length of chain
            'MISRES'     :     [],   # missing residues list
            'RES'        :     [],   # residues list
            'NUMS'       :     {},   # dict[pdbnum] = index of 'RES'
            'TYPE'       :  '\\N',   # RNA, DNA, Protein or Unknown
            'CHID'       : letter,   # letter identifier
            'MOL_ID'     :  '\\N',   # ID of molecule
            'LIGANDS'    :     [],   # ligands list
            'GARBAGE'    :  False,   # garbage chain or not
            'LUBRACKETS' :  '\\N',
            'BRACKETS'   :  '\\N',   # stems only
            'SLBRACKETS' :  '\\N',   # including wc/wb links
            'SCHEME'     :  '\\N',
            'SIGN'       :  '\\N',
            'BOUND'      :      0,
            'BOUNDWCWB'  :      0,
            'LIGBOUND'   :      0,
            'DIAGRAM'    :  '\\N'}

def Molecule(Mol_id):
                                        # See PDB Format Documentation
    return {'ID'                                  : Mol_id,       
            'MOLECULE'                            :     '',
            'CHAIN'                               :     '',
            'FRAGMENT'                            :     '',
            'SYNONYM'                             :     '',
            'EC'                                  :     '',
            'ENGINEERED'                          :     '',
            'MUTATION'                            :     '',
            'OTHER_DETAILS'                       :     '',
            'SYNTHETIC'                           :     '',
            'ORGANISM_SCIENTIFIC'                 :     '',
            'ORGANISM_COMMON'                     :     '',
            'ORGANISM_TAXID'                      :     '',
            'STRAIN'                              :     '',
            'VARIANT'                             :     '',
            'CELL_LINE'                           :     '',
            'ATCC'                                :     '',
            'ORGAN'                               :     '',
            'TISSUE'                              :     '',
            'CELL'                                :     '',
            'ORGANELLE'                           :     '',
            'SECRETION'                           :     '',
            'CELLULAR_LOCATION'                   :     '',
            'PLASMID'                             :     '',
            'GENE'                                :     '',
            'EXPRESSION_SYSTEM'                  :     '',
            'EXPRESSION_SYSTEM_COMMON'           :     '',
            'EXPRESSION_SYSTEM_TAXID'            :     '',
            'EXPRESSION_SYSTEM_STRAIN'           :     '',
            'EXPRESSION_SYSTEM_VARIANT'          :     '',
            'EXPRESSION_SYSTEM_CELL_LINE'        :     '',
            'EXPRESSION_SYSTEM_ATCC_NUMBER'      :     '',
            'EXPRESSION_SYSTEM_ORGAN'            :     '',
            'EXPRESSION_SYSTEM_TISSUE'           :     '',
            'EXPRESSION_SYSTEM_CELL'             :     '',
            'EXPRESSION_SYSTEM_ORGANELLE'        :     '',
            'EXPRESSION_SYSTEM_CELLULAR_LOCATION':     '',
            'EXPRESSION_SYSTEM_VECTOR_TYPE'      :     '',
            'EXPRESSION_SYSTEM_VECTOR'           :     '',
            'EXPRESSION_SYSTEM_PLASMID'          :     '',
            'EXPRESSION_SYSTEM_GENE'             :     ''}

def MisResidue():

    return {'ID'        :     0,   # ID of residue
            'TYPE'      :    '',   # RNA, DNA, Protein, Unknown, Metal or Water
            'PDBNUM'    :    '',   # Number + icode from PDB file
            'NAME'      :    '',   # Name of residue
            'CHAIN'     :    '',   # letter identifier of chain
            'ATOMS'     :    [],   # atoms list
            'FLOAT'     :     0,   # float equivalent of 'PDBNUM'
            'MISS'      :  True,
            'WING'      :  None,
            'OLDWING'   :  None,
            'FSTEMS'    :     0,
            'THREAD'    :  None,
            'ZIP'       : '\\N',
            'LUMULT'    : '\\N',
            'MULT'      : '\\N',
            'BPS'       :     0,
            'BRACKETS'  :   '-',   # stems only
            'SLBRACKETS':   '-'}   # including wc/wb links

def Residue(atom):

    icode = atom['RESNUM'][-1]
    if icode == ' ': icode = ''
    num   = int(atom['RESNUM'][:-1])

    lu    = atom['CHAIN'] + '.' + atom['RESNAME'] + '.' + str(num) + '.' + icode 
    Float = pdbnum_to_float(atom['RESNUM'])
        
    return {'ID'            :               0,   # ID of residue
            'TYPE'          :              '',   # RNA, DNA, Protein, Unknown, Metal or Water
            'PDBNUM'        :  atom['RESNUM'],   # Number + icode from PDB file
            'NAME'          : atom['RESNAME'],   # Name of residue
            'CHAIN'         :   atom['CHAIN'],   # letter identifier of chain
            'ATOMS'         :              [],   # atoms list
            'FLOAT'         :           Float,   # float equivalent of 'PDBNUM'
            'DSSR'          :              lu,   # chain.nucl.pdbnum.icode (like DSSR)
            'MISS'          :           False,
            'WING'          :            None,
            'OLDWING'       :            None,
            'FSTEMS'        :               0,
            'THREAD'        :            None,
            'ZIP'           :           '\\N',
            'LUMULT'        :           '\\N',
            'MULT'          :           '\\N',
            'BPS'           :               0,
            'BRACKETS'      :             '.',   # stems only
            'SLBRACKETS'    :             '.'}   # including wc/wb links

def ResidueCIF(line,inds): # Residue from mmcif._pdbx_poly_seq_scheme

    icode = line[inds[10]]
    if icode == '.': icode = ''

    lu    = line[inds[9]] + '.' + line[inds[3]] + '.' + line[inds[5]] + '.' + icode
    if not icode: icode = ' '
    Float = pdbnum_to_float(line[inds[5]]+icode)

    return {'ID'        :                     0,   # ID of residue
            'CIFID'     :    int(line[inds[2]]),   # seq_id from mmCIF
            'TYPE'      :                    '',   # RNA, DNA, Protein, Unknown, Metal or Water
            'PDBNUM'    :   line[inds[5]]+icode,   # Number + icode from PDB file
            'NAME'      :         line[inds[3]],   # Name of residue
            'CHAIN'     :         line[inds[9]],   # letter identifier of chain
            'ATOMS'     :                    [],   # atoms list
            'FLOAT'     :                 Float,   # float equivalent of 'PDBNUM'
            'DSSR'      :                    lu,   # chain.nucl.pdbnum.icode (like DSSR)
            'MISS'      :                 False,
            'WING'      :                  None,
            'OLDWING'   :                  None,
            'FSTEMS'    :                     0,
            'THREAD'    :                  None,
            'ZIP'       :                 '\\N',
            'LUMULT'    :                 '\\N',
            'MULT'      :                 '\\N',
            'BPS'       :                     0,
            'BRACKETS'  :                   '.',   # stems only
            'SLBRACKETS':                   '.'}   # including wc/wb links

def ResidueCIFAtom(atom): # Residue from mmcif._pdbx_poly_seq_scheme

    Float = pdbnum_to_float(atom['RESNUM'])
    lu    = atom['CHAIN'] + '.' + atom['RESNAME'] + '.' + atom['RESNUM'][:-1] + '.' + atom['RESNUM'][-1].replace(' ','')

    return {'ID'        :               0,   # ID of residue
            'CIFID'     :   atom['CIFID'],   # seq_id from mmCIF
            'TYPE'      :    atom['TYPE'],   # RNA, DNA, Protein, Unknown, Metal or Water
            'PDBNUM'    :  atom['RESNUM'],   # Number + icode from PDB file
            'NAME'      : atom['RESNAME'],   # Name of residue
            'CHAIN'     :   atom['CHAIN'],   # letter identifier of chain
            'ATOMS'     :              [],   # atoms list
            'FLOAT'     :           Float,   # float equivalent of 'PDBNUM'
            'DSSR'      :              lu,   # chain.nucl.pdbnum.icode (like DSSR)
            'MISS'      :           False,
            'WING'      :            None,
            'OLDWING'   :            None,
            'FSTEMS'    :               0,
            'THREAD'    :            None,
            'ZIP'       :           '\\N',
            'LUMULT'    :           '\\N',
            'MULT'      :           '\\N',
            'BPS'       :               0,
            'BRACKETS'  :             '.',   # stems only
            'SLBRACKETS':             '.'}   # including wc/wb links

def LuA_minor(lines): # parsing from DSSR

    num     = int(lines[0][:5])
    first   = cut_double_spaces(lines[0][5:-1]).split(' ')
    if not first[0]: first = first[1:]

    first   = first[:1] + first[1].split('\t') + first[2:]
    first   = first[:2] + first[2].split('|')  + first[3:]
    typ     = first[0][5:]
    look    = first[1]
    nucl    = '.'.join(first[2].split('.')[2:])
    pair    = first[3].split(',')
    pair    = '.'.join(pair[0].split('.')[2:]) + '.' + '.'.join(pair[1].split('.')[2:])
    bonds1  = lines[1].split(':')[1][2:-2].split(',')
    if not bonds1[0]: bonds1 = []
    bonds2  = lines[2].split(':')[1][2:-2].split(',')
    if not bonds2[0]: bonds2 = []

    return {'ID'    :    num,
            'TYPE'  :    typ,
            'CLASS' :   look,
            'NUCL'  :   nucl, 
            'PAIR'  :   pair,
            'BONDS1': bonds1,
            'BONDS2': bonds2}

def LuAtomBaseCap(line):

    cap = line.strip().split()
    cap[2] = cap[2].split('@')


    return {'ID'       : int(cap[0]),
            'TYPE'     : cap[1],
            'ATOM'     : cap[2][0],
            'ATOMBASE' : '.'.join(cap[2][1].split('.')[2:]),
            'BASE'     : '.'.join(cap[3].split('.')[2:]),
            'RISE'     : float(cap[4]),}

def LuBrackets(lines): # parsing from DSSR

    first  = cut_double_spaces(lines[0][:-1]).split(' ')

    ch = first[0].split('-')[-1]

    if 'whole' in lines[0]:

        chain = 'ALL'
        size  = int(first[1][4:].replace('*',''))
        Break = True

    else:

        chain  = ch
        size   = int(first[2][4:])
        Break  = '*' in lines[0]

    seq    = lines[1][:-1]
    dia    = lines[2][:-1]

    return {'CHAIN'  :  chain,
            'LENGTH' :   size,
            'BREAK'  :  Break,
            'SEQ'    :    seq,
            'DIAGRAM':    dia}

def LuBulge(lines): # parsing from DSSR

    #print(lines)##
    num   = int(lines[0][:5])

    first = cut_double_spaces(lines[0][5:-1]).split(' ')

    size  = int(first[1][4:-1]) - 4
    look  = first[2][1:-2]
    left  = int(first[5][2:first[5].find(',')])
    right = int(first[5][first[5].rfind('#')+1:-1])

    i     = 3
    if look.startswith('0,'): i = 4

    line = cut_double_spaces(lines[i][:-1]).split(' ')

    seq   = line[2]
    nucls = line[3].split(',')
    for i in range(len(nucls)): nucls[i] = '.'.join(nucls[i].split('.')[2:])

    return {'ID'    :   num,
            'LENGTH':  size,
            'TYPE'  :  look,
            'PREV'  :  left,
            'NEXT'  : right,
            'SEQ'   :   seq,
            'NUCLS' : nucls}

def LuHairpin(lines): # parsing from DSSR

    first  = cut_double_spaces(lines[0][:-1]).split(' ')[1:]
    second = cut_double_spaces(lines[2][:-1]).split(' ')[1:]
    num    = int(first[0])
    first  = first[1:]

    nts   = int(first[3][1:-2])
    seq   = second[1][1:-1]

    closing = int(first[6][2:-1])
    
    second = second[2]
    second = second.split(',')
    for i in range(len(second)): second[i] = '.'.join(second[i].split('.')[2:])
    nucls = second[1:-1]

    return {'ID'     :     num,
            'LENGTH' :     nts,
            'SEQ'    :     seq,
            'CLOSING': closing,
            'NUCLS'  :   nucls}

def LuHelix(lines): # parsing from DSSR

    steps = {}

    first = lines[0][2:].split(' ')
    first[0] = first[0].replace('*','')
    num = int(first[0][6:first[0].find('[')])
    bps = int(first[1][4:])
    strand1    = lines[1][6:].split(' ')[1][:-1]
    bp_type    = cut_double_spaces(lines[2][7:]).split(' ')[1][:-1]
    strand2    = lines[3][6:].split(' ')[1][:-1]
    helix_form = cut_double_spaces(lines[4][6:]).split(' ')[1][:-1]

    pairs = []

    for line in lines[10:]:

        if line[3] != ' ':

            pair_info = cut_double_spaces(line[5:]).split(' ')
            nucls = [pair_info[0].split('.')[2:],pair_info[1].split('.')[2:]]
            pairs.append('.'.join(nucls[0]) + '.' + '.'.join(nucls[1]))

        elif 'step-pars' in line: steps[pairs[-1]] = cut_double_spaces(line.split('[')[1])[:-2]

    return {'ID'    :        num,
            'SIZE'  :        bps,
            'SEQ1'  :    strand1,
            'SEQ2'  :    strand2,
            'BPTYPE':    bp_type,
            'FORM'  : helix_form,
            'PAIRS' :      pairs,
            'STEMS' :          0,
            'LONES' :          0},steps

def LuInternal(lines): # parsing from DSSR

    num = int(lines[0][:5].replace('x',''))

    first = cut_double_spaces(lines[0][5:-1]).split(' ')
    if not first[0]: first = first[1:]
    
    sym   = first[0][0] == 's'

    size  = int(first[3][4:-1]) - 4
    look  = first[4][1:-2]
    left  = int(first[7][2:first[7].find(',')])
    right = int(first[7][first[7].rfind('#')+1:-1])

    line1 = cut_double_spaces(lines[3]).split(' ')
    line2 = cut_double_spaces(lines[4]).split(' ')

    seq1   = line1[2]
    seq2   = line2[2]
    
    lnucls = line1[3][:-1].split(',')
    rnucls = line2[3][:-1].split(',')
    for i in range(len(lnucls)): lnucls[i] = '.'.join(lnucls[i].split('.')[2:])
    for i in range(len(rnucls)): rnucls[i] = '.'.join(rnucls[i].split('.')[2:])

    return {'ID'    :    num,
            'LENGTH':   size,
            'TYPE'  :   look,
            'SYM'   :    sym,
            'PREV'  :   left,
            'NEXT'  :  right,
            'LSEQ'  :   seq1,
            'RSEQ'  :   seq2,
            'LNUCLS': lnucls,
            'RNUCLS': rnucls}

def LuJunction(lines): # parsing from DSSR
    
    num = int(lines[0][:5].replace('*','').replace('x',''))

    first = cut_double_spaces(lines[0][5:-1]).split(' ')
    if not first[0]: first = first[1:]

    threadsnum = int(first[0][:-4])
    size       = int(first[2][4:-1]) - threadsnum*2
    look       = first[3][1:-2]

    closing = first[6][1:-1].split(',')
    for i in range(len(closing)): closing[i] = int(closing[i][1:])

    seqs   = []
    nuclss = []

    for line in lines[3:]:

        line      = cut_double_spaces(line[:-1]).split(' ')
        length    = int(line[1][4:])
        seq,nucls = '', []

        if length:

            seq   = line[2]
            nucls = line[3].split(',')
            for i in range(len(nucls)): nucls[i] = '.'.join(nucls[i].split('.')[2:])

        seqs.append(seq)
        nuclss.append(nucls)

    return {'ID'     :        num,
            'LENGTH' :       size,
            'TYPE'   :       look,
            'THREADS': threadsnum,
            'CLOSING':    closing,
            'SEQS'   :       seqs,
            'NUCLS'  :     nuclss}

def LuK_turn(lines): # parsing from DSSR

    num = int(lines[0][:5])

    neither = 'neither' in lines[0]
    w       = '[w]'     in lines[0]

    first = cut_double_spaces(lines[0][5:-1]).split(' ')

    Type = {'Normal':0,'REVERSE':1,'Undecided':2, '[w]':3}[first[0]]

    helix = '\\N'

    iloop = int(first[-5][first[-5].find('#')+1:-1])

    second = cut_double_spaces(lines[1][6:-1]).split(' ')

    stem1 = int(second[0][second[0].find('#')+1:second[0].find('[')])

    if w: stem2 = int(second[3][second[3].find('#')+1:second[3].find('[')])
    else: stem2 = int(second[4][second[4].find('#')+1:second[4].find('[')])

    if w: pair = '\\N'
    else:
        pair = second[3][:-1].split(',')
        pair = '.'.join(pair[0].split('.')[2:]) + '.' + '.'.join(pair[1].split('.')[2:])

    strand1 = lines[2].strip().split()[3]
    strand2 = lines[3].strip().split()[3]

    return {'ID'   :   num,
            'PAIR' :  pair,
            'HELIX': helix,
            'STEM1': stem1,
            'STEM2': stem2,
            'ILOOP': iloop,
            'TYPE' :  Type,
            'STRAND1': strand1,
            'STRAND2': strand2}

def LuKissing(line): # parsing from DSSR

    num   = int(line[:5])
    line  = cut_double_spaces(line[5:-1]).split(' ')
    kiss  = int(line[1][1:])
    loop1 = int(line[5][1:])
    loop2 = int(line[7][1:])

    return {'ID'      :   num,
            'KISS'    :  kiss,
            'HAIRPIN1': loop1,
            'HAIRPIN2': loop2}
    
def LuLone(line): # parsing from DSSR

    line = cut_double_spaces(line).split(' ')

    helix = None
    otherhel = []

    if line[0] == 'n/a':
        pass        
    elif ',' not in line[0]:
        helix = int(line[0][2:-1])
    else:
        helices = line[0][1:-1].split(',')
        helix = int(helices[0][1:])
        otherhel = [int(i[1:]) for i in helices[1:]]
    
    num  = int(line[1])
    line = line[2:]

    bps = 1

    nucls = [line[0].split('.')[2:],line[1].split('.')[2:]]

    strand1    = "5'-"+nucls[0][1]+"-3'"
    strand2    = "3'-"+nucls[1][1]+"-5'"
    helix_form = '.'

    pair = '.'.join(nucls[0])+'.'+'.'.join(nucls[1])
    
    return {'ID'          :        num,
            'DSSRID'      :        num,
            'HELIX'       :      helix,
            'OTHERHELICES':   otherhel,
            'LENGTH'      :        bps,
            'SEQ1'        :    strand1,
            'SEQ2'        :    strand2,
            'FORM'        : helix_form,
            'PAIRS'       :    [pair,]}                                

def LuMult(line): # parsing from DSSR

    while line[0] == ' ': line = line[1:]
    while '  ' in line: line = line.replace('  ',' ')

    line  = line[:-1].split(' ')

    num   = int(line[0])
    size  = int(line[1][4:].replace('*',''))
    nucls = line[3].split(',')
    plan  = float(line[4][10:])

    for i in range(len(nucls)): nucls[i] = '.'.join(nucls[i].split('.')[2:])

    seq   = line[2]

    return {'ID'        :   num,
            'SIZE'      :  size,
            'NUCLS'     : nucls,
            'SEQ'       :   seq,
            'PLANARITY' :  plan}
    
def LuNon_loop(line): # parsing from DSSR

    num   = int(line[:5])
    line  = cut_double_spaces(line[5:-1]).split(' ')
    Break = False
    if '*' in line[0]: Break = True
    size  = int(line[0][4:].replace('*',''))
    seq   = ''
    if size: seq = line[1]
    start = ''
    end   = ''

    if   size == 1:  start = end = '.'.join(line[2].split('.')[2:])

    elif size != 0:

        nucls = line[2].split(',')
        start = '.'.join(nucls[0].split('.')[2:])
        end   = '.'.join(nucls[-1].split('.')[2:])

    return {'ID'    :   num,
            'BREAK' : Break,
            'LENGTH':  size,
            'SEQ'   :   seq,
            'START' : start,
            'END'   :   end}

def LuNon_pair(line): # parsing from DSSR

    num = int(line[:5])
    #print(line)
    line = cut_double_spaces(line[5:-1]).split(' ')
    if not line[0]: line = line[1:]

    nucl1 = '.'.join(line[0].split('.')[2:])
    nucl2 = '.'.join(line[1].split('.')[2:])
    #print(nucl1,nucl2)
    stack, hb = 0, 0

    if line[2] == 'stacking:':
        stack = 2
        if len(line)>5 and line[5][:7]=='H-bonds': hb = 5

    elif line[2][:9]=='interBase':
        hb    = 3
    else:
        hb    = 2

    if stack: stacking = line[stack+1]
    else    : stacking = '\\N'

    if hb:
        hbondsnum = int(line[hb][8:-2])
        hbonds    = line[hb+1][1:-1]
    else:
        hbondsnum = 0
        hbonds    = '\\N'

    overlap = line[2]

    mindist = ''
    angle = ''

    for k in line:
        if k.startswith('interBase-angle'):
            angle = k.split('=')[1]
        if k.startswith('min_baseDist'):
            mindist = k.split('=')[1]

    return {'ID'       :       num,
            'NUCL1'    :     nucl1,
            'NUCL2'    :     nucl2,
            'STACKING' :  stacking,
            'HBONDSNUM': hbondsnum,
            'HBONDS'   :    hbonds,
            'MINDIST'  :   mindist,
            'ANGLE'    :     angle}

def LuStem(lines): # parsing from DSSR

    steps = {}

    first    = lines[0][2:].split(' ')
    nums     = first[0].split('[')
    stem_num = int(nums[0][5:])

    if nums[1]!=']': helix = int(nums[1][1:-1])
    else           : helix = None

    parallel = first[-1][:-1] == 'parallel'

    other = [] # if more than one helix!

    if parallel: flag = -2
    else       : flag = -1

    bps = int(first[flag][4:])
    for line in first[1:flag]: other.append(int(line.replace('#','').replace(']*','').replace(',','')))

    strand1    = lines[1][6:].split(' ')[1][:-1]
    strand2    = lines[3][6:].split(' ')[1][:-1]
    helix_form = cut_double_spaces(lines[4][6:]).split(' ')[1][:-1]

    pairs = []

    for line in lines[10:]:

        if line[3] != ' ':

            pair_info = cut_double_spaces(line[5:]).split(' ')
            nucls = [pair_info[0].split('.')[2:],pair_info[1].split('.')[2:]]
            pairs.append('.'.join(nucls[0]) + '.' + '.'.join(nucls[1]))

        elif 'step-pars' in line: steps[pairs[-1]] = cut_double_spaces(line.split('[')[1])[:-2]

    return {'ID'          :   stem_num,
            'DSSRID'      :   stem_num,
            'HELIX'       :      helix,
            'OTHERHELICES':      other,
            'LENGTH'      :        bps,
            'SEQ1'        :    strand1,
            'SEQ2'        :    strand2,
            'FORM'        : helix_form,
            'PAIRS'       :      pairs,
            'PARALLEL'    :   parallel},steps

def split_by_minus(line: str):    

    firstnum = int(line.split('.')[4])
    if firstnum >= 0:
        br = line.find('-')
    else:
        br = line.find('-',line.find(str(firstnum))+1)

    return [line[:br],line[br+1:]]

def LuU_turn(line): # parsing from DSSR
    #print(line)
    num   = int(line[:5].replace('*',''))
    line  = cut_double_spaces(line[5:-1]).split(' ')
    if not line[0]: line = line[1:]
    nucls = split_by_minus(line[0])
    nucl1 = '.'.join(nucls[0].split('.')[2:])
    nucl2 = '.'.join(nucls[1].split('.')[2:])
    bonds = line[2].split(',')
    if not bonds[0]: bonds = []

    return {'ID'   :   num,
            'NUCL1': nucl1,
            'NUCL2': nucl2,
            'BONDS': bonds}
    
def LuZipper(line): # parsing from DSSR

    num   = int(line[:5])
    line  = cut_double_spaces(line[5:-1]).split(' ')
    if not line[0]: line = line[1:]
    size  = int(line[0][4:])
    nucls = line[2].split(',')
    for i in range(len(nucls)): nucls[i] = '.'.join(nucls[i].split('.')[2:])
    seq   = line[1]

    return {'ID'    :   num,
            'LENGTH':  size,
            'SEQ'   :   seq,
            'NUCLS' : nucls}










