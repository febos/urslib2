import os
import glob
import time

def Into_models(filepath, outfolder):

    right_order = ['_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol',
                   '_atom_site.label_atom_id', '_atom_site.label_alt_id', '_atom_site.label_comp_id',
                   '_atom_site.label_asym_id', '_atom_site.label_entity_id', '_atom_site.label_seq_id',
                   '_atom_site.pdbx_PDB_ins_code', '_atom_site.Cartn_x', '_atom_site.Cartn_y',
                   '_atom_site.Cartn_z', '_atom_site.occupancy', '_atom_site.B_iso_or_equiv',
                   '_atom_site.Cartn_x_esd', '_atom_site.Cartn_y_esd', '_atom_site.Cartn_z_esd',
                   '_atom_site.occupancy_esd', '_atom_site.B_iso_or_equiv_esd', '_atom_site.pdbx_formal_charge',
                   '_atom_site.auth_seq_id', '_atom_site.auth_comp_id', '_atom_site.auth_asym_id',
                   '_atom_site.auth_atom_id', '_atom_site.pdbx_PDB_model_num']

    def Fixed(rows,lin):
        lin = lin.split()
        newline = []
        for i in range(len(right_order)):
            if right_order[i] in rows: newline.append(lin[rows.index(right_order[i])])
        return ' '.join(newline)+' \n'

    title = []
    filename = os.path.basename(filepath)
    models = []
    
    if outfolder[-1] != '/':
        outfolder += '/'    
    
    with open(filepath, 'r') as ciffile:

        hat = True # True - we are stil in headers (before ATOM/HETATM lines)
        current = -1
        can = False
        atom_rows = []
        
        for line in ciffile:

            if can and line[0] == '#': can = False
            if line.startswith('_atom_site.'):
                atom_rows.append(line[:-2])
                can = True

            if can and not line.startswith('_atom_site.'):

                if hat:
                    kk = atom_rows.index('_atom_site.pdbx_PDB_model_num')
                    hat = False

                curmodel = int(line.split()[kk])

                if curmodel != current:
                    current = curmodel
                    if curmodel==0: curmodel=1
                    models.append(open(outfolder + filename + str(curmodel), 'w'))
                    for i in title:
                        models[-1].write(i)
                    if len(atom_rows)!=len(right_order):
                        models[-1].write(' \n'.join([x for x in right_order if x in atom_rows])+' \n')
                    else: models[-1].write(' \n'.join(right_order)+' \n') # _atom_site. Headers

                if atom_rows == right_order: models[-1].write(line)
                else                       : models[-1].write(Fixed(atom_rows,line))

            else:
                if hat:
                    if not can: title.append(line)
                else:
                    for m in models:
                        m.write(line)

    for m in models: m.close()

    print(filename+' is successfully divided into models.',end='\n')


def All(infolder='',outfolder=''):

    good_input = 0

    while not good_input:
        
        if os.path.exists(infolder) and os.path.exists(outfolder): good_input = 1

        elif infolder or outfolder:

            print('Bad paths. Try again...')
            infolder  = input('Please enter the path to folder with cif files: ')
            outfolder = input('Please enter the path to output folder: ')

        else:

            infolder  = input('Please enter the path to folder with cif files: ')
            outfolder = input('Please enter the path to output folder: ')

    Time = time.time()

    if infolder[-1] != '/':
        infolder += '/'  

    files = glob.glob(infolder+'*.cif')

    files.sort()

    total     = len(files)
    counter = 1

    for f in files:

        print(str(counter)+'/'+str(total),end=' ')
        Into_models(f,outfolder)
        counter += 1

    Time = time.time() - Time    
    print('Time: %s min %s sec'%(int(Time//60),int(Time%60)))
    

if __name__ == '__main__':

    All()










        
