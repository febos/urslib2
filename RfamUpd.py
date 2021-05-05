
import urllib.request
import gzip,os



def update(folder):

    print('Beginning Rfam download with urllib...')

    ftp = 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/'

    files = [['family.txt.gz','rfam_family.txt.gz'],
             ['pdb_full_region.txt.gz','rfam_pdb.txt.gz']]

    for file in files: urllib.request.urlretrieve(ftp+file[0], folder+file[1])

    print('Unpacking...')

    for file in [x[1] for x in files]:
        
        outp = open(folder+file.replace('.gz',''),'wb')
        
        with gzip.open(folder+file,'rb') as f:
            content = f.read()
            outp.write(content)

        outp.close()
        os.remove(folder+file)
            
    print('Complete.')




if __name__=='__main__':

    folder = ''
    
    update(folder)

