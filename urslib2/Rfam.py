
import urllib.request
import gzip, os



def UpdateRfam(folder,
               rfam_ftp = 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/'):

    print('Beginning Rfam download with urllib...')

    files = ['family.txt.gz', 'pdb_full_region.txt.gz']

    for file in files: urllib.request.urlretrieve(os.path.join(rfam_ftp,file),
                                                  os.path.join(folder,file))

    print('Unpacking...')

    for file in files:
        
        outp = open(os.path.join(folder,file.replace('.gz','')),'wb')
        
        with gzip.open(os.path.join(folder,file),'rb') as f:
            content = f.read()
            outp.write(content)

        outp.close()
        os.remove(os.path.join(folder,file))
            
    print('Complete.')



def GetRfamInfo(path_to_rfam_files = '',
                update = False): 
    '''
    Parameters:
        path_to_rfam_files - path to family.txt and pdb_full_region.txt (by default - working directory)
        update - True/False - download new family.txt and pdb_full_region.txt
    Returns:
        dictionary {PDB-ID: {CHAIN: RFAM_INFO}}
    '''
    path_to_family = os.path.join(path_to_rfam_files,'family.txt')
    path_to_pdb    = os.path.join(path_to_rfam_files,'pdb_full_region.txt')

    if not os.path.exists(path_to_family) or \
       not os.path.exists(path_to_pdb) or \
       update == True:

        UpdateRfam(path_to_rfam_files)

    rfam = {}


    return rfam    

    
    

    

    


if __name__=='__main__':

    folder = ''
    
    GetRfamInfo(folder)

