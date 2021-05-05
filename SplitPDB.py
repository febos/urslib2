import os
import glob
import time

def Into_models(filepath, outfolder):

    title = []
    filename = os.path.basename(filepath)

    if outfolder[-1] != '/':
        outfolder += '/'    
    
    with open(filepath, 'r') as pdbfile:

        hat = True # True - we are stil in headers (before MODEL and ATOM/HETATM lines)
        
        for line in pdbfile:
            
            if line[:5] == 'MODEL':

                hat = False
                number = int(line[10:15])
                model = open(outfolder + filename + str(number), 'w')
                for i in title:
                    model.write(i)
                model.write(line)
                
            elif hat: # if there is no MODEL line in pdb-file: title == whole pdb-file
            
                title.append(line)

                if line == 'END                                                                             \n' or\
                   line == 'END                                                                             ':

                    model = open(outfolder + filename + str(1), 'w')
                    for i in title:
                        model.write(i)
                    model.close()
            
            elif line[:6] == 'ENDMDL':

                model.write(line)
                model.write('END                                                                             \n')
                model.close()
                
            elif not model.closed:
                
                model.write(line)    

    print(filename+' is successfully divided into models.',end='\n')


def All(infolder='',outfolder=''):

    good_input = 0

    while not good_input:
        
        if os.path.exists(infolder) and os.path.exists(outfolder): good_input = 1

        elif infolder or outfolder:

            print('Bad paths. Try again...')
            infolder  = input('Please enter the path to folder with pdb files: ')
            outfolder = input('Please enter the path to output folder: ')

        else:

            infolder  = input('Please enter the path to folder with pdb files: ')
            outfolder = input('Please enter the path to output folder: ')

    Time = time.time()

    if infolder[-1] != '/':
        infolder += '/'  

    files = glob.glob(infolder+'*.pdb')

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










        
