#!/usr/bin/env python3


import os
import shutil
import glob

the_folder=''
folder_to_save='movie_no_Log_err_Quintic_8x8x1_dedner=1.0_eta=1.36'
snapshot_key='OrszagTangVortex_'
parameter_file='OrszagTangVortex.yml'
statistics_file='statistics.txt'

#getting snapshot addresses
addressbook=glob.glob(the_folder+'*.hdf5')
snapshots=sorted([addressbook[i] for i in range(len(addressbook)) if snapshot_key in os.path.splitext(addressbook[i])[0]])


#creating folder for movie storage
try:
    os.mkdir(folder_to_save)
except:
    print(f'the folder already exists')
    
print('copying config and statistics')
shutil.copyfile(parameter_file,folder_to_save+'/'+parameter_file)
shutil.copyfile(statistics_file,folder_to_save+'/'+statistics_file)

#executing plotSolution.py for all snapshots
for i in range(len(snapshots)):
    print('')
    print(f'Processing snapshot {i}')
    filename=os.path.splitext(snapshots[i])[0]
    ext='.png'
    file_to_open = folder_to_save+'/' + filename+ext
    cwd = os.path.join(os.getcwd(),'plotSolution_Nikyta_2.py '+snapshots[i]+' '+file_to_open)
    os.system('{} {}'.format('python', cwd))

#running statistics
cwd= os.path.join(os.getcwd(),'statistics_reader.py '+folder_to_save)
os.system('{} {}'.format('python', cwd))
