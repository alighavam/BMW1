# creates target files for the BMW1 experiment

import numpy as np
import pandas as pd
import random
import os

def make_train_tgt(targetPath, sampleSeed):

    # Find the last train subject name:
    lastSubj = 'suow01'
    for file in os.listdir(targetPath):
        if file.endswith('.tgt'):
            if file[7:12] == 'train':
                if lastSubj < file[0:6]:
                    lastSubj = file[0:6]
    
    # create new train subject based on the last subject:
    newSubj = 'suwo{:02d}'.format(int(lastSubj[-2:])+1)            

    # read sample .tgt file, shuffle and make new .tgt file:
    for file in os.listdir(targetPath):
        if file.endswith('.tgt'):
            if file[0:12] == sampleSeed+'_train':
                sampleTarget = pd.read_csv(targetPath + '/' + file, sep='\t')
                sampleTargetShuffled = sampleTarget.sample(frac=1)
                newTargetName = newSubj + '_' + file[7:]
                sampleTargetShuffled.to_csv(targetPath + '/' + newTargetName, sep='\t', index=False)

def make_startTime(tStart,nTotal,nRest,iti,itiRest):

    # Set the length of the array
    length = nTotal

    # Initialize the array with zeros
    trial_time_array = np.zeros(length, dtype=int)

    # Set the first element
    trial_time_array[0] = tStart

    # Set the increment for regular elements
    regular_increment = iti

    # Set the increment for rest elements
    rest_increment = itiRest

    # Set the indices for random elements
    random_indices = np.random.choice(np.arange(1,length), nRest, replace=False)

    # Generate the array
    for i in range(1, length):
        if i in random_indices:
            trial_time_array[i] = trial_time_array[i - 1] + rest_increment
        else:
            trial_time_array[i] = trial_time_array[i - 1] + regular_increment

    return trial_time_array

def make_scan_tgt(targetPath, sampleSeed):

    # Find the last scan subject name:
    lastSubj = 'suow01' # you don't need to change this. This is just an initial value.
    for file in os.listdir(targetPath):
        if file.endswith('.tgt'):
            if file[7:11] == 'scan':
                if lastSubj < file[0:6]:
                    lastSubj = file[0:6]
    
    # create new train subject based on the last subject:
    newSubj = 'suwo{:02d}'.format(int(lastSubj[-2:])+1)            

    # read sample .tgt file, shuffle and make new .tgt file:
    for file in os.listdir(targetPath):
        if file.endswith('.tgt'):
            if file[0:11] == sampleSeed+'_scan':
                sampleTarget = pd.read_csv(targetPath + '/' + file, sep='\t')
                sampleTargetShuffled = sampleTarget.sample(frac=1)

                # make startTime vector:
                startTime = make_startTime(8245,96,9,7055,12070)

                # add the new startTime vector to the dataframe:
                sampleTargetShuffled.startTime = startTime

                # save the target file:
                newTargetName = newSubj + '_' + file[7:]
                sampleTargetShuffled.to_csv(targetPath + '/' + newTargetName, sep='\t', index=False)


dir_path = os.getcwd()
sampleSeed = 'suwo01'
make_train_tgt(dir_path + '/target/bimanualMRuwo_target/', sampleSeed)
make_scan_tgt(dir_path + '/target/bimanualMRuwo_target/', sampleSeed)

