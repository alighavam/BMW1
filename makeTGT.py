# creates target files for the BMW1 experiment

import numpy as np
import pandas as pd
import os

def make_train_tgt(targetPath, sampleSubjName):

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
            if file[0:12] == sampleSubjName+'_train':
                sampleTarget = pd.read_csv(targetPath + '/' + file, sep='\t')
                sampleTargetShuffled = sampleTarget.sample(frac=1)
                newTargetName = newSubj + '_' + file[7:]
                sampleTargetShuffled.to_csv(targetPath + '/' + newTargetName, sep='\t', index=False)


dir_path = os.getcwd()
sampleSubjName = 'suwo01'
sampleTarget = make_train_tgt(dir_path + '/target/bimanualMRuwo_target/', sampleSubjName)
# print(sampleTarget.u_or_b)
