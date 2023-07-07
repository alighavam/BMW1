function varargout = bmw1_imana(what,varargin)
% Function for preprocessing of BMW1 data.
% Ali Ghavampour 2023 - Diedrichsen Lab
%
% Project originally started 2015 by Atsushi Yokoi & Diogo Duarte

filePath = '/Users/aghavampour/Desktop/Projects/BimanualWrist/BMW1/';
addpath(strcat(filePath,'tools/'))
addpath(strcat(filePath,'functions/'))
addpath(genpath('dataframe-2016.1/'))

%% Directory specification
baseDir = '/Users/aghavampour/Desktop/Projects/BimanualWrist/BMW1/data';

fieldmapsDir    = fullfile(baseDir, 'fieldmaps');dircheck(fieldmapsDir);
behaviourDir    = fullfile(baseDir, 'behavioural_data');dircheck(behaviourDir);
analyzeDir 		= fullfile(baseDir, 'analyze');dircheck(analyzeDir);
anatomicalDir   = fullfile(baseDir, 'anatomicals');dircheck(anatomicalDir);
imagingDirRaw   = fullfile(baseDir, 'imaging_data_raw');dircheck(imagingDirRaw);
imagingDir      = fullfile(baseDir, 'imaging_data');dircheck(imagingDir);
freesurferDir   = fullfile(baseDir, 'surfaceFreesurfer');dircheck(freesurferDir)
caretDir        = fullfile(baseDir, 'surfaceCaret');dircheck(caretDir);
regDir          = fullfile(baseDir, 'RegionOfInterest');dircheck(regDir);
figDir          = fullfile(baseDir, 'figures');dircheck(figDir);
BGDir           = fullfile(anatomicalDir,'basal_ganglia');dircheck(BGDir);
statsDir        = fullfile(baseDir,'stats'); dircheck(statsDir);
physioDir       = fullfile(baseDir,'physio_data'); dircheck(physioDir);

dirNames        = {'baseDir','behavDir','analyzeDir','imagingDir','imagingDirRaw',...
                   'anatomicalDir','freesurferDir','caretDir','regDir',...
                   'physioDir','figDir','statsDir'};


%% subject info
subj_name = {'suwo04'};
loc_AC = {-[0,0,0]};

nRun = 10;
run             = {...
                    {{'1','2','3','4'},{'5','6','7','8','9','10'}}...
                  };

sess            = {'sess1','sess2'};
subj_sess       = {sess};
subj_runs       = {run};

%% Experiment specific parameters
use3D   = 0;            % in certain envinronment, using 3D image makes analysis faster
startTR = 9; % 8245 ms was start time
startTime = 8245;
numDummys = 3; % discard these volumes
nTR     = 740; %270-1;
numTRs = 740;
TR      = 1000; %2720;
trialTime   = 7000;
planTime    = 2000;
nSlice = 32;

%% MAIN OPERATION =========================================================

switch(what)
    case 'FMAP:makefieldmap'
        prefixepi  = '';
        prefixfieldmap  = '';
        sn      = varargin{1};
        sess    = varargin{2};
        subfolderRawdata    = sprintf('sess%d',sess);
        subfolderFieldmap   = sprintf('sess%d',sess);
        et1 = 4.92;
        et2 = 7.38;
        tert = 90*0.7/2;
        
%         spmj_makefieldmap(baseDir,subj_name{sn},subj_runs{sn}{sess}, ...
%                           'prefix','', ...
%                           'image',nTR-numDummys, ...
%                           'subfolderRawdata',subfolderRawdata, ...
%                           'subfolderFieldmap',subfolderFieldmap);
        makefieldmap(baseDir, subj_name{sn}, subj_runs{sn}{sess},...
                    'prefixepi',prefixepi,...
                    'prefixfieldmap',prefixfieldmap,...
                    'use3D',use3D,'image', nTR-numDummys,...
                    'subfolderRawdata',subfolderRawdata,...
                    'subfolderFieldmap',subfolderFieldmap);
end

















