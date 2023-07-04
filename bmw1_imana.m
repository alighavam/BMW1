function varargout=bmw1_imana(what,varargin)
%%
%
%
%
%
%
% Atsushi Yokoi & Diogo Duarte (2015)
% Ali Ghavampour (2023)

% use spm8
%selectSPMversion('spm8');

% change default
DefaultFigureUnit = get(0,'defaultFigureUnit');
%set(0,'DefaultFigureUnit','centimeter');

%% Change spm version to 12
% This code relies manly on spm 12. This is because of "fast" aquisition
% option within the 1st-level glm design/estimation.
% When defining Basal Ganglia ROIs, you may need to use spm 8
% when running segmentation.
if ~selectSPMversion('spm12');
    warning('Failed to set path for spm12!\nPlease check if spm12 is installed at right directory of your computer.');
    %     return;
end
spm_get_defaults('mat.format','-v7.3');

%% Get runtime info (not finished yet)
% This is to take the record of each function execution, recording input
% arguments. Ideally shold be attached to each resultant data/figures.
% Any smarter ida about how to implement this?
persistent this
if isempty(this)
    this.function   = mfilename;
    this.fullpath   = mfilename('fullpath');
    try
    this.spmversion = spm('version');
    catch
        this.spmversion = spm('ver');
    end
    %[this.flist,this.plist] = matlab.codetools.requiredFilesAndProducts(this.function,'toponly');
    % this takes long...
end
this.case       = what;
this.varargin   = varargin;

%% Directory specification
baseDir = '/Volumes/G_Thunderbolt/Yokoi_Research/data/BimanualWrist_MR/bmw1_uwo/'; % Atsushi

fieldmapsDir    = fullfile(baseDir, 'fieldmaps');dircheck(fieldmapsDir);
behaviourDir    = fullfile(baseDir, 'behavioural_data');dircheck(behaviourDir);
analyzeDir 		= fullfile(baseDir, 'analyze');dircheck(analyzeDir);
dicomDir        = [baseDir '/imaging_data_dicom']; dircheck(dicomDir);
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
    'dicomDir','anatomicalDir','freesurferDir','caretDir','regDir',...
    'physioDir','figDir','statsDir'};

%% Setup freesurfer
% not sure if this is necessary
setenv('SUBJECTS_DIR', freesurferDir);

%% Subject info (must be updated at every inclusion of new subject)
%
%
subj_name = {'suwo01','suwo02','suwo03'};
DicomFolderName = {{'2016_11_08_suwo01_v2',...
    '2016_12_06_SUWO01_v3'};
    {'2016_11_07_suwo02',...
    '2016_11_08_suwo02_v2'};
    {'2016_11_07_SUWO03',...
    '2016_11_09_suwo03_v2'}
    };
DicomName = {{'2016_11_08_SUWO01_V2.MR.DIEDRICHSEN_BIMANUAL',...
    '2016_12_06_SUWO01_V3.MR.DIEDRICHSEN_BIMANUAL'};
    {'2016_11_07_SUWO02.MR.DIEDRICHSEN_BIMANUAL',...
    '2016_11_08_SUWO02_V2.MR.DIEDRICHSEN_BIMANUAL'};
    {'2016_11_07_SUWO03.MR.DIEDRICHSEN_BIMANUAL',...
    '2016_11_09_SUWO03_V2.MR.DIEDRICHSEN_BIMANUAL'}
    };
NiiRawName = {{'161108140208STD131221107524367007',...
    '161206151253STD131221107524367007'},...
    {'161107115943STD131221107524367007',...
    '161108095900STD131221107524367007'},...
    {'161107134410STD131221107524367007',...
    '161109142145STD131221107524367007'},...
    };
fscanNum    = {{[2:6],[6:10]},...
    {[2:6],[2:4,6,7]},...
    {[2:4,8,9],[2:6]}};
fmapNum     = {{[8,9],[11,12]},...
    {[7,8],[8,9]},...
    {[10,11],[7,8]}};
T1Num     = {{[7],[]},...
    {[9],[]},...
    {[12],[]}};
T2Num     = {{[10],[]},...
    {[10],[]},...
    {[5],[]}};
loc_AC      = {-[93,136,126], -[91,134,134], -[88,145,136]};
% The values of loc_AC should be acquired manually prior to the preprocessing
%   Step 1: get .nii file of anatomical data by running "spmj_tar2nii(TarFileName,NiiFileName)"
%   Step 2: open .nii file with MRIcron and manually find AC and read the xyz coordinate values
%           (note: there values are not [0 0 0] in the MNI coordinate)
%   Step 3: set those values into loc_AC (subtract from zero)

subj_name_behav = {'suwo01','suwo02','suwo03'};

nRun = 10;
run             = {{'1','2','3','4','5'},{'6','7','8','9','10'}};%
sess            = {'sess1','sess2'};%
subj_sess       = {sess,sess,sess}; %
subj_runs       = {run,run,run};

%% Experiment specific parameters
use3D   = 0;            % in certain envinronment, using 3D image makes analysis faster
startTR = 9; % ? 8245 ms was start time
startTime = 8245;
numDummys = 3; % discard these volumes
nTR     = 740; %270-1;
numTRs = 740;
TR      = 1000; %2720;
trialTime   = 7000;
planTime    = 2000;
nSlice = 32;

%% Freesurfer & ROI parameters
%mysetEnv(); % set environmental variables (neccessary for use freesurfer, fls, and caret)
atlasA      = {'x'};
atlasname   = {'fsaverage_sym'};
hem         = {'lh','rh'};
hemName     = {'LeftHem','RightHem'};
hemName_s   = {'Left','Right'};
numregions_surf = 8;
numregions_BG   = 4;
numregions      = numregions_surf+numregions_BG;
regSide     = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 2];
regType     = [1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 9 10 11 12 9 10 11 12];
regname         = {'S1','M1','PMd','PMv','SMA','V12','SPLa','SPLp','CaudateN' 'Pallidum', 'Putamen' 'Thalamus'};
regname_surf    = {'S1','M1','PMd','PMv','SMA','V12','SPLa','SPLp'};
regname_BG      = {'CaudateN' 'Pallidum', 'Putamen' 'Thalamus'};

%% Analysis parameters
Contrasts      = {''};
nConditions = 6+6+36; % 6 directions per each unimanual and combination of them for bimanual 8+8+64 was for pilot data;
nDirections = 6;
Directions = [0:60:300];

%% GLM def
glmName     = {'GLM_firstlevel_6direction_fasthpf','GLM_firstlevel_6direction_fastnohpf',...
    'GLM_firstlevel_6direction_rwlshpf','GLM_firstlevel_6direction_rwlsnohpf'};
glmtypes    = {'glm_design_6direction_fasthpf','glm_design_6direction_fastnohpf',...
    'glm_design_6direction_rwlshpf','glm_design_6direction_rwlsnohpf'};

varargout = {};
if ismac
    mysetEnv(); % Atsushi
elseif isunix
    setenvLinux();
end

%% Plot parameters
% color
purple = [255 102 178]/255;
gold = [255 178 102]/255;
darkpurple = [204 0 102]/255;
darkgold = [204 102 0]/255;
red = [255 102 102]/255;
blue = [102 102 255]/255;
darkred = [204 0 0]/255;
darkblue = [0 0 204]/255;

%% Main operation
switch(what)
    case 'dir'                  % Show directory info
        
        for d=1:numel(dirNames);
            D.(dirNames{d}) = eval(dirNames{d});
        end
        
        varargout = {D};
    case 'demographic'
        sn = 1:9;
        cd(behavDir); D = [];
        for s = sn
            isfile = [];
            % search directory which contains subjectinfo.txt
            infoDir = dir(sprintf('*%s*',subj_name{s}));
            for f=1:numel(infoDir);
                isfile(f) = exist(fullfile(infoDir(f).name,'subjectinfo.txt'),'file');
            end
            
            if any(isfile)
                S = dload(fullfile(behavDir,infoDir(find(isfile>0,1,'first')).name,'subjectinfo.txt'));
            end
            S.subj= s;
            S.exp = {'bmw1_uwo'};
            D = addstruct(D,S);
        end
        
        varargout = {D};
    case 'check_ima_file'       % Check file numbers to which run to be taken (excluding incomplete runs due to crash etc.)
        s = varargin{1};
        nimage = {0};
        % Get .ima file info
        F = Dir(fullfile(dicomDir,subj_name{s},[DicomName{s},'*.IMA']));
        for f=1:length(F);
            fparts = strsplit(F(f).name,'.');
            nscan = str2double(fparts{4});
            if nscan>length(nimage)
                nimage{end+1} = 0;
            end
            nimage{nscan} = nimage{nscan}+1;
        end
        disp('.IMA file information')
        fprintf('scan no.\timage size\n');
        disp('----------------------------')
        for i=1:numel(nimage)
            fprintf('%4.4d\t%d\n',i,nimage{i});
        end
        varargout = {nimage};
    case 'check_physio_file'    % Check file numbers to which run to be taken (excluding incomplete runs due to crash etc.)
        s = varargin{1};
        
        % Get .resp and .puls file info
        R = dir(fullfile(dicomDir,subj_name{s},'Physio_log_*.resp'));
        P = dir(fullfile(dicomDir,subj_name{s},'Physio_log_*.puls'));
        if isempty(R)||isempty(P)
            warning('Physio log file missing.');
            return;
        end
        for r=1:length(R);
            fparts = strsplit(R(r).name,'_');
            lognum_resp(r) = str2double(fparts{4});
            % append full path
            R(r).name = fullfile(dicomDir,subj_name{s},R(r).name);
        end
        for r=1:length(P);
            fparts = strsplit(P(r).name,'_');
            lognum_puls(r) = str2double(fparts{4});
            % append full path
            P(r).name = fullfile(dicomDir,subj_name{s},P(r).name);
        end
        disp('Physio log file information')
        fprintf('resp log no.\tpuls log no.\n');
        disp('----------------------------')
        for i=1:length(R)
            fprintf('%4.4d\t%4.4d\n',lognum_resp(i),lognum_puls(i));
        end
        varargout = {R,P};
        
    case 'preprocess'           % All preprocessing by just one go (AC coordinates (loc_AC) are prerequisite)
        sn = varargin{1};
        for s=sn
            bmw1_uwo_imana('dicom_import',s);
            bmw1_uwo_imana('make_4dNifti',s);
            bmw1_uwo_imana('realign_unwarp',s);  % functional
            bmw1_uwo_imana('move_images',s);      % functional
            bmw1_uwo_imana('copy_images_to_sess1',s);      % functional
            bmw1_uwo_imana('reslice_LPI',s);      % anatomical
            bmw1_uwo_imana('centre_AC',s);        % anatomical
            bmw1_uwo_imana('segmentation',s);     % antomical
            bmw1_uwo_imana('meanimage_bias_correction',s,1);           % functional
            bmw1_uwo_imana('coreg',s);           % functional
            bmw1_uwo_imana('make_samealign',s);  % functional
            bmw1_uwo_imana('check_samealign',s); % functional
            bmw1_uwo_imana('link_processed_images',s); % functional
            bmw1_uwo_imana('make_maskImage',s);  % anatomical/functional
        end
    case 'dicom_import'         % Import dicom file into Nifti file
        % Possibly do this for a run per loop to preserve memory
        s = varargin{1};
        ss = varargin{2};
        type = varargin{3};
        
        cd(fullfile(dicomDir,subj_name{s},DicomFolderName{s}{ss}));
        
        switch type
            case 'functional'
                scanNum = fscanNum;
                postfix = 'func';
                numTRs = nTR;
            case 'anatomical'
                scanNum = T1Num;
                postfix = 'anat';
                numTRs = 176;
            case 'fieldmap'
                scanNum = fmapNum;
                postfix = 'fmap';
                numTRs = 46;
        end
        
        Nseries = numel(scanNum{s}{ss});
        if Nseries==0
            return;
        end
        for r = 1:Nseries
            % Get DICOM FILE NAMES
            scanNumber  = scanNum{s}{ss}(r);
            DIR         = dir(sprintf('%s.%4.4d.*.IMA',DicomName{s}{ss},scanNumber));
            Names       = vertcat(DIR.name);
            
            % Check if volume number is correct and convert to nii
            if (numTRs~=length(DIR))&&(2*numTRs~=length(DIR))
                warning('Number of volume (%d) doesn''t match!\n',length(DIR));
                keyboard();
            end;
            
            % Get dicom headers
            HDR = spm_dicom_headers(Names,1);
            dirname{r} = sprintf('series%2.2d.%s',scanNumber,postfix); % create new directory and cd
            dircheck(dirname{r});
            cd(fullfile(dicomDir,subj_name{s},DicomFolderName{s}{ss},dirname{r}));
            %SA 01/20/2016: need to change dir b/c cannot push to dir in spm_dicom_convert for spm12b (check with Joern about this..)
            
            % Do spm_dicom_convert
            spm_dicom_convert(HDR,'all','flat','nii');%,dirname{r});  % Convert the data to nifti
            cd .. %SA 01/20/2016
            
            % Create dicom folders to move all .IMA files for each run? this may make physio stuff easier
            
        end;
    case 'make_4dNifti'         % Make 4D niftis out of your raw data files
        s = varargin{1};
        ss = varargin{2};
        
        % Make new directory
        subjDir = fullfile(imagingDirRaw,subj_name{s});dircheck(subjDir);
        sessDir = fullfile(subjDir, subj_sess{s}{ss}); dircheck(sessDir);
        
        % Merge .nii files into 4D Nifti data
        for i=1:length(fscanNum{s}{ss})
            scanNum = fscanNum{s}{ss}(i);
            outfilename = fullfile(sessDir,sprintf('%s_run_%s.nii',subj_name{s},subj_runs{s}{ss}{i}));
            k=0; P = {};
            for j = (numDummys+1):nTR
                k=k+1;
                P{k}=fullfile(dicomDir,subj_name{s},DicomFolderName{s}{ss},sprintf('series%2.2d.func',scanNum),...
                    sprintf('f%s-%4.4d-%5.5d-%6.6d-01.nii',NiiRawName{s}{ss},scanNum, j, j));
            end;
            spm_file_merge(char(P),outfilename);
            fprintf('%d\n',i);
        end;
    case 'move_fmap_anatomical'        % Move field map and anatomical images
        s = varargin{1};
        ss = varargin{2};
        
        % Dicom directory
        subjDir_dicom = fullfile(dicomDir,subj_name{s},DicomFolderName{s}{ss});
        
        % Move .nii files
        subjDir_fmap = fullfile(fieldmapsDir,subj_name{s});dircheck(subjDir_fmap);
        sessDir_fmap = fullfile(subjDir_fmap, subj_sess{s}{ss}); dircheck(sessDir_fmap);
        for i=1:length(fmapNum{s}{ss})
            scanNum = fmapNum{s}{ss}(i);
            switch i
                case 1
                    source = fullfile(subjDir_dicom, sprintf('series%2.2d.fmap',scanNum),...
                        sprintf('%s_magnitude.nii',subj_name{s}));
                    
                    dest = fullfile(sessDir_fmap,sprintf('%s_magnitude.nii',subj_name{s}));
                    
                    % make 4D-nifti
                    P{1} = fullfile(subjDir_dicom, sprintf('series%2.2d.fmap',scanNum),...
                        sprintf('s%s-%4.4d-%5.5d-%6.6d-01.nii',NiiRawName{s}{ss},scanNum, 1, 1));
                    P{2} = fullfile(subjDir_dicom, sprintf('series%2.2d.fmap',scanNum),...
                        sprintf('s%s-%4.4d-%5.5d-%6.6d-02.nii',NiiRawName{s}{ss},scanNum, 1, 1));
                    spm_file_merge(char(P),source);
                case 2
                    source = fullfile(subjDir_dicom, sprintf('series%2.2d.fmap',scanNum),...
                        sprintf('s%s-%4.4d-%5.5d-%6.6d-02.nii',NiiRawName{s}{ss},scanNum, 1, 1));
                    
                    dest = fullfile(sessDir_fmap,sprintf('%s_phase.nii',subj_name{s}));
            end
            movefile(source,dest);
            
            fprintf('%d\n',i);
        end;
        % Move .nii files
        subjDir_ana = fullfile(anatomicalDir,subj_name{s});dircheck(subjDir_ana);
        if length(T1Num{s}{ss})==0
            return;
        end
        for i=1:length(T1Num{s}{ss})
            scanNum = T1Num{s}{ss}(i);
            
            source = fullfile(subjDir_dicom, sprintf('series%2.2d.anat',scanNum),...
                sprintf('s%s-%4.4d-%5.5d-%6.6d-01.nii',NiiRawName{s}{ss},scanNum, 1, 176));
            
            dest = fullfile(subjDir_ana,sprintf('%s_anatomical_raw.nii',subj_name{s}));
            
            movefile(source,dest);
            
            fprintf('%d\n',i);
        end;
    case 'make_fmap'                % Create field map image (VDM)
        prefixepi  = ''; %'a' for slice timing-corrected data;
        prefixfieldmap  = '';
        sn      = varargin{1};
        sess    = varargin{2};
        subfolderRawdata    = sprintf('sess%d',sess);
        subfolderFieldmap   = sprintf('sess%d',sess);
        
        %spmj_makefieldmap(baseDir, subj_name{sn}, subj_runs{sn}{sess},...
        makefieldmap(baseDir, subj_name{sn}, subj_runs{sn}{sess},...
            'prefixepi',prefixepi,...
            'prefixfieldmap',prefixfieldmap,...
            'use3D',use3D,'image', nTR-numDummys,...
            'subfolderRawdata',subfolderRawdata,...
            'subfolderFieldmap',subfolderFieldmap);
    case 'realign_unwarp'      % Do spm_realign and spm_unwarp
        sn = varargin{1};
        %sess = varargin{2};
        prefix  =''; % 'a' for slice-time-corrected data
        prefixfieldmap = '';
        for s=sn
            if (1)
                subfolderRawdata    = subj_sess{s}; %{'sess1'};
                subfolderFieldmap   = subj_sess{s}; %{'sess1'};
                subj_runs_               = {{'_1','_2','_3','_4','_5'},{'_6','_7','_8','_9','_10'}};%;%subj_runs{s}; %run;
                
                spmj_realign_unwarp_sess(baseDir, subj_name{s}, subj_runs_, nTR,...
                    'prefix',prefix, 'use3D', use3D, 'subfolderRawdata',subfolderRawdata,...
                    'subfolderFieldmap',subfolderFieldmap);
            else
                for ss=sess
                    subfolderFieldmap = sprintf('sess%d',ss);
                    subfolderRawdata  = sprintf('sess%d',ss);
                    
                    % specify target image (last functional volume of the session)
                    if use3D
                        targetImg = fullfile(imagingDirRaw, subj_name{s} ,subfolderRawdata, [prefix subj_name{s} ,'_run_',run{ss}{end},'_',num2str(nTR-numDummys),'.nii']);
                    else
                        targetImg= fullfile(imagingDirRaw, subj_name{s} ,subfolderRawdata, [prefix subj_name{s}, '_run_',run{ss}{end}, '.nii,' ,num2str(nTR-numDummys)]);
                    end;
                    
                    realign_unwarp(baseDir, subj_name{s}, run{ss}, 1, nTR-numDummys,...
                        'prefix',prefix, 'use3D', use3D,...
                        'prefixfieldmap',prefixfieldmap,...
                        'subfolderFieldmap',subfolderFieldmap,...
                        'subfolderRawdata',subfolderRawdata,...
                        'targetImg', targetImg);
                end
                
            end
        end
    case 'realign'                  % Do spm_realign only
        sn = varargin{1};
        sess = varargin{2};
        prefix  =''; % 'a' for slice-time-corrected data
        for s=sn
            for ss=sess
                subfolderRawdata  = sprintf('sess%d',ss);
                
                % specify target image (last functional volume of the session)
                if use3D
                    targetImg = fullfile(imagingDirRaw, subj_name{s} ,subfolderRawdata, [prefix subj_name{s} ,'_run_',run{ss}{end},'_',num2str(nTR-numDummys),'.nii']);
                else
                    targetImg= fullfile(imagingDirRaw, subj_name{s} ,subfolderRawdata, [prefix subj_name{s}, '_run_',run{ss}{end}, '.nii,' ,num2str(nTR-numDummys)]);
                end;
                
                realign(baseDir, subj_name{s}, run{ss}, 1, nTR-numDummys,...
                    'prefix',prefix, 'use3D', use3D,...
                    'subfolderRawdata',subfolderRawdata,...
                    'targetImg', targetImg);
            end
        end
    case 'move_images'          % Move images created by realine_unwarp into imaging_data
        sn = varargin{1};
        prefix='u';
        vararginoptions(varargin(2:end), {'prefix'});
        
        for s=sn
            for ss=1:length(subj_runs{s})
                dircheck(fullfile(imagingDir,subj_name{s}, subj_sess{s}{ss}));
                for r=1:length(subj_runs{s}{ss});
                    source = fullfile(imagingDirRaw,subj_name{s}, subj_sess{s}{ss}, sprintf('%s%s_run_%s.nii',prefix,subj_name{s},subj_runs{s}{ss}{r}));
                    dest = fullfile(imagingDir,subj_name{s}, subj_sess{s}{ss}, sprintf('%s%s_run_%s.nii',prefix,subj_name{s},subj_runs{s}{ss}{r}));
                    try
                        movefile(source,dest);
                    catch message
                        warning(message.message);
                    end
                    source = fullfile(imagingDirRaw,subj_name{s}, subj_sess{s}{ss}, sprintf('rp_%s_run_%s.txt',subj_name{s},subj_runs{s}{ss}{r}));
                    dest = fullfile(imagingDir,subj_name{s}, subj_sess{s}{ss}, sprintf('rp_%s_run_%s.txt',subj_name{s},subj_runs{s}{ss}{r}));
                    try
                        movefile(source,dest);
                    catch message
                        warning(message.message);
                    end
                end;
                source = fullfile(imagingDirRaw,subj_name{s}, subj_sess{s}{ss}, sprintf('mean%s%s_run_%s.nii',prefix,subj_name{s}, subj_runs{s}{ss}{1}));
                dest = fullfile(imagingDir,subj_name{s}, subj_sess{s}{ss}, sprintf('%smeanepi_%s.nii',prefix,subj_name{s}));
                try
                    movefile(source,dest);
                catch message
                    warning(message.message);
                end
            end
        end
    case 'reslice_LPI'          % Reslice anatomical image within LPI coordinate systems
        sn  = varargin{1};
        
        % (1) Reslice anatomical image to set it within LPI co-ordinate frames
        source  = fullfile(anatomicalDir,subj_name{sn},[subj_name{sn}, '_anatomical_raw','.nii']);
        dest    = fullfile(anatomicalDir,subj_name{sn},[subj_name{sn}, '_anatomical','.nii']);
        spmj_reslice_LPI(source,'name', dest);
        
        % (2) In the resliced image, set translation to zero
        V               = spm_vol(dest);
        dat             = spm_read_vols(V);
        V.mat(1:3,4)    = [0 0 0];
        spm_write_vol(V,dat);
    case 'centre_AC'            % Re-center AC
        sn      = varargin{1};
        
        img    = fullfile(anatomicalDir,subj_name{sn},[subj_name{sn}, '_anatomical','.nii']);
        V               = spm_vol(img);
        dat             = spm_read_vols(V);
        V.mat(1:3,4)    = loc_AC{sn};
        spm_write_vol(V,dat);
    case 'segmentation'         % Segmentation + Normalization
        sn      = varargin{1};
        usespm8 = 1;
        vararginoptions(varargin(2:end),{'usespm8'});
        
        if usespm8
            selectSPMversion('spm8');
            spm fmri
            spm_jobman('initcfg');
        end
        ver = spm('ver');
        switch ver
            case 'SPM8'
                spmj_segmentation(fullfile(anatomicalDir,subj_name{sn},[subj_name{sn}, '_anatomical','.nii']));
            case 'SPM12'
                SPMhome=fileparts(which('spm.m'));
                J=[];
                for s=sn
                    J.channel.vols = {fullfile(anatomicalDir,subj_name{s},[subj_name{s},'_anatomical.nii,1'])};
                    J.channel.biasreg = 0.001;
                    J.channel.biasfwhm = 60;
                    J.channel.write = [0 0];
                    J.tissue(1).tpm = {fullfile(SPMhome,'tpm/TPM.nii,1')};
                    J.tissue(1).ngaus = 1;
                    J.tissue(1).native = [1 0];
                    J.tissue(1).warped = [0 0];
                    J.tissue(2).tpm = {fullfile(SPMhome,'tpm/TPM.nii,2')};
                    J.tissue(2).ngaus = 1;
                    J.tissue(2).native = [1 0];
                    J.tissue(2).warped = [0 0];
                    J.tissue(3).tpm = {fullfile(SPMhome,'tpm/TPM.nii,3')};
                    J.tissue(3).ngaus = 2;
                    J.tissue(3).native = [1 0];
                    J.tissue(3).warped = [0 0];
                    J.tissue(4).tpm = {fullfile(SPMhome,'tpm/TPM.nii,4')};
                    J.tissue(4).ngaus = 3;
                    J.tissue(4).native = [1 0];
                    J.tissue(4).warped = [0 0];
                    J.tissue(5).tpm = {fullfile(SPMhome,'tpm/TPM.nii,5')};
                    J.tissue(5).ngaus = 4;
                    J.tissue(5).native = [1 0];
                    J.tissue(5).warped = [0 0];
                    J.tissue(6).tpm = {fullfile(SPMhome,'tpm/TPM.nii,6')};
                    J.tissue(6).ngaus = 2;
                    J.tissue(6).native = [0 0];
                    J.tissue(6).warped = [0 0];
                    J.warp.mrf = 1;
                    J.warp.cleanup = 1;
                    J.warp.reg = [0 0.001 0.5 0.05 0.2];
                    J.warp.affreg = 'mni';
                    J.warp.fwhm = 0;
                    J.warp.samp = 3;
                    J.warp.write = [0 0];
                    matlabbatch{1}.spm.spatial.preproc=J;
                    spm_jobman('run',matlabbatch);
                end;
        end
    case 'meanimage_bias_correction'                                         % correct bias for mean image of 1st session
        prefix  = 'u';
        if use3D
            postfix = '_1';
        else
            postfix = '';
        end;
        sn      = varargin{1};
        sess    = varargin{2};
        vararginoptions(varargin(3:end), {'prefix', 'postfix'});
        for s=sn
           for  ss=sess
               P{1}    = fullfile(baseDir, 'imaging_data',subj_name{s},subj_sess{s}{ss},[prefix, 'meanepi_',  subj_name{s}, '.nii']);
               spmj_bias_correct(P);
           end
        end
    case 'coreg'                                                      % coregister rbumean image to anatomical image for each session
        prefix = 'u'; % if we use non-unwarped image, use 'r' instead.
        sn = varargin{1};
        sess = 1;
        vararginoptions(varargin(3:end), {'prefix','sess'});
        
        % (1) Manually seed the functional/anatomical registration
        % - Do "coregtool" on the matlab command window
        % - Select anatomical image and bmeanepi image to overlay
        % - Manually adjust bmeanepi image and save result as rbmeanepi
        %   image
        set(0,'DefaultFigureUnit','pixels');
        cd(fullfile(anatomicalDir,subj_name{sn}));
        coregtool;
        keyboard();
        
        % (2) Run automated co-registration to register bias-corrected meanimage to anatomical image
        J.ref = {fullfile(anatomicalDir,subj_name{sn},[subj_name{sn}, '_anatomical','.nii'])};
        J.source = {fullfile(imagingDir,subj_name{sn},subj_sess{sn}{sess},sprintf('rb%smeanepi_%s.nii', prefix, subj_name{sn}))}; 
        %J.source = {fullfile(imagingDir,subj_name{sn},subj_sess{sn}{sess},sprintf('b%smeanepi_%s.nii', prefix, subj_name{sn}))}; 
        J.other = {''};
        J.eoptions.cost_fun = 'nmi';
        J.eoptions.sep = [4 2];
        J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        J.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estimate=J;
        spm_jobman('run',matlabbatch);
        
        % (3) Check alignment manually
        coregtool;
        keyboard();
    case 'make_samealign'                                       % align to first image (rbmeanepi_* of first session)
        prefix  = 'u';
        sn      = varargin{1};
        target   = 1;
        vararginoptions(varargin(2:end), {'prefix','target'});
        
        cd(fullfile(imagingDir,subj_name{sn}));
        
        nSess = numel(subj_runs{sn});
        Q={};
        for k=1:nSess            
            for r=1:numel(subj_runs{sn}{k})
                runname = subj_runs{sn}{k}{r};
                for i=1:nTR
                    if use3D
                        Q{end+1}    = [fullfile(imagingDir,subj_name{sn},sprintf('%s',subj_sess{sn}{k}),[prefix, subj_name{sn},'_run_',runname,'_',num2str(i),'.nii'])];
                    else
                        Q{end+1}    = [fullfile(imagingDir,subj_name{sn},sprintf('%s',subj_sess{sn}{k}),[prefix, subj_name{sn},'_run_',runname,'.nii,',num2str(i)])];
                    end;
                end;
            end;            
        end;
        P{1} = [fullfile(imagingDir,subj_name{sn},sprintf('%s',subj_sess{sn}{target}),['rb' prefix 'meanepi_' subj_name{sn} '.nii'])];
        spmj_makesamealign_nifti(char(P),char(Q));
    case 'check_samealign'                                                   % check alignment to first image (rbmeanepi_* of first session)
        prefix  = 'u';
        sn      = varargin{1};
        target    = 1;
        vararginoptions(varargin(2:end), {'prefix','target'});
        
        cd(fullfile(imagingDir,subj_name{sn}));
        
        nSess = numel(subj_runs{sn});
        Q={};
        for k=1:nSess            
            for r=1:numel(subj_runs{sn}{k})
                runname = subj_runs{sn}{k}{r};
                for i=1:nTR
                    if use3D
                        Q{end+1}    = [fullfile(imagingDir,subj_name{sn},sprintf('%s',subj_sess{sn}{k}),[prefix, subj_name{sn},'_run_',runname,'_',num2str(i),'.nii'])];
                    else
                        Q{end+1}    = [fullfile(imagingDir,subj_name{sn},sprintf('%s',subj_sess{sn}{k}),[prefix, subj_name{sn},'_run_',runname,'.nii,',num2str(i)])];
                    end;
                end;
            end;            
        end;
        P{1} = [fullfile(imagingDir,subj_name{sn},sprintf('%s',subj_sess{sn}{target}),['rb' prefix 'meanepi_' subj_name{sn} '.nii'])];
        spmj_checksamealign(char(P),char(Q))
    case 'link_processed_images'                            % make symbolic link for sess2 imaging data into first session directory
        sn         = varargin{1};
        prefix     = 'u';
        sourceSess = 2;
        destSess = 1;
        for s=sn;
            % go to source session directory (in our case, sess2)
            for r = 1:numel(subj_runs{s}{sourceSess})
                runname = subj_runs{s}{sourceSess}{r};
                if use3D
                    source = fullfile(imagingDir,subj_name{s},sprintf('%s',subj_sess{s}{sourceSss}),...
                        [prefix, subj_name{s},'_run_',runname,'_*.nii']);
                    dest = fullfile(imagingDir,subj_name{s},sprintf('%s',subj_sess{s}{destSss}),...
                        [prefix, subj_name{s},'_run_',runname,'_*.nii']);
                else
                    source = fullfile(imagingDir,subj_name{s},sprintf('%s',subj_sess{s}{sourceSess}),...
                        [prefix, subj_name{s},'_run_',runname,'.nii']);
                    dest = fullfile(imagingDir,subj_name{s},sprintf('%s',subj_sess{s}{destSess}),...
                        [prefix, subj_name{s},'_run_',runname,'.nii']);
                end;
                % copy
                %suc=copyfile(source,dest);
                [suc,message] = system(sprintf('ln -s %s %s',source,dest));
                if ~suc
                    fprintf('Successfully created a symbolic link from: %s\n',source);
                    fprintf('\t to: %s\n',dest);
                else
                    warning('link failed: %s',message);
                end
            end;
        end;
    case 'make_maskImage'       % Make the mask image
        s  = varargin{1};
        prefix = 'u';
        vararginoptions(varargin(2:end), {'prefix'});
        
        for sn=s
            cd(fullfile(imagingDir,subj_name{sn}));
            
            % noskull mask for 1st level glm
            % -------------------------------
            nam={};
            nam{1}  = fullfile(imagingDir,subj_name{sn}, 'sess1', ['rb' prefix 'meanepi_' subj_name{sn} '.nii']);
            nam{2}  = fullfile(anatomicalDir, subj_name{sn}, ['c1' subj_name{sn}, '_anatomical.nii']);
            nam{3}  = fullfile(anatomicalDir, subj_name{sn}, ['c2' subj_name{sn}, '_anatomical.nii']);
            nam{4}  = fullfile(anatomicalDir, subj_name{sn}, ['c3' subj_name{sn}, '_anatomical.nii']);
            spm_imcalc_ui(nam, 'rmask_noskull.nii', 'i1>1 & (i2+i3+i4)>0.2')
            
            source = fullfile(imagingDir,subj_name{sn}, 'rmask_noskull.nii');
            dest = fullfile(anatomicalDir,subj_name{sn},'rmask_noskull.nii');
            movefile(source,dest);
            
            % gray matter mask for covariance estimation
            % ------------------------------------------
            nam={};
            nam{1}  = fullfile(imagingDir,subj_name{sn}, 'sess1', ['rb' prefix 'meanepi_' subj_name{sn} '.nii']);
            nam{2}  = fullfile(anatomicalDir, subj_name{sn}, ['c1' subj_name{sn}, '_anatomical.nii']);
            spm_imcalc_ui(nam, 'rmask_gray.nii', 'i1>1 & i2>0.4')
            
            source = fullfile(imagingDir,subj_name{sn}, 'rmask_gray.nii');
            dest = fullfile(anatomicalDir,subj_name{sn},'rmask_gray.nii');
            movefile(source,dest);
        end;
    case '_copy_processed_images'                                             % copy all preprocessed functional images to first session folder
        sn         = varargin{1};
        prefix     = 'ua';
        for s=sn;
            for sess=1:numel(subj_sess{s});
                for r = 1:numel(subj_runs{sn}{sess})
                    runname = subj_runs{sn}{sess}{r};
                    if use3D
                        source = fullfile(imagingDir,subj_name{s},sprintf('%s',subj_sess{s}{sess}),...
                            [prefix, subj_name{s},'_run',runname,'_*.nii']);
                    else
                        source = fullfile(imagingDir,subj_name{s},sprintf('%s',subj_sess{s}{sess}),...
                            [prefix, subj_name{s},'_run',runname,'.nii']);
                    end;
                    
                    source = fullfile(imagingDir,subj_name{s},sprintf('%s',subj_sess{s}{sess}),'*');
                    dest = fullfile(imagingDir,subj_name{s},'sess1_2');
                    dircheck(dest);
                    
                    copyfile(source,dest);
                    %movefile(source,dest);
                    fprintf('Copied file from: %s\n',source);
                    fprintf(' to: %s\n',dest);
                end;
            end;
        end;
    case '_move_processed_images'                                             % move all preprocessed functional images to first session folder and delete empty folders
        sn = varargin{1};
        
        for s=sn;
            for sess=1:numel(subj_sess{s});
                source = fullfile(imagingDir,subj_name{s},sprintf('%s',subj_sess{s}{sess}),'*');
                dest = fullfile(imagingDir,subj_name{s},'sess1');
                
                if strcmp(subj_sess{s}{sess},'sess1')==0
                    %movefile(source,dest);
                    fprintf('Moved file from: %s\n',source);
                    fprintf(' to: %s\n',dest);
                    
                    emptyfolder = fullfile(imagingDir,subj_name{s},sprintf('%s',subj_sess{s}{sess}));
                    
                    rmdir(emptyfolder,'s'); % sometime there is .DS_Store inside
                    fprintf('Deleted %s\n',emptyfolder);
                end
            end;
        end;
        
    case 'glm_do'                   % Do design and estimate in one go
        sn = varargin{1};
        glm = varargin{2};
        
        %close all
        %spm fmri
        %spm_jobman('initcfg');
        
        for g = glm
            for s = sn
                %glmtype = sprintf('glm_design_%d',g);
                glmtype = glmtypes{g};
                bmw1_uwo_imana(glmtype,s);
                bmw1_uwo_imana('glm_estimate',s,g);
                bmw1_uwo_imana('glm_contrast',s,g);
            end
            try
                bmw1_uwo_imana('ROI_timeseries_get',sn,g);
            catch message
                warning(message.message);
            end
        end
    case 'glm_design_6direction_fasthpf'    % Desing 1st-level GLM, with 'fast' option of spm12 with high-pass filter
        glm     = 1;
        sn      = varargin{1};
        Phrf    = [4 10];%[3.389  16  1   1    6  2.8101    0.1727]';
        prefix  = 'u';
        usefit  = 0;
        S       = [];
        
        % (0) Adjusting hrf parameters
        % ------------------------------------------------------------------------------------------
        % load fitted parameters (see ROI_timeseries_fit2)
        if usefit==1
            load(fullfile(regDir,sprintf('ROI_timeseries_fit2_%d_%s.mat',10,'12')));
            Para = getrow(T,T.SN==sn);
            Phrf = Para.P(Para.fit);
            clear T Ts
        end
        % (0) check spm version
        % ------------------------------------------------------------------------------------------
        if isempty(strfind(which('spm.m'),'spm12'));
            if ~(selectSPMversion('spm12'))
                warning('Failed to set path to spm12. Ending operation...');
                return;
            end
        end
        
        % (1) Setting experiment parameters
        % ------------------------------------------------------------------------------------------
        run = subj_runs{sn};
        
        try % load .dat file
            T = dload(fullfile(behaviourDir,[subj_name{sn}],['BimanualWrist_MR_UWO_', subj_name_behav{sn},'_scan.dat']));
            T.u_or_b    = T.Uni_or_Bi;
            T.hand      = T.Hand;
            T.targetAngleL = T.targetAngle_L;
            T.targetAngleR = T.targetAngle_R;
        catch % if there is no .dat file, use .tgt files
            T = [];
            for r = 1:nRun
                T_ = dload(fullfile(behaviourDir,[subj_name_behav{sn}],[subj_name_behav{sn},'_r',num2str(r),'.tgt']));
                T_.BN = repmat(r,size(T_.startSlice));
                T = addstruct(T,T_);
            end
        end
        
        DirL         = unique(T.targetAngleL);
        DirR         = unique(T.targetAngleR);
        nDirL        = length(DirL);
        nDirR        = length(DirL);
        
        %T.startTR = T.startSlice / 32 - (startTR)*ones(size(T.startSlice));    % the fact SPM start TR=0 is dealt here
        T.startTR = T.startTime/TR - numDummys;
        
        % Set delay and duration parameters according to individual behavioural data
        avrgRT      = nanmean(T.RT,1);      % reaction time
        avrgMT      = nanmean(T.MT,1);      % movement time
        holdTime    = nanmean(T.time2plan);
        delay       = 0;%1.5;%(holdTime+avrgRT)/TR; % Delay (for some reason delay is long)
        dur         = (avrgRT+avrgMT+holdTime)/TR;%avrgMT/TR*0.5;        % Duration
        
        % (2) Setting SPM GLM parameters
        % ------------------------------------------------------------------------------------------
        glmDir = fullfile(baseDir,glmName{glm});
        J.dir = {fullfile(glmDir,subj_name{sn})};
        if (~exist(J.dir{1},'dir'))
            mkdir(J.dir{1});
        end;
        J.timing.units = 'secs'; %'scans';
        J.timing.RT = 1.00; %2.72;
        J.timing.fmri_t = 16;
        J.timing.fmri_t0 = 1;
        
        % (3) Marking individual sessions in behavioural data
        % ------------------------------------------------------------------------------------------
        T.sess = zeros(length(T.BN),1);
        nRun = numel(subj_runs{sn}{1});
        for i=1:numel(run)
            sessR = [1:numel(subj_runs{sn}{i})] + (i-1)*nRun;
            nRun = numel(subj_runs{sn}{i});
            sessIdx = ismember(T.BN,sessR);
            T.sess(sessIdx) = i;
        end;
        
        sessCount = 1;
        runCount = 1;
        for k = 1:length(run) % session
            for r = 1:numel(run{k}) % run
                for i = 1:nTR-numDummys % image
                    if use3D
                        N{i} = [fullfile(imagingDir,subj_name{sn},sprintf('sess%d',k),[prefix,subj_name{sn},'_run_',run{k}{r},'_',num2str(i),'.nii'])];
                    else
                        N{i} = [fullfile(imagingDir,subj_name{sn},sprintf('sess%d',k),[prefix,subj_name{sn},'_run_',run{k}{r},'.nii,',num2str(i)])];
                    end;
                end;
                J.sess(sessCount).scans = N;
                J.sess(sessCount).cond = [];
                
                % get data for spesific block (=run)
                %blockIdx    = find(T.BN==runCount);% & T.sess==k);
                blockIdx    = find(T.BN==str2double(run{k}{r}));% & T.sess==k);
                R           = getrow(T,blockIdx);
                
                %--- Trial Type
                % Unimanua left
                for direction = 1:nDirL
                    index = (R.u_or_b==0)&(R.hand==0)&(R.targetAngleL==DirL(direction));
                    
                    J.sess(sessCount).cond(end+1).name = sprintf('UL%d/Run%d',direction, runCount);
                    J.sess(sessCount).cond(end).onset = [R.startTR(index) + delay];
                    J.sess(sessCount).cond(end).duration =  dur;
                    J.sess(sessCount).cond(end).tmod = 0;
                    J.sess(sessCount).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    J.sess(sessCount).cond(end).orth = 0;
                    
                    % saving run information in order to construct session
                    % specific contrasts outside of the glm
                    S_.sn           = sn;
                    S_.run          = runCount;%r;
                    S_.numEvents    = sum(index);
                    S_.hand         = 0; % left
                    S_.regType      = 1;
                    S_.u_or_b       = 0;
                    S_.dirL         = DirL(direction);
                    S_.dirR         = -1;
                    S_.movType      = direction;
                    
                    S               = addstruct(S,S_);
                end
                % Unimanual right
                for direction = 1:nDirR
                    index = (R.u_or_b==0)&(R.hand==1)&(R.targetAngleR==DirR(direction));
                    
                    J.sess(sessCount).cond(end+1).name = sprintf('UR%d/Run%d',direction, runCount);
                    J.sess(sessCount).cond(end).onset = [R.startTR(index) + delay];
                    J.sess(sessCount).cond(end).duration =  dur;
                    J.sess(sessCount).cond(end).tmod = 0;
                    J.sess(sessCount).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    J.sess(sessCount).cond(end).orth = 0;
                    
                    % saving run information in order to construct session
                    % specific contrasts outside of the glm
                    S_.sn           = sn;
                    S_.run          = runCount;%r;
                    S_.numEvents    = sum(index);
                    S_.hand         = 1; % right
                    S_.regType      = 1;  %
                    S_.u_or_b       = 0;
                    S_.dirL         = -1;
                    S_.dirR         = DirR(direction);
                    S_.movType      = nDirL + direction;
                    
                    S               = addstruct(S,S_);
                end
                % Bimanual
                for direction_ = 1:nDirL
                    for direction = 1:nDirR
                        index = (R.u_or_b==1)&(R.targetAngleL==DirL(direction_))&(R.targetAngleR==DirR(direction));
                        
                        J.sess(sessCount).cond(end+1).name = sprintf('B%d_%d/Run%d',direction_, direction, runCount);
                        J.sess(sessCount).cond(end).onset = [R.startTR(index) + delay];
                        J.sess(sessCount).cond(end).duration =  dur;
                        J.sess(sessCount).cond(end).tmod = 0;
                        J.sess(sessCount).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                        J.sess(sessCount).cond(end).orth = 0;
                        
                        % saving run information in order to construct session
                        % specific contrasts outside of the glm
                        S_.sn           = sn;
                        S_.run          = runCount;%r;
                        S_.numEvents    = sum(index);
                        S_.hand         = 2; % bimanual
                        S_.regType      = 1;  % Sequence regressor
                        S_.u_or_b       = 1;
                        S_.dirL         = DirL(direction_);
                        S_.dirR         = DirR(direction);
                        S_.movType      = (nDirL+nDirR) + (direction_-1)*nDirR + direction;
                        
                        S               = addstruct(S,S_);
                    end
                end
                %--- Maybe add error regrssor here (for another glm case)
                
                J.sess(sessCount).multi = {''};
                J.sess(sessCount).regress = struct('name', {}, 'val', {});
                J.sess(sessCount).multi_reg = {''};
                J.sess(sessCount).hpf = 128;
                sessCount = sessCount+1;
                runCount = runCount+1;
            end
        end;
        
        J.fact              = struct('name', {}, 'levels', {});
        J.bases.hrf.derivs  = [0 0];
        J.bases.hrf.params  = Phrf; % change hrf parameters
        J.volt              = 1;
        J.global            = 'None';
        J.mask              = {fullfile(anatomicalDir,subj_name{sn},'rmask_noSkull.nii')};
        J.cvi_mask          = {fullfile(anatomicalDir,subj_name{sn},'rmask_gray.nii')};
        J.cvi               =  'fast';
        J.mthresh           = 0.05;
        
        % (4) Specify 1st-level glm
        % ------------------------------------------------------------------------------------------
        spm_rwls_run_fmri_spec(J);
        
        save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','S');
        varargout = {glm};
    case 'glm_design_6direction_fastnohpf'  % Desing 1st-level GLM, with 'fast' option of spm12 without high-pass filter
        glm     = 2;
        sn      = varargin{1};
        Phrf    = [4 10];%[3.389 16  1   1    6  2.8101    0.1727]';
        prefix  = 'u';
        usefit  = 0;
        S       = [];
        
        % (0) Adjusting hrf parameters
        % ------------------------------------------------------------------------------------------
        % load fitted parameters (see ROI_timeseries_fit2)
        if usefit==1
            load(fullfile(regDir,sprintf('ROI_timeseries_fit2_%d_%s.mat',glm,'12')));
            Para = getrow(T,T.SN==sn);
            Phrf = Para.P(Para.fit);
            clear T Ts
        end
        
        % (0) check spm version
        % ------------------------------------------------------------------------------------------
        if isempty(strfind(which('spm.m'),'spm12'));
            if ~(selectSPMversion('spm12'))
                warning('Failed to set path to spm12. Ending operation...');
                return;
            end
        end
        
        % (1) Setting experiment parameters
        % ------------------------------------------------------------------------------------------
        run = subj_runs{sn};
        
        try % load .dat file
            T = dload(fullfile(behaviourDir,[subj_name{sn}],['BimanualWrist_MR_UWO_', subj_name_behav{sn},'_scan.dat']));
            T.u_or_b    = T.Uni_or_Bi;
            T.hand      = T.Hand;
            T.targetAngleL = T.targetAngle_L;
            T.targetAngleR = T.targetAngle_R;
        catch % if there is no .dat file, use .tgt files
            T = [];
            for r = 1:nRun
                T_ = dload(fullfile(behaviourDir,[subj_name_behav{sn}],[subj_name_behav{sn},'_r',num2str(r),'.tgt']));
                T_.BN = repmat(r,size(T_.startSlice));
                T = addstruct(T,T_);
            end
        end
        
        DirL         = unique(T.targetAngleL);
        DirR         = unique(T.targetAngleR);
        nDirL        = length(DirL);
        nDirR        = length(DirL);
        
        %T.startTR = T.startSlice / 32 - (startTR-1)*ones(size(T.startSlice));    % the fact SPM start TR=0 is dealt here
        T.startTR = T.startTime/TR-numDummys;
        
        % Set delay and duration parameters according to individual behavioural data
        avrgRT      = nanmean(T.RT,1);      % reaction time
        avrgMT      = nanmean(T.MT,1);      % movement time
        holdTime    = nanmean(T.time2plan);
        delay       = 0;%1.5%(holdTime+avrgRT)/TR; % Delay from task onset
        dur         = (avrgRT+avrgMT+holdTime)/TR;%avrgMT/TR*0.5;% Duration
        
        % (2) Setting SPM GLM parameters
        % ------------------------------------------------------------------------------------------
        glmDir = fullfile(baseDir,glmName{glm});
        J.dir = {fullfile(glmDir,subj_name{sn})};
        if (~exist(J.dir{1},'dir'))
            mkdir(J.dir{1});
        end;
        J.timing.units = 'secs';%'scans';
        J.timing.RT = 1.00;%2.72;
        J.timing.fmri_t = 16;
        J.timing.fmri_t0 = 1;
        
        % (3) Marking individual sessions in behavioural data
        % ------------------------------------------------------------------------------------------
        T.sess = zeros(length(T.BN),1);
        nRun = numel(subj_runs{sn}{1});
        for i=1:numel(run)
            sessR = [1:numel(subj_runs{sn}{i})] + (i-1)*nRun;
            nRun = numel(subj_runs{sn}{i});
            sessIdx = ismember(T.BN,sessR);
            T.sess(sessIdx) = i;
        end;
        
        sessCount = 1;
        runCount = 1;
        for k = 1:length(run) % session
            for r = 1:numel(run{k}) % run
                for i = 1:nTR-numDummys % image
                    if use3D
                        N{i} = [fullfile(imagingDir,subj_name{sn},sprintf('sess%d',k),[prefix,subj_name{sn},'_run_',run{k}{r},'_',num2str(i),'.nii'])];
                    else
                        N{i} = [fullfile(imagingDir,subj_name{sn},sprintf('sess%d',k),[prefix,subj_name{sn},'_run_',run{k}{r},'.nii,',num2str(i)])];
                    end;
                end;
                J.sess(sessCount).scans = N;
                J.sess(sessCount).cond = [];
                
                % get data for spesific block (=run)
                %blockIdx    = find(T.BN==runCount);% & T.sess==k);
                blockIdx    = find(T.BN==str2double(run{k}{r}));% & T.sess==k);
                R           = getrow(T,blockIdx);
                
                %--- Trial Type
                % Unimanua left
                for direction = 1:nDirL
                    index = (R.u_or_b==0)&(R.hand==0)&(R.targetAngleL==DirL(direction));
                    
                    J.sess(sessCount).cond(end+1).name = sprintf('UL%d/Run%d',direction, runCount);
                    J.sess(sessCount).cond(end).onset = [R.startTR(index) + delay];
                    J.sess(sessCount).cond(end).duration =  dur;
                    J.sess(sessCount).cond(end).tmod = 0;
                    J.sess(sessCount).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    J.sess(sessCount).cond(end).orth = 0;
                    
                    % saving run information in order to construct session
                    % specific contrasts outside of the glm
                    S_.sn           = sn;
                    S_.run          = runCount;%r;
                    S_.numEvents    = sum(index);
                    S_.hand         = 0; % left
                    S_.regType      = 1;
                    S_.u_or_b       = 0;
                    S_.dirL         = DirL(direction);
                    S_.dirR         = -1;
                    S_.movType      = direction;
                    S_.delayTime    = delay;
                    S_.duration     = dur;
                    
                    S               = addstruct(S,S_);
                end
                % Unimanual right
                for direction = 1:nDirR
                    index = (R.u_or_b==0)&(R.hand==1)&(R.targetAngleR==DirR(direction));
                    
                    J.sess(sessCount).cond(end+1).name = sprintf('UR%d/Run%d',direction, runCount);
                    J.sess(sessCount).cond(end).onset = [R.startTR(index) + delay];
                    J.sess(sessCount).cond(end).duration =  dur;
                    J.sess(sessCount).cond(end).tmod = 0;
                    J.sess(sessCount).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    J.sess(sessCount).cond(end).orth = 0;
                    
                    % saving run information in order to construct session
                    % specific contrasts outside of the glm
                    S_.sn           = sn;
                    S_.run          = runCount;%r;
                    S_.numEvents    = sum(index);
                    S_.hand         = 1; % right
                    S_.regType      = 1;  %
                    S_.u_or_b       = 0;
                    S_.dirL         = -1;
                    S_.dirR         = DirR(direction);
                    S_.movType      = nDirL + direction;
                    S_.delayTime    = delay;
                    S_.duration     = dur;
                    
                    S               = addstruct(S,S_);
                end
                % Bimanual
                for direction_ = 1:nDirL
                    for direction = 1:nDirR
                        index = (R.u_or_b==1)&(R.targetAngleL==DirL(direction_))&(R.targetAngleR==DirR(direction));
                        
                        J.sess(sessCount).cond(end+1).name = sprintf('B%d_%d/Run%d',direction_, direction, runCount);
                        J.sess(sessCount).cond(end).onset = [R.startTR(index) + delay];
                        J.sess(sessCount).cond(end).duration =  dur;
                        J.sess(sessCount).cond(end).tmod = 0;
                        J.sess(sessCount).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                        J.sess(sessCount).cond(end).orth = 0;
                        
                        % saving run information in order to construct session
                        % specific contrasts outside of the glm
                        S_.sn           = sn;
                        S_.run          = runCount;%r;
                        S_.numEvents    = sum(index);
                        S_.hand         = 2; % bimanual
                        S_.regType      = 1;  % Sequence regressor
                        S_.u_or_b       = 1;
                        S_.dirL         = DirL(direction_);
                        S_.dirR         = DirR(direction);
                        S_.movType      = (nDirL+nDirR) + (direction_-1)*nDirR + direction;
                        S_.delayTime    = delay;
                        S_.duration     = dur;
                        
                        S               = addstruct(S,S_);
                    end
                end
                %--- Maybe add error regrssor here (for another glm case)
                
                J.sess(sessCount).multi = {''};
                J.sess(sessCount).regress = struct('name', {}, 'val', {});
                J.sess(sessCount).multi_reg = {''};
                J.sess(sessCount).hpf = inf;%128;
                sessCount = sessCount+1;
                runCount = runCount+1;
            end
        end;
        
        J.fact              = struct('name', {}, 'levels', {});
        J.bases.hrf.derivs  = [0 0];
        J.bases.hrf.params  = Phrf; % change hrf parameters
        J.volt              = 1;
        J.global            = 'None';
        J.mask              = {fullfile(anatomicalDir,subj_name{sn},'rmask_noSkull.nii')};
        J.cvi_mask          = {fullfile(anatomicalDir,subj_name{sn},'rmask_gray.nii')};
        J.cvi               =  'fast';
        J.mthresh           = 0.05;
        
        % (4) Specify 1st-level glm
        % ------------------------------------------------------------------------------------------
        spm_rwls_run_fmri_spec(J);
        
        save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','S');
        varargout = {glm};
    case 'glm_design_6direction_rwlshpf'    % Desing 1st-level GLM, with 'rwls' option of spm12 with high-pass filter
        glm     = 3;
        sn      = varargin{1};
        Phrf    = [4 10];%[3.389 16  1   1    6  2.8101    0.1727]';
        prefix  = 'u';
        usefit  = 0;
        S       = [];
        
        % (0) Adjusting hrf parameters
        % ------------------------------------------------------------------------------------------
        % load fitted parameters (see ROI_timeseries_fit2)
        if usefit==1
            load(fullfile(regDir,sprintf('ROI_timeseries_fit2_%d_%s.mat',10,'12')));
            Para = getrow(T,T.SN==sn);
            Phrf = Para.P(Para.fit);
            clear T Ts
        end
        
        % (0) check spm version
        % ------------------------------------------------------------------------------------------
        if isempty(strfind(which('spm.m'),'spm12'));
            if ~(selectSPMversion('spm12'))
                warning('Failed to set path to spm12. Ending operation...');
                return;
            end
        end
        
        % (1) Setting experiment parameters
        % ------------------------------------------------------------------------------------------
        run = subj_runs{sn};
        
        try % load .dat file
            T = dload(fullfile(behaviourDir,[subj_name{sn}],['BimanualWrist_MR_UWO_', subj_name_behav{sn},'_scan.dat']));
            T.u_or_b    = T.Uni_or_Bi;
            T.hand      = T.Hand;
            T.targetAngleL = T.targetAngle_L;
            T.targetAngleR = T.targetAngle_R;
        catch % if there is no .dat file, use .tgt files
            T = [];
            for r = 1:nRun
                T_ = dload(fullfile(behaviourDir,[subj_name_behav{sn}],[subj_name_behav{sn},'_r',num2str(r),'.tgt']));
                T_.BN = repmat(r,size(T_.startSlice));
                T = addstruct(T,T_);
            end
        end
        
        DirL         = unique(T.targetAngleL);
        DirR         = unique(T.targetAngleR);
        nDirL        = length(DirL);
        nDirR        = length(DirL);
        
        %T.startTR = T.startSlice / 32 - (startTR)*ones(size(T.startSlice));    % the fact SPM start TR=0 is dealt here
        T.startTR = T.startTime/TR-numDummys;
        
        % Set delay and duration parameters according to individual behavioural data
        avrgRT      = nanmean(T.RT,1);      % reaction time
        avrgMT      = nanmean(T.MT,1);      % movement time
        holdTime    = nanmean(T.time2plan);
        delay       = 0;%1.5;%(holdTime+avrgRT)/TR; % Delay (for some reason delay is long)
        dur         = (avrgRT+avrgMT+holdTime)/TR;%avrgMT/TR*0.5;        % Duration
        
        % (2) Setting SPM GLM parameters
        % ------------------------------------------------------------------------------------------
        glmDir = fullfile(baseDir,glmName{glm});
        J.dir = {fullfile(glmDir,subj_name{sn})};
        if (~exist(J.dir{1},'dir'))
            mkdir(J.dir{1});
        end;
        J.timing.units = 'secs';%'scans';
        J.timing.RT = 1.00;%2.72;
        J.timing.fmri_t = 16;
        J.timing.fmri_t0 = 1;
        
        % (3) Marking individual sessions in behavioural data
        % ------------------------------------------------------------------------------------------
        T.sess = zeros(length(T.BN),1);
        nRun = numel(subj_runs{sn}{1});
        for i=1:numel(run)
            sessR = [1:numel(subj_runs{sn}{i})] + (i-1)*nRun;
            nRun = numel(subj_runs{sn}{i});
            sessIdx = ismember(T.BN,sessR);
            T.sess(sessIdx) = i;
        end;
        
        sessCount = 1;
        runCount = 1;
        for k = 1:length(run) % session
            for r = 1:numel(run{k}) % run
                for i = 1:nTR-numDummys % image
                    if use3D
                        N{i} = [fullfile(imagingDir,subj_name{sn},sprintf('sess%d',k),[prefix,subj_name{sn},'_run_',run{k}{r},'_',num2str(i),'.nii'])];
                    else
                        N{i} = [fullfile(imagingDir,subj_name{sn},sprintf('sess%d',k),[prefix,subj_name{sn},'_run_',run{k}{r},'.nii,',num2str(i)])];
                    end;
                end;
                J.sess(sessCount).scans = N;
                J.sess(sessCount).cond = [];
                
                % get data for spesific block (=run)
                %blockIdx    = find(T.BN==runCount);% & T.sess==k);
                blockIdx    = find(T.BN==str2double(run{k}{r}));% & T.sess==k);
                R           = getrow(T,blockIdx);
                
                %--- Trial Type
                % Unimanua left
                for direction = 1:nDirL
                    index = (R.u_or_b==0)&(R.hand==0)&(R.targetAngleL==DirL(direction));
                    
                    J.sess(sessCount).cond(end+1).name = sprintf('UL%d/Run%d',direction, runCount);
                    J.sess(sessCount).cond(end).onset = [R.startTR(index) + delay];
                    J.sess(sessCount).cond(end).duration =  dur;
                    J.sess(sessCount).cond(end).tmod = 0;
                    J.sess(sessCount).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    J.sess(sessCount).cond(end).orth = 0;
                    
                    % saving run information in order to construct session
                    % specific contrasts outside of the glm
                    S_.sn           = sn;
                    S_.run          = runCount;%r;
                    S_.numEvents    = sum(index);
                    S_.hand         = 0; % left
                    S_.regType      = 1;
                    S_.u_or_b       = 0;
                    S_.dirL         = DirL(direction);
                    S_.dirR         = -1;
                    S_.movType      = direction;
                    
                    S               = addstruct(S,S_);
                end
                % Unimanual right
                for direction = 1:nDirR
                    index = (R.u_or_b==0)&(R.hand==1)&(R.targetAngleR==DirR(direction));
                    
                    J.sess(sessCount).cond(end+1).name = sprintf('UR%d/Run%d',direction, runCount);
                    J.sess(sessCount).cond(end).onset = [R.startTR(index) + delay];
                    J.sess(sessCount).cond(end).duration =  dur;
                    J.sess(sessCount).cond(end).tmod = 0;
                    J.sess(sessCount).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    J.sess(sessCount).cond(end).orth = 0;
                    
                    % saving run information in order to construct session
                    % specific contrasts outside of the glm
                    S_.sn           = sn;
                    S_.run          = runCount;%r;
                    S_.numEvents    = sum(index);
                    S_.hand         = 1; % right
                    S_.regType      = 1;  %
                    S_.u_or_b       = 0;
                    S_.dirL         = -1;
                    S_.dirR         = DirR(direction);
                    S_.movType      = nDirL + direction;
                    
                    S               = addstruct(S,S_);
                end
                % Bimanual
                for direction_ = 1:nDirL
                    for direction = 1:nDirR
                        index = (R.u_or_b==1)&(R.targetAngleL==DirL(direction_))&(R.targetAngleR==DirR(direction));
                        
                        J.sess(sessCount).cond(end+1).name = sprintf('B%d_%d/Run%d',direction_, direction, runCount);
                        J.sess(sessCount).cond(end).onset = [R.startTR(index) + delay];
                        J.sess(sessCount).cond(end).duration =  dur;
                        J.sess(sessCount).cond(end).tmod = 0;
                        J.sess(sessCount).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                        J.sess(sessCount).cond(end).orth = 0;
                        
                        % saving run information in order to construct session
                        % specific contrasts outside of the glm
                        S_.sn           = sn;
                        S_.run          = runCount;%r;
                        S_.numEvents    = sum(index);
                        S_.hand         = 2; % bimanual
                        S_.regType      = 1;  % Sequence regressor
                        S_.u_or_b       = 1;
                        S_.dirL         = DirL(direction_);
                        S_.dirR         = DirR(direction);
                        S_.movType      = (nDirL+nDirR) + (direction_-1)*nDirR + direction;
                        
                        S               = addstruct(S,S_);
                    end
                end
                %--- Maybe add error regrssor here (for another glm case)
                
                J.sess(sessCount).multi = {''};
                J.sess(sessCount).regress = struct('name', {}, 'val', {});
                J.sess(sessCount).multi_reg = {''};
                J.sess(sessCount).hpf = 128;
                sessCount = sessCount+1;
                runCount = runCount+1;
            end
        end;
        
        J.fact              = struct('name', {}, 'levels', {});
        J.bases.hrf.derivs  = [0 0];
        J.bases.hrf.params  = Phrf; % change hrf parameters
        J.volt              = 1;
        J.global            = 'None';
        J.mask              = {fullfile(anatomicalDir,subj_name{sn},'rmask_noSkull.nii')};
        J.cvi_mask          = {fullfile(anatomicalDir,subj_name{sn},'rmask_gray.nii')};
        J.cvi               =  'wls';
        J.mthresh           = 0.05;
        
        % (4) Specify 1st-level glm
        % ------------------------------------------------------------------------------------------
        spm_rwls_run_fmri_spec(J);
        
        save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','S');
        varargout = {glm};
    case 'glm_design_6direction_rwlsnohpf'  % Desing 1st-level GLM, with 'rwls' option of spm12 without high-pass filter
        glm     = 4;
        sn      = varargin{1};
        Phrf    = [4 10];%[3.389 16  1   1    6  2.8101    0.1727]';
        prefix  = 'u';
        usefit  = 0;
        S       = [];
        
        % (0) Adjusting hrf parameters
        % ------------------------------------------------------------------------------------------
        % load fitted parameters (see ROI_timeseries_fit2)
        if usefit==1
            load(fullfile(regDir,sprintf('ROI_timeseries_fit2_%d_%s.mat',10,'12')));
            Para = getrow(T,T.SN==sn);
            Phrf = Para.P(Para.fit);
            clear T Ts
        end
        
        % (0) check spm version
        % ------------------------------------------------------------------------------------------
        if isempty(strfind(which('spm.m'),'spm12'));
            if ~(selectSPMversion('spm12'))
                warning('Failed to set path to spm12. Ending operation...');
                return;
            end
        end
        
        % (1) Setting experiment parameters
        % ------------------------------------------------------------------------------------------
        run = subj_runs{sn};
        
        try % load .dat file
            T = dload(fullfile(behaviourDir,[subj_name{sn}],['BimanualWrist_MR_UWO_', subj_name_behav{sn},'_scan.dat']));
            T.u_or_b    = T.Uni_or_Bi;
            T.hand      = T.Hand;
            T.targetAngleL = T.targetAngle_L;
            T.targetAngleR = T.targetAngle_R;
        catch % if there is no .dat file, use .tgt files
            T = [];
            for r = 1:nRun
                T_ = dload(fullfile(behaviourDir,[subj_name_behav{sn}],[subj_name_behav{sn},'_r',num2str(r),'.tgt']));
                T_.BN = repmat(r,size(T_.startSlice));
                T = addstruct(T,T_);
            end
        end
        
        DirL         = unique(T.targetAngleL);
        DirR         = unique(T.targetAngleR);
        nDirL        = length(DirL);
        nDirR        = length(DirL);
        
        %T.startTR = T.startSlice / 32 - (startTR)*ones(size(T.startSlice));    % the fact SPM start TR=0 is dealt here
        T.startTR = T.startTime/TR-numDummys;
        
        % Set delay and duration parameters according to individual behavioural data
        avrgRT      = nanmean(T.RT,1);      % reaction time
        avrgMT      = nanmean(T.MT,1);      % movement time
        holdTime    = nanmean(T.time2plan);
        delay       = 0;%;1.5;%(holdTime+avrgRT)/TR; % Delay (for some reason delay is long)
        dur         = (avrgRT+avrgMT+holdTime)/TR;%avrgMT/TR*0.5;        % Duration
        
        % (2) Setting SPM GLM parameters
        % ------------------------------------------------------------------------------------------
        glmDir = fullfile(baseDir,glmName{glm});
        J.dir = {fullfile(glmDir,subj_name{sn})};
        if (~exist(J.dir{1},'dir'))
            mkdir(J.dir{1});
        end;
        J.timing.units = 'secs';%'scans';
        J.timing.RT = 1.00;%2.72;
        J.timing.fmri_t = 16;
        J.timing.fmri_t0 = 1;
        
        % (3) Marking individual sessions in behavioural data
        % ------------------------------------------------------------------------------------------
        T.sess = zeros(length(T.BN),1);
        nRun = numel(subj_runs{sn}{1});
        for i=1:numel(run)
            sessR = [1:numel(subj_runs{sn}{i})] + (i-1)*nRun;
            nRun = numel(subj_runs{sn}{i});
            sessIdx = ismember(T.BN,sessR);
            T.sess(sessIdx) = i;
        end;
        
        sessCount = 1;
        runCount = 1;
        for k = 1:length(run) % session
            for r = 1:numel(run{k}) % run
                for i = 1:nTR-numDummys % image
                    if use3D
                        N{i} = [fullfile(imagingDir,subj_name{sn},sprintf('sess%d',k),[prefix,subj_name{sn},'_run_',run{k}{r},'_',num2str(i),'.nii'])];
                    else
                        N{i} = [fullfile(imagingDir,subj_name{sn},sprintf('sess%d',k),[prefix,subj_name{sn},'_run_',run{k}{r},'.nii,',num2str(i)])];
                    end;
                end;
                J.sess(sessCount).scans = N;
                J.sess(sessCount).cond = [];
                
                % get data for spesific block (=run)
                %blockIdx    = find(T.BN==runCount);% & T.sess==k);
                blockIdx    = find(T.BN==str2double(run{k}{r}));% & T.sess==k);
                R           = getrow(T,blockIdx);
                
                %--- Trial Type
                % Unimanua left
                for direction = 1:nDirL
                    index = (R.u_or_b==0)&(R.hand==0)&(R.targetAngleL==DirL(direction));
                    
                    J.sess(sessCount).cond(end+1).name = sprintf('UL%d/Run%d',direction, runCount);
                    J.sess(sessCount).cond(end).onset = [R.startTR(index) + delay];
                    J.sess(sessCount).cond(end).duration =  dur;
                    J.sess(sessCount).cond(end).tmod = 0;
                    J.sess(sessCount).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    J.sess(sessCount).cond(end).orth = 0;
                    
                    % saving run information in order to construct session
                    % specific contrasts outside of the glm
                    S_.sn           = sn;
                    S_.run          = runCount;%r;
                    S_.numEvents    = sum(index);
                    S_.hand         = 0; % left
                    S_.regType      = 1;
                    S_.u_or_b       = 0;
                    S_.dirL         = DirL(direction);
                    S_.dirR         = -1;
                    S_.movType      = direction;
                    
                    S               = addstruct(S,S_);
                end
                % Unimanual right
                for direction = 1:nDirR
                    index = (R.u_or_b==0)&(R.hand==1)&(R.targetAngleR==DirR(direction));
                    
                    J.sess(sessCount).cond(end+1).name = sprintf('UR%d/Run%d',direction, runCount);
                    J.sess(sessCount).cond(end).onset = [R.startTR(index) + delay];
                    J.sess(sessCount).cond(end).duration =  dur;
                    J.sess(sessCount).cond(end).tmod = 0;
                    J.sess(sessCount).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    J.sess(sessCount).cond(end).orth = 0;
                    
                    % saving run information in order to construct session
                    % specific contrasts outside of the glm
                    S_.sn           = sn;
                    S_.run          = runCount;%r;
                    S_.numEvents    = sum(index);
                    S_.hand         = 1; % right
                    S_.regType      = 1;  %
                    S_.u_or_b       = 0;
                    S_.dirL         = -1;
                    S_.dirR         = DirR(direction);
                    S_.movType      = nDirL + direction;
                    
                    S               = addstruct(S,S_);
                end
                % Bimanual
                for direction_ = 1:nDirL
                    for direction = 1:nDirR
                        index = (R.u_or_b==1)&(R.targetAngleL==DirL(direction_))&(R.targetAngleR==DirR(direction));
                        
                        J.sess(sessCount).cond(end+1).name = sprintf('B%d_%d/Run%d',direction_, direction, runCount);
                        J.sess(sessCount).cond(end).onset = [R.startTR(index) + delay];
                        J.sess(sessCount).cond(end).duration =  dur;
                        J.sess(sessCount).cond(end).tmod = 0;
                        J.sess(sessCount).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                        J.sess(sessCount).cond(end).orth = 0;
                        
                        % saving run information in order to construct session
                        % specific contrasts outside of the glm
                        S_.sn           = sn;
                        S_.run          = runCount;%r;
                        S_.numEvents    = sum(index);
                        S_.hand         = 2; % bimanual
                        S_.regType      = 1;  % Sequence regressor
                        S_.u_or_b       = 1;
                        S_.dirL         = DirL(direction_);
                        S_.dirR         = DirR(direction);
                        S_.movType      = (nDirL+nDirR) + (direction_-1)*nDirR + direction;
                        
                        S               = addstruct(S,S_);
                    end
                end
                %--- Maybe add error regrssor here (for another glm case)
                
                J.sess(sessCount).multi = {''};
                J.sess(sessCount).regress = struct('name', {}, 'val', {});
                J.sess(sessCount).multi_reg = {''};
                J.sess(sessCount).hpf = inf;%128;
                sessCount = sessCount+1;
                runCount = runCount+1;
            end
        end;
        
        J.fact              = struct('name', {}, 'levels', {});
        J.bases.hrf.derivs  = [0 0];
        J.bases.hrf.params  = Phrf; % change hrf parameters
        J.volt              = 1;
        J.global            = 'None';
        J.mask              = {fullfile(anatomicalDir,subj_name{sn},'rmask_noSkull.nii')};
        J.cvi_mask          = {fullfile(anatomicalDir,subj_name{sn},'rmask_gray.nii')};
        J.cvi               =  'wls';
        J.mthresh           = 0.05;
        
        % (4) Specify 1st-level glm
        % ------------------------------------------------------------------------------------------
        spm_rwls_run_fmri_spec(J);
        
        save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','S');
        varargout = {glm};
    case 'glm_estimate'             % Estimate 1st-level GLM
        sn=varargin{1};
        glmType=varargin{2};
        
        switch glmType
            case {1,2,3,4}
                % Check spm version
                if isempty(strfind(which('spm.m'),'spm12'));
                    if ~(selectSPMversion('spm12'))
                        warning('Failed to set path to spm12. Ending operation...');
                        return;
                    end
                end
                glmDir = fullfile(baseDir,glmName{glmType},subj_name{sn});
                cd(glmDir);
                load('SPM.mat');
                spm_rwls_spm(SPM);
        end
        % for checking
        % spm_rwls_resstats(SPM)
    case 'glm_contrast'                     % Make 1st-level GLM constast
        sn = varargin{1};
        glm = varargin{2};
        
        switch (glm)
            case {1,2,3,4}
                type = 'glm_contrast_6direction';
            otherwise
                error('invalid value for glm.');
        end
        bmw1_uwo_imana(type,sn,glm);
    case 'glm_contrast_6direction'          % Make 1st-level GLM contrast
        % 10 contrasts including t and F
        
        sn      = varargin{1};
        glmType = varargin{2};
        glmDir  = fullfile(baseDir,glmName{glmType});
        nConditions = 48;
        
        for s = sn
            %             C = [];
            %             Names = [];
            cd(fullfile(glmDir, subj_name{s}));
            load('SPM.mat');
            T   = load('SPM_info.mat');
            SPM = rmfield(SPM,'xCon');
            nrun    = numel(SPM.nscan);
            tt      = repmat([1:nConditions],1,nrun);
            %sess    = unique(T.sess)';
            
            % calc bimanual contrast for intrinsic,extrinsic,asymmetric
            % movements
            extrinsic = [0:60:300];
            intrinsic = [180,120,60,0,300,240];
            for t=1:length(T.hand)
                if T.u_or_b(t,1)==1
                    idx = find(extrinsic==T.dirL(t,1));
                    if T.dirR(t,1)==extrinsic(idx)
                        BiCat(t,1) = 1;
                    elseif T.dirR(t,1)==intrinsic(idx)
                        BiCat(t,1) = 2;
                    else
                        BiCat(t,1) = 3;
                    end
                else
                    BiCat(t,1) = 0;
                end
            end
            
            
            % Note: the rest phase is not explicitly included in design matrix
            
            % (1) t-Contrast movement against rest for each movement type
            % (1:48)
            for c = 1:nConditions
                con = zeros(1,size(SPM.xX.X,2));
                con(tt==c) = 1;
                con = con/sum(con);
                SPM.xCon(c) = spm_FcUtil('Set',sprintf('MovType%d',c), 'T', 'c',con',SPM.xX.xKXs);
            end
            
            % (2) t-Contrast movement for any movement against rest (49)
            con=zeros(1,size(SPM.xX.X,2));
            con(1:length(tt))=1;
            con=con/sum(con);
            SPM.xCon(end+1)=spm_FcUtil('Set','Movement', 'T', 'c',con',SPM.xX.xKXs);
            
            % (3) F-Contrast every movement against rest (F-test) (50)
            nvar = size(SPM.xX.X,2);
            con = zeros(nConditions,nvar);
            for i = 1:nConditions
                idx = tt==i;
                con(i,idx)=1/sum(idx);
            end;
            SPM.xCon(end+1)=spm_FcUtil('Set','Seq', 'F', 'c',con',SPM.xX.xKXs);
            
            % (4) t-Contrast unimanual left movements against rest (51)
            con=zeros(1,size(SPM.xX.X,2));
            con(T.hand==0 & T.u_or_b==0) = 1;
            con = con/sum(con);
            SPM.xCon(end+1) = spm_FcUtil('Set','Uni Left', 'T', 'c',con',SPM.xX.xKXs);
            
            % (5) t-Contrast unimanual right movements against rest (52)
            con=zeros(1,size(SPM.xX.X,2));
            con(T.hand==1 & T.u_or_b==0) = 1;
            con = con/sum(con);
            SPM.xCon(end+1) = spm_FcUtil('Set','Uni Right', 'T', 'c',con',SPM.xX.xKXs);
            
            % (6) t-Contrast bimanual movements against rest (53)
            con=zeros(1,size(SPM.xX.X,2));
            con(T.u_or_b==1) = 1;
            con = con/sum(con);
            SPM.xCon(end+1) = spm_FcUtil('Set','Bimanual', 'T', 'c',con',SPM.xX.xKXs);
            
            % (7) t-Contrast: uni L>R (54)
            con=zeros(1,size(SPM.xX.X,2));
            con(T.u_or_b==0&T.hand==0) = 1;
            con(T.u_or_b==0&T.hand==1) = -1;
            con(con==1) = con(con==1)/sum(con==1);
            con(con==-1) = con(con==-1)/sum(con==-1);
            SPM.xCon(end+1) = spm_FcUtil('Set','UniL>UniR', 'T', 'c',con',SPM.xX.xKXs);
            
            % (8) t-Contrast: bi>uni (55)
            con=zeros(1,size(SPM.xX.X,2));
            con(T.u_or_b==1) = 1;
            con(T.u_or_b==0) = -1;
            con(con==1) = con(con==1)/sum(con==1);
            con(con==-1) = con(con==-1)/sum(con==-1);
            SPM.xCon(end+1) = spm_FcUtil('Set','Bi>Uni', 'T', 'c',con',SPM.xX.xKXs);
            
            % (9) t-Contrast: bi-asymmetric>bi-symmetric-extrinsic (56)
            con=zeros(1,size(SPM.xX.X,2));
            con(BiCat==3) = 1;
            con(BiCat==1) = -1;
            con(con==1) = con(con==1)/sum(con==1);
            con(con==-1) = con(con==-1)/sum(con==-1);
            SPM.xCon(end+1) = spm_FcUtil('Set','BiAsym>BiSymEx', 'T', 'c',con',SPM.xX.xKXs);
            
            % (10) t-Contrast: bi-asymmetric>bi-symmetric-intrinsic (57)
            con=zeros(1,size(SPM.xX.X,2));
            con(BiCat==3) = 1;
            con(BiCat==2) = -1;
            con(con==1) = con(con==1)/sum(con==1);
            con(con==-1) = con(con==-1)/sum(con==-1);
            SPM.xCon(end+1) = spm_FcUtil('Set','BiAsym>BiSymIn', 'T', 'c',con',SPM.xX.xKXs);
            
            %             C(end+1,:) = con;
            %             Names{end+1} = 'PassiveMultiDigit';
            %             Nums = 1:length(Names);
            %
            %             %____do the constrasts
            %             SPM = spmj_setContrasts(C,Names,Nums,SPM,optimize); % use optimal contrast
            %             SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
            
            
            % (End) Generate contrasts
            SPM=spm_contrasts(SPM,[1:length(SPM.xCon)]);
            save('SPM.mat','SPM', '-v7.3');
        end
        
    case 'surf_freesurfer'              % SURFACE PREPROCESS 1: Run recon-all on freesurfer
        sn = varargin{1};
        %mysetEnv(freesurferDir);
        %freesurfer_reconall(freesurferDir,subj_name{i},fullfile(anatomicalDir,subj_name{i},[subj_name{i} '_anatomical.nii']));
        fprintf('When this did not work, copy&paste following in terminal...\n');
        fprintf('--------------------------------------------------------------\n\n')
        fprintf('export FREESURFER_HOME=''/Applications/freesurfer''\n');
        fprintf('export ANATOMICAL_DIR=''%s''\n',anatomicalDir);
        fprintf('source $FREESURFER_HOME/SetUpFreeSurfer.sh\n');
        fprintf('source $FREESURFER_HOME/FreeSurferEnv.sh\n');
        fprintf('export SUBJECTS_DIR=''%s''\n',freesurferDir);
        fprintf('recon-all -s %s -i $ANATOMICAL_DIR/%s/%s_anatomical.nii -all -cw256\n\n',subj_name{sn},subj_name{sn},subj_name{sn});
    case 'surf_xhemireg'                % SURFACE PREPROCESS 2: cross-Register surfaces left / right hem
        sn=varargin{1};
        %for i=sn
        %freesurfer_registerXhem({subj_name{i}},freesurferDir,'hemisphere',[1 2]); % For debug... [1 2] orig
        %end;
        if ~exist(fullfile(freesurferDir,'fsaverage_sym'),'dir');
            command = sprintf('ln -s /Applications/freesurfer/subjects/fsaverage_sym %s',freesurferDir);
            system(command);
            disp('Setup symbolic link for ''fsaverage_sym''!');
        end
        fprintf('When this did not work, copy&paste following in terminal...\n');
        fprintf('--------------------------------------------------------------\n\n')
        fprintf('export FREESURFER_HOME=''/Applications/freesurfer''\n');
        fprintf('export ANATOMICAL_DIR=''%s''\n',anatomicalDir');
        fprintf('source $FREESURFER_HOME/SetUpFreeSurfer.sh\n');
        fprintf('source $FREESURFER_HOME/FreeSurferEnv.sh\n');
        fprintf('export SUBJECTS_DIR=''%s''\n',freesurferDir);
        fprintf('surfreg --s %s --t fsaverage_sym --lh\n',subj_name{sn});
        fprintf('xhemireg --s %s\n',subj_name{sn});
        fprintf('surfreg --s %s --t fsaverage_sym --lh --xhemi\n\n',subj_name{sn});
    case 'surf_map_ico'                 % SURFACE PREPROCESS 3: Align to the new atlas surface (map icosahedron)
        sn=varargin{1};
        for i=sn
            freesurfer_mapicosahedron_xhem(subj_name{i},freesurferDir,'smoothing',1,'hemisphere',[1:2]);
        end;
    case 'surf_make_caret'              % SURFACE PREPROCESS 4: Translate into caret format
        sn=varargin{1};
        for i=sn
            caret_importfreesurfer(['x' subj_name{i}],freesurferDir,caretDir);
        end;
        
    case 'MVA_search_surf'          % Define surface-based search light using volume data (efficient!)
        sn = varargin{1};
        glm = varargin{2};
        atlas = 1;
        refDir= fullfile(baseDir,glmName{glm});
        
        for s=sn
            for h=1:2
                caret_subjDIR = fullfile(caretDir,[atlasA{atlas},subj_name{s}],hemName{h});
                coord_pial= caret_load(fullfile(caret_subjDIR, [hem{h} '.PIAL.coord']));
                coord_white= caret_load(fullfile(caret_subjDIR, [hem{h} '.WHITE.coord']));
                topo= caret_load(fullfile(caret_subjDIR, [hem{h} '.CLOSED.topo']));
                surf(h).c1=coord_white.data';
                surf(h).c2=coord_pial.data';
                surf(h).f=topo.data';
            end;
            clear coord_pial coord_white coord_caret;
            volDef= spm_vol(fullfile(refDir, subj_name{s},'mask.nii'));
            volDef.mask=spm_read_vols(volDef);
            [LI,voxmin,voxmax,voxel,node,surfindx,depth,radvox]= lmva_voxelselection_surf(surf, [6 160],volDef);
            LI=LI';
            save(fullfile(refDir,subj_name{s}, 'vol_roi_160vox.mat'), 'LI','voxmin','voxmax','voxel','radvox');
            save(fullfile(refDir,subj_name{s}, 'vol_surf.mat'), 'voxel','node','surfindx','depth');
            
            V=volDef;
            X=zeros(V.dim);
            depth=depth+1;
            depth(isnan(depth))=1;
            X(voxel)=depth;
            V=volDef;
            V.dt=[4 0];
            V.pinfo=[3/100 0 0]';
            V.fname=fullfile(refDir,subj_name{s},'vol_roi_depth.nii');
            spm_write_vol(V,X);
        end;
        varargout = {};
    case 'MVA_dist_var_raw'         % Calc G and pattern distance (prewhiten data using covariance matrix estimated from raw time series data)
        sn              = varargin{1};
        glm             = 1;
        SVD             = 0;
        split 			= 0; % split data into two and estimate G and dist from each half
        suffix          = '';
        
        vararginoptions({varargin{2:end}},...
            {'glm', 'suffix', 'SVD'});
        
        ldafunction = @sh1_calcDistVar_raw;
        
        lG 		= nConditions^2; % length of G matrix
        ldist 	= nConditions*(nConditions-1)/2; % length of distance
        
        home = cd;
        for s = sn
            glmDirSubj=fullfile(baseDir,glmName{glm},subj_name{s});
            cd(glmDirSubj);
            nii_out = {};
            
            % Define output files
            % dist -> 4D
            for elem = 1:ldist
                nii_out{end+1} = fullfile(glmDirSubj,[subj_name{s},'_dist_raw', suffix, '.nii,', num2str(elem)]);
            end
            % sigma -> 4D
            for elem = 1:lG
                nii_out{end+1} = fullfile(glmDirSubj,[subj_name{s},'_sigma_raw', suffix, '.nii,', num2str(elem)]);
            end
            
            % load search light definition, etc.
            load(fullfile(glmDirSubj,'SPM.mat'));
            T   = load(fullfile(glmDirSubj,'SPM_info.mat'));
            %SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
            SL  = load(fullfile(glmDirSubj,['vol_roi_160vox.mat']));
            idx = find(T.regType==1);
            
            params = {SPM, T.run(idx)', T.movType(idx)',...
                'split', split,...
                'SVD',SVD};
            
            % run searchlight analysis
            lmva_spm(SL,SPM.xY.P,nii_out,ldafunction,'params',params,'isNP',1);
        end % s
        
        cd(home);
        varargout = {nii_out};
    case 'MVA_dist_var_raw_summary' % searchlight analysis (only save summary value)
        sn              = varargin{1};
        glm             = 1;
        SVD             = 0;
        split 			= 0; % split data into two and estimate G and dist from each half
        suffix          = '';
        nDirL           = 8;
        nDirR           = 8;
        
        vararginoptions({varargin{2:end}},...
            {'glm', 'suffix', 'SVD','nDirL','nDirR'});
        
        ldafunction = @bmw1_calcDistVar_raw;
        
        suffix = sprintf('_glm%d',glm);
        home = cd;
        for s = sn
            glmDirSubj=fullfile(baseDir,glmName{glm},subj_name{s});
            cd(glmDirSubj);
            nii_out = {};
            
            % Define output files
            % mean distance -> 4D
            for elem = 1:6
                nii_out{end+1} = fullfile(glmDirSubj,[subj_name{s},'_dist_raw', suffix, '.nii,', num2str(elem)]);
            end
            
            % load search light definition, etc.
            load(fullfile(glmDirSubj,'SPM.mat'));
            T   = load(fullfile(glmDirSubj,'SPM_info.mat'));
            %SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
            SL  = load(fullfile(glmDirSubj,['vol_roi_160vox.mat']));
            idx = find(T.regType==1);
            
            params = {SPM, T.run(idx)', T.movType(idx)',...
                'split', split,...
                'SVD',SVD,'summary',1,...
                'nDirL',nDirL,'nDirR',nDirR};
            
            % run searchlight analysis
            lmva_spm(SL,SPM.xY.P,nii_out,ldafunction,'params',params,'isNP',1);
        end % s
        
        cd(home);
        varargout = {nii_out};
    case 'MVA_summary_raw'          % **new** searchlight analysis (only save summary value)
        sn  = varargin{1};
        glm = 10;
        blocksize= 5e07; % default value is 7e05
        
        vararginoptions({varargin{2:end}},...
            {'glm'});
        
        rsafunction = @bmw1_SL;
        Ndata = [4,3,7,7,7,49];
        resultNames = {'meanDist','uniLR','phiSummary','phiL','phiR','phiB'};
        
        home = cd;
        for s = sn
            glmDirSubj=fullfile(baseDir,glmName{glm},subj_name{s});
            cd(glmDirSubj);
            nii_out = {};
            
            % Define output files in 4D nifty file
            %varargout = {[dUL,dUR,dBL,dBR,dUex,dUin,dUun,...
            %  pL,pR,pBL,pBR,pLR,pUBL,pUBR,...
            %  vec(phiL)',vec(phiR)',vec(phiB)']}
            for f=1:numel(resultNames)
                for elem = 1:Ndata(f)
                    nii_out{end+1} = fullfile(glmDirSubj,...
                        sprintf('%s_%s.glm%d.nii,%d',subj_name{s},resultNames{f},glm, elem));
                end
            end
            % load search light definition, etc.
            load(fullfile(glmDirSubj,'SPM.mat'));
            T   = load(fullfile(glmDirSubj,'SPM_info.mat'));
            %SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
            SL  = load(fullfile(glmDirSubj,['vol_roi_160vox.mat']));
            
            params = {SPM, T.movType(T.regType==1)};
            %bmw1_SL(Y, SPM, condition, varargin)
            % run searchlight analysis
            %lmva_spm(SL,SPM.xY.P,nii_out,ldafunction,'params',params,'isNP',1);
            rsa.runSearchlight(SL,SPM.xY.P,nii_out,rsafunction,...
                'optionalParams',params,'idealBlock',blocksize);
            
        end % s
        
        cd(home);
        varargout = {nii_out};
        
    case 'surf_copy_searchlight'    % Copy search linght definition
        sn         = varargin{1};
        glmsource  = varargin{2};
        glmdest    = varargin{3};
        
        Slfiles = {'vol_roi_160vox.mat','vol_surf.mat','vol_roi_depth.nii'};
        source  = fullfile(baseDir,glmName{glmsource});
        dest    = fullfile(baseDir,glmName{glmdest});
        
        for s=sn;
            for f=1:3;
                SLsource    = fullfile(source,subj_name{s},Slfiles{f});
                SLdest      = fullfile(dest,subj_name{s},Slfiles{f});
                copyfile(SLsource,SLdest);
                fprintf('Copied file from: %s\n',SLsource);
                fprintf(' to: %s\n',SLdest);
            end;
        end;
    case 'surf_map_con'             % map contrast (.img) onto surface (.metric)
        % map volume images to metric file and save them in individual surface folder
        sn  = varargin{1};
        glm = varargin{2};
        
        hemisphere = [1:2];
        atlas = 1;
        fname = 'movement';
        vararginoptions({varargin{3:end}},{'atlas','hemisphere','fname'});
        
        Contrasts = {'all_movement','unimanual_left','unimanual_right','bimanual',...
            'UL>UR','B>U','BAsym>BSymEx','BAsym>BSymIn'};
        switch glm
            case {1,6}
                % Contrast no. (glm 1)
                %   1:80;      Individual movement type vs rest (t-contrast)
                %   81  ;      All movement vs rest (t-contrast)
                %   82  ;      All movement vs rest (F-contrast)
                %   83  ;       Unimanual left vs rest (t-contrast)
                %   84  ;       Unimanual right vs rest (t-contrast)
                %   85  ;       Bimanual vs rest (t-contrast)
                Conuse = [81 83 84 85];%
            case {3,4,7}
                % Contrast no. (glm 3,4)
                %   1:24;      Individual movement type vs rest (t-contrast)
                %   25  ;      All movement vs rest (t-contrast)
                %   26  ;      All movement vs rest (F-contrast)
                %   27  ;       Unimanual left vs rest (t-contrast)
                %   28  ;       Unimanual right vs rest (t-contrast)
                %   29  ;       Bimanual vs rest (t-contrast)
                Conuse = [25 27 28 29];
            case {2,8,9,10,11,12}
                % Contrast no. (glm 8)
                %   1:48;      Individual movement type vs rest (t-contrast)
                %   49  ;      All movement vs rest (t-contrast)
                %   50  ;      All movement vs rest (F-contrast)
                %   51  ;       Unimanual left vs rest (t-contrast)
                %   52  ;       Unimanual right vs rest (t-contrast)
                %   53  ;       Bimanual vs rest (t-contrast)
                %   54  ;       UniL>UniR (t-contrast)
                %   55  ;       Bi>Uni (t-contrast)
                %   56  ;       Bi-asymmetric>Bi-symmetric(ex) (t-contrast)
                %   57  ;       Bi-asymmetric>Bi-symmetric(in) (t-contrast)
                Conuse = [49 51 52 53 54 55 56 57];
            otherwise
                warning('glm should be smaller than 8');
                return;
        end
        
        
        for s = sn
            % Set .img file name
            glmDir = fullfile(baseDir,glmName{glm});
            images = [];
            column_name = [];
            for j = 1:numel(Contrasts)
                i = Conuse(j);
                if i==82
                    images{end+1} = fullfile(glmDir,subj_name{s},sprintf('spmF_%04d.nii',i));
                    column_name{end+1} = Contrasts{j};
                else
                    %images{end+1} = fullfile(glmDir,subj_name{s},sprintf('spmT_%04d.nii',i));
                    images{end+1} = fullfile(glmDir,subj_name{s},sprintf('con_%04d.nii',i));
                    column_name{end+1} = Contrasts{j};
                end
            end;
            images = images(~cellfun('isempty',images));
            
            % Caret stuffs
            for h=hemisphere
                caretSDir = fullfile(caretDir,[atlasA{atlas},subj_name{s}],hemName{h});
                specname = fullfile(caretSDir,[atlasA{atlas},subj_name{s} '.' hem{h}   '.spec']);
                
                % Load .coord file
                white = fullfile(caretSDir,[hem{h} '.WHITE.coord']);
                pial = fullfile(caretSDir,[hem{h} '.PIAL.coord']);
                C1 = caret_load(white);
                C2 = caret_load(pial);
                
                % Convert .image into .metric and save
                metric_out = fullfile(caretSDir,[subj_name{s}, sprintf('_Contrasts.glm%d.metric',glm)]);
                M = caret_vol2surf_own(C1.data, C2.data, images, 'ignore_zeros', 1);
                M.column_name = column_name;
                caret_save(metric_out, M);
                fprintf('Subj %d, Hem %d\n',s,h);
                
                % Smooth output .metric file (optional)
                % Load .topo file
                closed = fullfile(caretSDir,[hem{h} '.CLOSED.topo']);
                Out = caret_smooth(metric_out, 'coord', white, 'topo', closed);%,...
                %'algorithm','FWHM','fwhm',12);
                char(Out)
            end;
        end
    case 'surf_map_summary_distance'% Map summary distance (uniL,uniR,Bi,...)
        sn          = varargin{1};
        glm         = varargin{2};
        
        hemisphere  = [1:2];
        atlas       = 1;
        smooth      = 0;
        
        vararginoptions({varargin{3:end}},{'atlas','hemisphere','smooth','fnames','colnames'});
        
        fnames       = {{sprintf('meanDist.glm%d',glm),...
            sprintf('uniLRext.glm%d',glm),...
            sprintf('uniLRint.glm%d',glm)}};
        
        
        colnames     = {{'Unimanual left','Unimanual right',...
            'Bimanual left','Bimanual right',...
            'Uextrinsic','Uintrinsic'}};
        
        % convert .nii data into .metric data
        for f=1:numel(fnames)
            bmw1_imana('surf_map_any',sn,fnames{f},'glm',glm,'column_names',colnames{f},...
                'atlas',atlas,'hemisphere',hemisphere,'smooth',smooth);
        end
    case 'surf_map_summary_covhat'  % Map summary covariance estimate across hands or uni/bi
        sn          = varargin{1};
        glm         = varargin{2};
        
        hemisphere  = [1:2];
        atlas       = 1;
        smooth      = 0;
        
        vararginoptions({varargin{3:end}},{'atlas','hemisphere','smooth','fnames','colnames'});
        
        fnames       = {{sprintf('covLR.glm%d',glm),...
            sprintf('covULBL.glm%d',glm),...
            sprintf('covURBR.glm%d',glm)}};
        
        
        colnames     = {{'corr(ULUR)','corr(ULBL)','corr(URBR)'}};
        
        % convert .nii data into .metric data
        for f=1:numel(fnames)
            bmw1_imana('surf_map_any',sn,fnames{f},'glm',glm,'column_names',colnames{f},...
                'atlas',atlas,'hemisphere',hemisphere,'smooth',smooth);
        end
    case 'surf_groupmap_con'        % Make group statistical map
        sn      = varargin{1};
        glm     = varargin{2};
        conuse  = [1:8];
        
        Contrasts = {'all_movement','unimanual_left','unimanual_right','bimanual',...
            'UL_lt_UR','B_lt_U','BAsym_lt_BSymEx','BAsym_lt_BSymIn'};
        
        INname = repmat({sprintf('Contrasts.glm%d',glm)},1,length(conuse));
        
        OUTname = Contrasts(conuse);
        
        bmw1_imana('surf_group_make',sn,'INname',INname,'OUTname',OUTname,...
            'inputcol',[conuse],'replaceNaN',[conuse]);
        
        bmw1_imana('surf_group_smooth',sn,'SPMname',OUTname,'Smooth_iterations',...
            repmat(15,1,length(conuse)));
        
        bmw1_imana('surf_group_cSPM',sn,'SPMname',OUTname,'sqrtTransform',...
            zeros(1,length(conuse)),'SummaryName','.ConSummary.metric');
    case 'surf_groupmap_meandist'   % Make group statistical map
        sn      = varargin{1};
        glm     = varargin{2};
        conuse  = [1:6];
        
        Contrasts = {'UniLeft','UniRight',...
            'BiLeft','BiRight',...
            'Uex','Uin'};
        
        INname = repmat({sprintf('meanDist.glm%d',glm)},1,length(conuse));
        
        OUTname = Contrasts(conuse);
        
        bmw1_imana('surf_group_make',sn,'INname',INname,'OUTname',OUTname,...
            'inputcol',[conuse],'replaceNaN',[conuse]);
        
        bmw1_imana('surf_group_smooth',sn,'SPMname',OUTname,'Smooth_iterations',...
            repmat(15,1,length(conuse)));
        
        bmw1_imana('surf_group_cSPM',sn,'SPMname',OUTname,'sqrtTransform',...
            ones(1,length(conuse)),'SummaryName','.distSummary.metric');
    case 'surf_groupmap_covhat'     % Make group statistical map
        sn      = varargin{1};
        glm     = varargin{2};
        conuse  = [1:3];
        
        Contrasts = {'corrULUR','corrULBL','corrURBR'};
        
        INname = repmat({sprintf('covLR.glm%d',glm)},numel(Contrasts),1);
        
        OUTname = Contrasts(conuse);
        
        bmw1_imana('surf_group_make',sn,'INname',INname,'OUTname',OUTname,...
            'inputcol',[conuse],'replaceNaN',[conuse]);
        
        bmw1_imana('surf_group_smooth',sn,'SPMname',OUTname,'Smooth_iterations',...
            repmat(15,1,length(conuse)));
        
        bmw1_imana('surf_group_cSPM',sn,'SPMname',OUTname,'sqrtTransform',...
            zeros(1,length(conuse)),'SummaryName','.corrSummary.metric');
        
    case 'surf_map_any'             % map searchlight result (.nii) onto the surface (.metric)
        % map volume images to metric file and save them in individual surface folder
        sn          = varargin{1};
        fname       = varargin{2};
        hemisphere  = [1:2];
        atlas       = 1;
        glm         = 1;
        column_names= {};
        smooth      = 0;
        colname     = {};
        vararginoptions({varargin{3:end}},{'atlas','hemisphere','glm','column_names','smooth'});
        
        for s = sn
            for h = hemisphere
                % Load .coord files
                caretSDir = fullfile(caretDir,[atlasA{atlas},subj_name{s}],hemName{h});
                specname = fullfile(caretSDir,[atlasA{atlas},subj_name{s} '.' hem{h}   '.spec']);
                white = fullfile(caretSDir,[hem{h} '.WHITE.coord']);
                pial = fullfile(caretSDir,[hem{h} '.PIAL.coord']);
                
                C1 = caret_load(white);
                C2 = caret_load(pial);
                
                % Setting input and output files
                images  = [];
                if ~iscell(fname)
                    
                    image   = fullfile(baseDir, glmName{glm},subj_name{s},[subj_name{s},'_',fname, '.nii']);
                    % Check image size
                    A = spm_vol(image);
                    imsize = length(A);
                    if imsize==1
                        images{1} = image;
                    else
                        for i = 1:imsize
                            images{end+1} = fullfile(baseDir, glmName{glm},subj_name{s},...
                                [subj_name{s},'_',fname, '.nii,', num2str(i)]);
                        end
                    end
                    metric_out  = fullfile(caretSDir,[subj_name{s},'_',fname, '.metric']);
                else
                    for f = 1:length(fname)
                        image = fullfile(baseDir, glmName{glm},subj_name{s},...
                            [subj_name{s},'_',fname{f}, '.nii']);
                        
                        A = spm_vol(image);
                        imsize = length(A);
                        if imsize==1
                            images{end+1} = image;
                        else
                            for i = 1:imsize
                                images{end+1} = fullfile(baseDir, glmName{glm},subj_name{s},...
                                    [subj_name{s},'_',fname{f}, '.nii,', num2str(i)]);
                            end
                        end
                        
                        if isempty(column_names);
                            colname{f} = fname{f};
                        else
                            colname = column_names;
                        end
                        
                    end
                    metric_out  = fullfile(caretSDir,[subj_name{s},'_',fname{1}, '.metric']); % use name of first image for metric output
                    
                end
                
                % Convert volume (.nii) data to surface (.metric) data
                M = caret_vol2surf_own(C1.data,C2.data,images,'ignore_zeros',1,'column_names',colname);
                caret_save(metric_out,M);
                
                if smooth
                    % Smooth output .metric file (optional)
                    % Load .topo file
                    closed = fullfile(caretSDir,[hem{h} '.CLOSED.topo']);
                    Out = caret_smooth(metric_out, 'coord', white, 'topo', closed);
                    char(Out)
                end
                fprintf('Subj %d, Hem %d\n',s,h);
            end
        end
    case 'surf_group_make'              % Make the group caret surface files
        sn = varargin{1};
        
        INname = repmat({'distanceSummary'},1,8);
        
        OUTname = {'ActiveSingleDigit','ActiveMultidigit',...
            'PassiveSingleDigit','PassiveMultiDigit',...
            'ActiveSingle-PassiveSingle','ActiveSingle-PassiveMulti',...
            'ActiveMulti-PassiveSingle','ActiveMulti-PassiveMulti'};
        
        inputcol    = [2:9];
        replaceNaN  = ones(1,8);
        
        vararginoptions(varargin(2:end),{'atlas','INname','OUTname','inputcol','replaceNaN'});
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
            
            %----loop over each input metric file and make a group metric file
            for j = 1:length(INname);
                
                %----define names of subj metric files
                for i = sn;
                    infilenames{j}{i}=[caretDir filesep 'x' subj_name{i} filesep hemName{h} filesep subj_name{i} '_' INname{j} '.metric'];
                end;
                
                %----name for group metric file in average surface folder
                outfilenames{j}=[surfaceGroupDir filesep hem{h} '.' OUTname{j} '.metric'];
                %----make the group metric
                fprintf('hem: %i  image: %i \n', h,j);
                caret_metricpermute(infilenames{j},'outfilenames',outfilenames(j),'inputcol',inputcol(j),'replaceNaNs',replaceNaN(j));
            end;
        end;
    case 'surf_group_smooth'            % Smooth the group maps
        sn = varargin{1};
        
        SPMname = {'ActiveSingleDigit','ActiveMultidigit',...
            'PassiveSingleDigit','PassiveMultiDigit',...
            'ActiveSingle-PassiveSingle','ActiveSingle-PassiveMulti',...
            'ActiveMulti-PassiveSingle','ActiveMulti-PassiveMulti'};
        Smooth_iterations = repmat(15,1,numel(SPMname));
        
        vararginoptions(varargin(2:end),{'SPMname','files','Smooth_iterations'});
        
        files=1:length(SPMname);
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym' filesep hemName{h}];
            cd(surfaceGroupDir)
            %----define name of coord and topology
            coordfile=[caretDir filesep 'fsaverage_sym'   filesep hemName{h} filesep hem{h} '.WHITE.coord'];
            topofile=[caretDir filesep 'fsaverage_sym'   filesep hemName{h} filesep hem{h} '.CLOSED.topo'];
            %----get the full directory name of the metric files and the smoothed metric files that we create below
            for i=files
                filename = [surfaceGroupDir filesep hem{h} '.' SPMname{i} '.metric']; % unsmoothed
                sfilename = caret_smooth(filename,'coord',coordfile,'topo',topofile,'iterations', Smooth_iterations(i));
            end;
            %----smooth the metric files and save them with the prefix 's'
        end;
    case 'surf_group_cSPM'
        sn = varargin{1};
        
        SPMname = {'ActiveSingleDigit','ActiveMultidigit',...
            'PassiveSingleDigit','PassiveMultiDigit',...
            'ActiveSingle-PassiveSingle','ActiveSingle-PassiveMulti',...
            'ActiveMulti-PassiveSingle','ActiveMulti-PassiveMulti'};
        
        sqrtTransform = [1 1 1 1 1 1 1 1]; % Should you take ssqrt before submitting?
        SummaryName = '.summary.metric';
        
        vararginoptions(varargin(2:end),{'SPMname','sqrtTransform','SummaryName'});
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h}];% filesep smoothdir];
            %----get the full directory name of the metric files and the smoothed metric files that we create below
            for i=1:length(SPMname);
                sfilenames{i}=[surfaceGroupDir filesep 's' hem{h} '.' SPMname{i} '.metric']; % smoothed
            end;
            %----loop over the metric files and calculate the cSPM of each with the smoothed metrics
            for i=1:length(SPMname);
                Data = caret_load(sfilenames{i});
                
                if sqrtTransform(i)
                    Data.data=ssqrt(Data.data);
                end;
                
                cSPM                            = caret_getcSPM('onesample_t','data',Data.data(:,sn));
                data(:,i)                       = cSPM.con(1).con; % mean
                data(:,i+length(SPMname))       = cSPM.con(1).Z; % T
                column_name{i}                  = ['mean_' SPMname{i}];
                column_name{i+length(SPMname)}  = ['T_' SPMname{i}];
                
                caret_savecSPM([surfaceGroupDir    filesep hem{h} '.' SPMname{i} '_stats.metric'],cSPM); %filesep smoothdir
                save([surfaceGroupDir  filesep   'cSPM_' SPMname{i} '.mat'],'cSPM');%filesep smoothdir
            end;
            C = caret_struct('metric','data',data,'column_name',column_name);
            caret_save([surfaceGroupDir  filesep hem{h} SummaryName],C); %filesep smoothdir
        end;
    case 'surf_nii_calc'            % Do calculation using .nii and save as .nii
        sn = varargin{1};
        glm= 10;
        mode = 'LR';
        vararginoptions(varargin(2:end),{'mode'});
        
        for s = sn % loop over subjects
            disp(subj_name{s});
            glmDirSubj = fullfile(baseDir, glmName{glm}, subj_name{s});
            
            switch mode
                case 'LR'
                    % calculate d(unrelated)-d(extrinsic) or d(unrelated)-d(intrinsic)
                    %--- Load .nii data for distance and sigma
                    for i=1:3
                        Vin{i} = fullfile(glmDirSubj,sprintf('%s_uniLR.glm%d.nii,%d',subj_name{s},glm,i));
                    end
                    Vout{1} = fullfile(glmDirSubj,sprintf('%s_uniLRext.glm%d.nii',subj_name{s},glm));
                    Vout{2} = fullfile(glmDirSubj,sprintf('%s_uniLRint.glm%d.nii',subj_name{s},glm));
                    
                    spm_imcalc(Vin,Vout{1},'i3-i1');
                    spm_imcalc(Vin,Vout{2},'i3-i2');
                case 'UB'
                    % calculate cov(u,b)/sqrt(var(u)var(b))
                    %--- Load .nii data for distance and sigma
                    for i=1:7
                        Vin{i} = fullfile(glmDirSubj,sprintf('%s_phiSummary.glm%d.nii,%d',subj_name{s},glm,i));
                    end
                    Vout{1} = fullfile(glmDirSubj,sprintf('%s_covLR.glm%d.nii',subj_name{s},glm));
                    Vout{2} = fullfile(glmDirSubj,sprintf('%s_covULBL.glm%d.nii',subj_name{s},glm));
                    Vout{3} = fullfile(glmDirSubj,sprintf('%s_covURBR.glm%d.nii',subj_name{s},glm));
                    
                    spm_imcalc(Vin,Vout{1},'(i5./(ssqrt(i1).*ssqrt(i2))).*(i1>1e-3&i2>1e-3)');
                    spm_imcalc(Vin,Vout{2},'(i6./(ssqrt(i1).*ssqrt(i3))).*(i1>1e-3&i2>1e-3)');
                    spm_imcalc(Vin,Vout{3},'(i7./(ssqrt(i2).*ssqrt(i4))).*(i1>1e-3&i2>1e-3)');
            end
        end
    case 'surf_metric_calc'
        
    case 'ROI_define_all'                                                   % Define surface and BG ROIs (CB ROIs to be implemented)
        sn = varargin{1};
        glm = varargin{2};
        
        % Define surface ROIs
        bmw1_uwo_imana('ROI_define_surf',sn,'glm',glm); % only for active voxels
        bmw1_uwo_imana('ROI_define_surf',sn,'goalNodes',-1,'glm',glm); % all voxels within ROI
        
        % Define BG ROIs
        for step = 1:3
            fprintf('ROI_BG_MNImask: step%d',step);
            if step==1; sn = 1:numel(subj_name); end;
            bmw1_uwo_imana('ROI_BG_MNImask',step,sn);
        end
        bmw1_uwo_imana('ROI_define_BG',sn);
        
        % Define CB ROIs
        % (not implemented yet)
    case 'ROI_makepaint_propatlas'            % ROI: Generates ROI paint file
        for h=1:2
            groupDir=[caretDir filesep 'fsaverage_sym' filesep hemName{h} ];
            cd(groupDir);
            M=caret_load([hem{h} '.propatlas.metric']);
            P=caret_load(['lateral.paint']);
            C=caret_load([hem{h} '.FLAT.coord']);
            
            % need to adjust?
            S1=sum(M.data(:,[1 2 3 4]),2);
            M1=sum(M.data(:,[7 8]),2);
            PM=M.data(:,9);
            [Prop,ROI]=max([S1 M1 PM],[],2);
            ROI(Prop<0.2)=0;
            ROI(ROI==1 & (P.data(:,2)==0 | C.data(:,2)<-9 | C.data(:,2)>20))=0; % Hand area M1
            ROI(ROI==2 & (P.data(:,2)==0 | C.data(:,2)<-9 | C.data(:,2)>20))=0; % Hand area S1
            ROI(ROI==3 & P.data(:,2)==0)=5; % SMA / PreSMA
            ROI(ROI==3 & C.data(:,2)<-9)=4; % ventral premotor
            ROI(ROI==3 & (C.data(:,2)>23 | C.data(:,1)<-40))=0; % PMd
            
            ROI(P.data(:,1)==1)=6;
            ROI(sum(M.data(:,11:12),2)>0.2)=7;
            
            % Save Paint file
            Paint=caret_struct('paint','data',ROI,'paintnames',regname,'column_name',{'ROI'});
            caret_save(['ROI.paint'],Paint);
        end;
    case 'ROI_define_surf'
        % the contrasts are mapped onto the surface map using surf_map_con
        % and here, the goalNodes number of most activated voxels are
        % picked for each of the anatomically defines regions.
        % contrast are saved in the metric file specified by the id 'fname'
        %
        % get region data for movement vs rest contrast using the first glm
        % bmw1_imana('ROI_define',1,'fname','movement','glm',1)
        %
        % get region data for movement vs rest contrast using the second glm
        % bmw1_imana('ROI_define',1,'fname','movement','glm',2)
        
        sn          = varargin{1};
        goalNodes   = -1; %800;
        atlas       = 1;
        glm         = 1;
        fname       = 'movement';
        vararginoptions({varargin{2:end}},{'goalNodes','atlas','glm','fname'});
        
        % check fsaverage_sym
        
        for s=sn            
            if exist(fullfile(regDir,subj_name{s}))
                cd(fullfile(regDir,subj_name{s}));
            else
                mkdir(fullfile(regDir,subj_name{s}));
                cd(fullfile(regDir,subj_name{s}));
            end
            
            glmDir = fullfile(baseDir,glmName{glm});
            for h=1:2
                caretSDir = fullfile(caretDir,[atlasA{atlas},subj_name{s}],hemName{h});
                glmFileName = fullfile(caretSDir,[subj_name{s}, sprintf('_Contrasts_%d.metric',glm)]);
                
                % borrowed from sl1 ROI.paint
                caretFSAVGDir = fullfile(caretDir,atlasname{atlas});
                if ~exist(caretFSAVGDir, 'dir');
                    home = pwd;
                    cd(baseDir); cd ../../;
                    source = fullfile(cd, 'Atlas_templates', 'fsaverage_sym');
                    dest = fullfile(caretDir, 'fsaverage_sym');
                    copyfile(source, dest);
                    cd(home);
                end
                C = caret_load(fullfile(caretDir,atlasname{atlas},hemName{h},['ROI.paint']));
                
                if goalNodes>0
                    %M = caret_load(fullfile(caretDir,[atlasA{atlas} subj_name{sn}],hemName{h},[subj_name{sn} sprintf('_func_%s_%d.metric',fname,glm)]));
                    M = caret_load(glmFileName);
                    tavrg = (M.data(:,1)); % t-contrast (movement vs rest)
                end
                
                caretSubjDir = fullfile(caretDir,['x' subj_name{s}]);
                %file = fullfile(glmDir,subj_name{s},'mask.img');
                file = fullfile(glmDir,subj_name{s},'mask.nii');
                
                for i = 1:numregions_surf
                    indx        = (C.data(:,1)==i);
                    numNodes    = sum(indx);
                    
                    R{i+(h-1)*numregions_surf}.type = 'surf_nodes';
                    R{i+(h-1)*numregions_surf}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                    R{i+(h-1)*numregions_surf}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                    R{i+(h-1)*numregions_surf}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                    R{i+(h-1)*numregions_surf}.linedef=[5,0,1];
                    R{i+(h-1)*numregions_surf}.image=file;
                    R{i+(h-1)*numregions_surf}.name=[subj_name{s} '_' regname{i} '_' hem{h}];
                    
                    if goalNodes<0
                        R{i+(h-1)*numregions_surf}.location=find(indx);
                        % R{i+(h-1)*numregions}.location=find(indx);
                    else
                        R{i+(h-1)*numregions_surf}.location=find(indx & tavrg>prctile(tavrg(indx),100*(1-goalNodes/numNodes)));
                        % R{i+(h-1)*numregions}.location=find(indx & tavrg(:,sn)>prctile(tavrg(indx,sn),100*(1-goalNodes/numNodes)));
                    end
                end;
            end;
            
            R = region_calcregions(R);
            if goalNodes<0
                save(sprintf('%s_regAll_%s.mat',subj_name{s},fname),'R'); % Because purely anatomical ROI is not affected by activation
            else
                save(sprintf('%s_reg%d_%s_%d.mat',subj_name{s},goalNodes,fname,glm),'R');
            end
            fprintf('%d\n',s);
        end
        
        varargout = {};
    case 'ROI_BG_MNImask'                     % Do segmentation for BG using FSL in MNI space
        %sl1_imana('ROI_BG_MNImask', 1)
        step = varargin{1};
        sn = varargin{2};
        
        switch (step)
            case 1 % run FSL routine
                for s= sn%1:numel(subj_name)
                    IN= fullfile(anatomicalDir, subj_name{s}, [subj_name{s}, '_anatomical.nii']);
                    outDir = fullfile(baseDir, 'basal_ganglia', 'FSL');
                    if ~exist(outDir,'dir')
                        mkdir(outDir);
                    end
                    OUT= fullfile(outDir, subj_name{s}, [subj_name{s}, '_BG.nii']);
                    %calc with FSL
                    comm=sprintf('run_first_all -i %s -o %s', IN, OUT);
                    fprintf('%s\n',comm);
                    system(comm);
                end
            case 2 % make the ROI images in subject space and do mni transform
                %         10 Left-Thalamus-Proper 40
                %         11 Left-Caudate 30
                %         12 Left-Putamen 40
                %         13 Left-Pallidum 40
                %         49 Right-Thalamus-Proper 40
                %         50 Right-Caudate 30
                %         51 Right-Putamen 40
                %         52 Right-Pallidum 40
                BGnumber= [11 13 12 10; 50 52 51 49];
                %'CaudateN' 'Pallidum', 'Putamen' 'Thalamus'
                for s= sn%1:numel(subj_name)
                    %----deform info for basal ganglia ROI to individual space
                    nam_def = fullfile(anatomicalDir,subj_name{s}, [subj_name{s},'_anatomical_seg_sn.mat']);
                    mniDir = fullfile(baseDir, 'basal_ganglia', 'FSL', 'MNI', subj_name{s});
                    if ~exist(mniDir,'dir')
                        mkdir(mniDir);
                    end
                    
                    for h=1:2
                        for i=1:numregions_BG
                            fprintf('Working on subj: %i region: %s \n', s, [regname{i+numregions_surf},'_',hem{h}])
                            
                            %----get names!
                            IN= fullfile(baseDir, 'basal_ganglia', 'FSL', subj_name{s},...
                                [subj_name{s},'_BG_all_fast_firstseg.nii']);
                            
                            OUT{i}= fullfile(baseDir, 'basal_ganglia', 'FSL', subj_name{s},...
                                [subj_name{s},'_',regname{i+numregions_surf},'_',hem{h},'.nii']);
                            
                            OUT_MNI{i}= fullfile(baseDir, 'basal_ganglia', 'FSL', 'MNI', subj_name{s},...
                                [subj_name{s},'_',regname{i+numregions_surf},'_',hem{h}, '.nii']);
                            
                            %----make subj specific ROI image; still in MNI space!
                            spm_imcalc_ui(IN,OUT{i},sprintf('i1==%d',BGnumber(h,i)));
                        end
                        %----do deformation
                        spmj_normalization_write(nam_def, OUT,'outimages',OUT_MNI);
                    end
                end
            case 3 %make the avrg mask image
                for h=1:2
                    for i=1:numregions_BG
                        for s = 1:numel(sn)%1:numel(subj_name)
                            IN{s} = fullfile(baseDir, 'basal_ganglia', 'FSL', 'MNI', subj_name{sn(s)},...
                                [subj_name{sn(s)},'_',regname{i+numregions_surf},'_',hem{h}, '.nii']);
                        end
                        outDir = fullfile(baseDir, 'basal_ganglia', 'FSL', 'MNI', 'avrg');
                        if ~exist(outDir, 'dir');
                            mkdir(outDir);
                        end
                        OUT = fullfile(outDir,...
                            ['avrg_',regname{i+numregions_surf},'_',hem{h}, '.nii']);
                        spmj_imcalc_mtx(IN,OUT,'mean(X)');
                    end
                end
        end
    case 'ROI_define_BG'                      % Direct specification of the BG ROIS: Joern
        %sl1_imana('ROI_define_BG', 1:16)
        sn      = varargin{1};
        %glm     = varargin{2};
        regtype = sprintf('regAll_movement.mat');
        
        vararginoptions(varargin(2:end), {'regtype'});
        
        % -------------------------------
        % Read in BG-regions in MNI space: define possibly different
        for h=1:2
            for i=1:numregions_BG
                %----make subj specific ROI image; still in MNI space!
                %IN=fullfile(baseDir, 'basal_ganglia', [regname{i+numregions_surf},'_',hem{h}, '.nii']);
                %BGreg{h,i}=region('roi_image',IN,255);
                % defined through FSL BG
                IN=fullfile(baseDir, 'basal_ganglia', 'FSL', 'MNI', 'avrg',...
                    ['avrg_',regname{i+numregions_surf},'_',hem{h}, '.nii']);
                BGreg{h,i}=region('image',IN,0.001);
                BGreg{h,i}=region_calcregions(BGreg{h,i});
            end;
        end;
        
        % -------------------------------
        % Transform each of the regions into the individual space
        cd(regDir);
        for s=sn
            nam_def = fullfile(anatomicalDir,subj_name{s}, [subj_name{s},'_anatomical_seg_sn.mat']);
            for c=1:2
                load(fullfile('./',subj_name{s},[subj_name{s},'_',regtype]));
                
                for h=1:2
                    for i = 1:numregions_BG
                        num = numregions_surf*2 + numregions_BG*(h-1) + i; % Append BG ROI info
                        R{num} = region_deformation(BGreg{h,i},nam_def);
                        R{num}.xyz_mni = BGreg{h,i}.data;             % Use the original MNI space coordinates
                        R{num}.name = [subj_name{s},'_',regname_BG{i},'_',hem{h}];
                    end;
                end; %hem
                save(fullfile('./',subj_name{s},[subj_name{s},'_',regtype]),'R');
            end %con
            fprintf('Working on subj: %d\n',s);
        end; %subj
        varargout={R};
        
    case 'ROI_timeseries_get'                       % Extract raw and estimated time series from ROIs
        sn = varargin{1};
        glm = varargin{2};
        ROI = 'all';
        pre=ceil(10*2.7);
        post=ceil(10*2.7);
        vararginoptions(varargin(3:end),{'ROI','pre','post'});
        
        glmDir = fullfile(baseDir,glmName{glm});
        T=[];
        for s=sn;
            subj = subj_name{s};
            fprintf('%s\n',subj);
            
            % load SPM.mat
            cd(fullfile(glmDir,subj));
            load SPM;
            
            % load ROI definition
            switch ROI
                case 'all'
                    load(fullfile(regDir,subj,[subj '_regAll_movement.mat']));
                otherwise
                    load(fullfile(regDir,subj,sprintf('%s_reg%s_movement_%d',subj,ROI,glm)));
            end
            
            % choose ROI?
            %R = {R{[2 9]}};
            
            % extract time series data
            [y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM,R);
            
            D = spmj_get_ons_struct(SPM);
            
            for r=1:size(y_raw,2)
                for i=1:size(D.block,1);
                    D.y_adj(i,:)=cut(y_adj(:,r),pre,round(D.ons(i))-1,post,'padding','nan')';
                    D.y_hat(i,:)=cut(y_hat(:,r),pre,round(D.ons(i))-1,post,'padding','nan')';
                    D.y_res(i,:)=cut(y_res(:,r),pre,round(D.ons(i))-1,post,'padding','nan')';
                end;
                Roinames = strsplit(R{r}.name,'_');
                Roi = Roinames{2};
                Side = Roinames{3};
                
                roidx = cellfun(@strcmp,regname,repmat({Roi},size(regname)));
                sidx = cellfun(@strcmp,hem,{Side,Side});
                
                len = size(D.event,1);
                
                D.SN        = ones(len,1)*s;
                D.region    = ones(len,1)*r;
                D.name      = repmat({R{r}.name},len,1);
                D.type      = D.event;
                D.regSide   = repmat(find(sidx),len,1);
                D.regType   = repmat(find(roidx),len,1);
                D.glmType   = repmat(glm,len,1);
                T           = addstruct(T,D);
            end;
        end;
        save(fullfile(regDir,sprintf('hrf_%s_%d.mat',ROI,glm)),'-struct','T');
        varargout={T};
    case 'ROI_timeseries_plot1'                     % Plot extracted time series
        s = varargin{1};
        glm = varargin{2};
        roi = varargin{3};
        ROI = 'all';
        
        T = load(fullfile(regDir,sprintf('hrf_%s_%d.mat',ROI,glm)));
        T = getrow(T,ismember(T.regType,roi));
        pre = ceil(10*2.7);
        post = ceil(10*2.7);
        switch glm
            case {1,2,3,4}
                types = {[1:6],[7:12],[13:48]};
        end
        
        % loop over trial type
        figure('name',sprintf('%s',glmName{glm}));
        typename = {'Unimanual left trials','Unimanual right tirlas','Bimanual trials'};
        for type = 1:3
            % Unimanual left trial
            for h = 1:2
                subset_type = ismember(T.type,types{type});
                subset      = T.regSide==h & ismember(T.SN,s) & subset_type;
                
                subplot(3,2,h+2*(type-1));
                traceplot([-pre:post],T.y_adj,'errorfcn','stderr',...
                    'split',[T.regType],'subset',subset,...
                    'leg',regname(roi),'leglocation','bestoutside'); % ,
                hold on;
                traceplot([-pre:post],T.y_hat,'linestyle','--',...
                    'split',[T.regType],'subset',subset,...
                    'linewidth',3); % ,
                drawline([-7 7 14]/(TR/1000),'dir','vert','linestyle','--');
                drawline([0],'dir','horz','linestyle','--');
                hold off;
                xlabel('TR');
                ylabel('activation');
                title(typename{type})
                drawline(0);
            end;
        end
    case 'ROI_timeseries_plot2'                     % Plot extracted time series
        s = varargin{1};
        glm = varargin{2};
        ROI = varargin{3};
        roi = varargin{4};
        
        lineColor = jet(numel(s));
        shadeColor= min(1,lineColor + 0.3);
        lineColor = mat2cell(lineColor,ones(numel(s),1),3);
        shadeColor= mat2cell(shadeColor,ones(numel(s),1),3);
        
        % load data
        T = load(fullfile(regDir,sprintf('hrf_%s_%d.mat',ROI,glm)));
        T = getrow(T,ismember(T.regType,roi));
        
        % elevate baseline for each subject for visibility
        baseline= 1.5*T.SN.*double(T.regType<9)+0.5*T.SN.*double(T.regType>8);
        T.y_adj = bsxfun(@minus,T.y_adj,baseline);
        T.y_hat = bsxfun(@minus,T.y_hat,baseline);
        
        pre = 10;
        post = 10;
        
        % trial types (uni left/uni right/bi)
        switch glm
            case 1
                types = {[1:8],[9:16],[16:80]};
            case 3
                types = {[1:4],[5:8],[9:24]};
            case 4
                types = {[1:4],[5:8],[9:24]};
            case 6
                types = {[1:8],[9:16],[16:80]};
            case 7
                types = {[1:4],[5:8],[9:24]};
            case {8,9,10}
                types = {[1:6],[7:12],[13:48]};
        end
        
        % loop over trial type and roi
        typename = {'Unimanual left trials','Unimanual right tirlas','Bimanual trials'};
        for reg=roi;
            figure('name',sprintf('%s-%s',glmName{glm},regname{reg}));
            for type = 1:3
                % Unimanual left trial
                for h = 1:2
                    subset_type = ismember(T.type,types{type});
                    subset      = T.regType==reg & T.regSide==h & ismember(T.SN,s) & subset_type;
                    
                    subplot(3,2,h+2*(type-1));
                    traceplot([-pre:post],T.y_adj,'errorfcn','stderr',...
                        'split',[T.SN],'subset',subset,...
                        'leg',subj_name(s),'leglocation','bestoutside',...
                        'linecolor',lineColor,'patchcolor',shadeColor); % ,
                    hold on;
                    traceplot([-pre:post],T.y_hat,'linestyle','--',...
                        'split',[T.SN],'subset',subset,...
                        'linewidth',3,'linecolor',lineColor); % ,
                    drawline([-7 7 14]/2.72,'dir','vert','linestyle','--');
                    %drawline([0],'dir','horz','linestyle','--');
                    hold off;
                    xlabel('TR');
                    ylabel('activation (a.u.)');
                    set(gca,'yticklabel',{});
                    title(typename{type})
                    drawline(0);
                end;
            end;
        end;
    case 'ROI_timeseries_fit'                       % Fit the hrf to the data extracted from ROIs
        sn  = varargin{1};
        glm = 10;
        regions = [2];
        hem = 1:2;
        fithrf = [1,2]'; % hrf parameter to be fitted
        LB     = [1  6   0.2 0.2  3  -3 0.1]';
        UB     = [11 25  4   4    12 4  1.5]';
        duration = 1;%1200/2720;%7500/2720;
        onsetshift = -5000/2720;%-2500/2720;
        
        vararginoptions({varargin{2:end}},{'glm','region','hem','fit','LB','UB'});
        
        P0      = [3.389  16  1   1    6  1.0140    0.4990]';
        %P0      = [6  16  1   1    6  2.0  1.0]';
        pre     = 4;
        post    = 12;
        
        T = []; Ts = [];
        for s=sn
            subj = subj_name{s};
            fprintf('%s\n',subj);
            
            cd(fullfile(baseDir,glmName{glm},subj_name{s}));
            
            load SPM;
            load(fullfile(regDir,subj,[subj '_regAll_movement.mat']));
            
            % Clearn default onset and duration
            if ismember(fithrf,[6,7])
                for r = 1:length(SPM.nscan)
                    for u=1:length(SPM.Sess(r).U)
                        SPM.Sess(r).U(u).dur = ones(size(SPM.Sess(r).U(u).dur))*duration; % 1
                        SPM.Sess(r).U(u).ons = SPM.Sess(r).U(u).ons+onsetshift; % return to TR at announce trial
                    end;
                    SPM.Sess(r).U=spm_get_ons(SPM,r);
                end;
            end
            
            for h = hem
                for reg = regions
                    if reg<9
                        roi = numregions_surf*(h-1)+reg;
                    elseif reg>=9
                        roi = 2*numregions_surf + numregions_BG*(h-1)+(reg-8);
                    end
                    fprintf('%s\n',R{roi}.name);
                    
                    % Get data
                    Y = region_getdata(SPM.xY.VY,R{roi});
                    
                    % Check Error before
                    Ypre        = spm_filter(SPM.xX.K,SPM.xX.W*Y);
                    Yres        = spm_sp('r',SPM.xX.xKXs,Ypre);
                    Epre        = sum(sum(Yres.^2))/numel(Yres(:));
                    
                    % Fit a common hrf
                    e = 1;
                    P0_ = P0;
                    iter = 1;
                    while ((e >= 1)&&(iter<10))%((e >= 1)&&(iter<5))
                        fprintf('%d th iteration...',iter);
                        
                        [P,SPM,Yhat,Yres] = spmj_fit_hrf(SPM,Y,...
                            'fit',fithrf,'LB',LB,'UB',UB,'P0',P0_);
                        
                        % update initial value
                        P0_(fithrf) = P0(fithrf)+0.1*rand*(UB(fithrf)-LB(fithrf));
                        iter = iter+1;
                        
                        % Check Error after
                        Epost       = sum(sum(Yres.^2))/numel(Yres(:));
                        e           = Epost/Epre;
                    end
                    
                    fprintf('Epost/Epre: %d\n',e);
                    
                    % Parameter values
                    T_.fit = fithrf';
                    T_.P = P(fithrf)';
                    T_.R2 = 1-(trace(Yres*Yres')/(trace(Yres*Yres')+trace(Yhat*Yhat')));
                    T_.Eratio = Epost/Epre;
                    T_.regSide = h;
                    T_.regType = reg;
                    T_.SN = s;
                    T_.hemname = {hemName{h}};
                    T_.regname = {regname{reg}};
                    T = addstruct(T,T_);
                    
                    % Timeseries
                    D = spmj_get_ons_struct(SPM);
                    
                    y_hat = mean(Yhat,2);
                    y_res = mean(Yres,2);
                    for i=1:size(D.block,1);
                        D.y_hat(i,:)=cut(y_hat,pre,round(D.ons(i)),post,'padding','nan')';
                        D.y_res(i,:)=cut(y_res,pre,round(D.ons(i)),post,'padding','nan')';
                        D.y_adj(i,:)=D.y_hat(i,:)+D.y_res(i,:);
                    end;
                    D.regType   = ones(size(D.event,1),1)*reg;
                    D.regSide   = ones(size(D.event,1),1)*h;
                    D.sn        = ones(size(D.event,1),1)*s;
                    Ts          = addstruct(Ts,D);
                end
            end
        end
        
        % plot
        figure('name',[what,'_parameters']);
        subplot(2,2,1)
        lineplot(T.regSide,T.P(:,1),'split',T.SN); title('onset')
        subplot(2,2,2)
        try;lineplot(T.regSide,T.P(:,2),'split',T.SN); title('duration');catch;end;
        subplot(2,2,3)
        lineplot(T.regSide,T.R2,'split',T.SN); title('R^2')
        subplot(2,2,4)
        lineplot(T.regSide,T.Eratio,'split',T.SN); title('Reduction of error by fitting')
        
        figure('name',[what,'_timeseries'])
        for h=hem
            subplot(1,2,h)
            title(hemName{h});
            for reg = regions
                subset = Ts.regType==reg&Ts.regSide==h;
                traceplot([-pre:post],Ts.y_adj,'errorfcn','stderr',...
                    'split',[Ts.regType],'leg',regname(regions),'subset',subset); % ,
                hold on;
                traceplot([-pre:post],Ts.y_hat,'linestyle','--',...
                    'split',[Ts.regType], 'linewidth',3,'subset',subset); % ,
                %hold off;
            end
            xlabel('TR');
            ylabel('activation');
            
            drawline(0);
        end
        
        % save data
        fitstr = sprintf('%d',fithrf);
        save(fullfile(regDir,sprintf('ROI_timeseries_fit_%d_%s.mat',glm,fitstr)),'T','Ts');
        
        varargout = {T,Ts};
    case 'ROI_timeseries_fit2'                      % Fit the hrf to the data extracted from ROIs
        sn  = varargin{1};
        glm = 9;
        regions = [2];
        hem = 1;
        fithrf = [1,2]'; % hrf parameter to be fitted
        LB     = [2  9   0.2 0.2  3  0 0.1]';
        UB     = [8 17  4   4    12 1.5  4.5]';
        
        Niter = 10;
        
        vararginoptions({varargin{2:end}},{'glm','region','hem','fit','LB','UB'});
        pre     = 4;
        post    = 12;
        
        Ts = [];T = [];
        for s=sn
            subj = subj_name{s};
            fprintf('%s\n',subj);
            
            switch glm
                case 8
                    manualDelay = -5000/2720;
                    delay       = 1+manualDelay;
                    dur         = (1000)/2720;
                    duration    = 1;
                    onsetshift  = delay;
                    P0          = [4.5  16  1   1    6  1.5    1]';
                    
                    D = dload(fullfile(behaviourDir,[subj_name{s}],[subj_name_behav{s},'.dat']));
                    avrgRT      = nanmean(D.RT,1);      % reaction time
                    avrgMT      = nanmean(D.MT,1);      % movement time
                    holdTime    = nanmean(D.time2plan);
                    delay       = 1.5;%(holdTime+avrgRT)/TR; % Delay (for some reason delay is long)
                    dur         = (avrgRT+avrgMT+holdTime-500)/TR;%avrgMT/TR*0.5;
                    %onsetshift  = -delay;
                    P0      = [4.5,16,1,1,6,0.8,dur]';
                case {9,10,11,12}
                    D = dload(fullfile(behaviourDir,[subj_name{s}],[subj_name_behav{s},'.dat']));
                    avrgRT      = nanmean(D.RT,1);      % reaction time
                    avrgMT      = nanmean(D.MT,1);      % movement time
                    holdTime    = nanmean(D.time2plan);
                    delay       = 1.5;%(holdTime+avrgRT)/TR; % Delay (for some reason delay is long)
                    dur         = (avrgRT+avrgMT+holdTime)/TR;%avrgMT/TR*0.5;
                    onsetshift  = -delay+1;
                    P0      = [4,10,1,1,6,0,0]'; % P0 should be given in second?
            end
            
            
            cd(fullfile(baseDir,glmName{glm},subj_name{s}));
            
            load SPM;
            load(fullfile(regDir,subj,[subj '_regAll_movement.mat']));
            
            % Get and tabulate design for pattern consistency
            E   = load('SPM_info.mat');
            indx= find(E.run<11);
            E   = getrow(E,indx);
            Condition = E.movType(indx,1);
            Partition = E.run(indx,1);
            
            % Clearn default onset and duration
            if (1)%ismember(fithrf,[6,7])
                for r = 1:length(SPM.nscan)
                    for u=1:length(SPM.Sess(r).U)
                        %SPM.Sess(r).U(u).dur = ones(size(SPM.Sess(r).U(u).dur))*(1/16);%*duration; % 1
                        SPM.Sess(r).U(u).ons = SPM.Sess(r).U(u).ons+onsetshift; % return to TR at announce trial
                    end;
                    SPM.Sess(r).U=spm_get_ons(SPM,r);
                end;
            end
            
            for h = hem
                for reg = regions
                    if reg<9
                        roi = numregions_surf*(h-1)+reg;
                    elseif reg>=9
                        roi = 2*numregions_surf + numregions_BG*(h-1)+(reg-8);
                    end
                    fprintf('%s\n',R{roi}.name);
                    
                    % Get data
                    Y = region_getdata(SPM.xY.VY,R{roi});
                    numVox = size(R{roi}.data,1);
                    
                    % cut out NaN or all-zero voxels
                    idx = sum(Y.*Y,1)==0;
                    if sum(idx)==0
                        fprintf('... done (%d voxels)\n',numVox);
                    else
                        fprintf('... done (%d voxels, %d were discarded since containing no data)\n',numVox,sum(idx));
                    end
                    Y = Y(:,~idx);
                    numVox = size(Y,2);
                    
                    % get pre-whitened beta and calc pattern consistency
                    B = mva_prewhiten_beta(Y,SPM);
                    R2pre = rsa_patternConsistency(B,Partition,Condition,'removeMean',1);
                    
                    % Check Error before
                    Ypre        = spm_filter(SPM.xX.K,SPM.xX.W*Y);
                    Yres        = spm_sp('r',SPM.xX.xKXs,Ypre);
                    Epre        = sum(sum(Yres.^2))/numel(Yres(:));
                    
                    % Fit a common hrf
                    e = 1; E = [];
                    cons = 0; Cons = [];
                    P0_ = P0; Para = [];
                    iter = 1;
                    fprintf('Fitting hrf...');
                    while ((iter<Niter))%((e >= 1)&&(iter<5))
                        fprintf('iteration no. %d\n',iter);
                        
                        
                        [P,SPM,Yhat,Yres] = spmj_fit_hrf(SPM,Y,...
                            'fit',fithrf,'LB',LB,'UB',UB,'P0',P0_);
                        Para{end+1} = P;
                        % update initial value
                        P0_(fithrf) = P0(fithrf)+0.1*rand*(UB(fithrf)-LB(fithrf));
                        iter = iter+1;
                        
                        % Check Error after
                        Epost       = sum(sum(Yres.^2))/numel(Yres(:));
                        e           = Epost/Epre
                        E           = [E,Epost];
                        % get pre-whitened beta and calc pattern consistency
                        B = mva_prewhiten_beta(Y,SPM);
                        R2post = rsa_patternConsistency(B,Partition,Condition,'removeMean',1);
                        cons   = R2post/R2pre
                        Cons   = [Cons,R2post];
                    end
                    
                    % find best one
                    bestindex1  = find(E==min(E));
                    bestindex2  = find(Cons==max(Cons));
                    bestindex   = unique([bestindex1,bestindex2]);
                    P = Para{bestindex(1)};
                    Epost = E(bestindex(1));
                    R2post = Cons(bestindex(1));
                    
                    disp('Fitted hrf parameters:');
                    fprintf('%2.2f\t',P);
                    fprintf('Epost/Epre: %1.3f\n',e);
                    fprintf('Pattern Consistency post/pre: %1.3f\n',cons);
                    
                    % Parameter values
                    T_.fit = fithrf';
                    T_.P = P';
                    T_.R2 = 1-(trace(Yres*Yres')/(trace(Yres*Yres')+trace(Yhat*Yhat')));
                    T_.Consistency = R2post;
                    T_.ConsistencyGain = R2post/R2pre;
                    T_.Eratio = Epost/Epre;
                    T_.regSide = h;
                    T_.regType = reg;
                    T_.SN = s;
                    T_.hemname = {hemName{h}};
                    T_.regname = {regname{reg}};
                    T = addstruct(T,T_);
                    
                    % Timeseries
                    D = spmj_get_ons_struct(SPM);
                    
                    y_hat = mean(Yhat,2);
                    y_res = mean(Yres,2);
                    for i=1:size(D.block,1);
                        D.y_hat(i,:)=cut(y_hat,pre,round(D.ons(i)),post,'padding','nan')';
                        D.y_res(i,:)=cut(y_res,pre,round(D.ons(i)),post,'padding','nan')';
                        D.y_adj(i,:)=D.y_hat(i,:)+D.y_res(i,:);
                    end;
                    D.regType   = ones(size(D.event,1),1)*reg;
                    D.regSide   = ones(size(D.event,1),1)*h;
                    D.sn        = ones(size(D.event,1),1)*s;
                    Ts          = addstruct(Ts,D);
                end
            end
        end
        
        % plot fitting result
        figure('name',sprintf('hrf fitting results: %s%s',regname{regions(reg)},hemName{hem(h)}),'position',[5 5 25 25]);
        for h=1:numel(hem)
            for reg=1:numel(regions)
                subplot(2,2,1)
                myboxplot([],T.P(:,fithrf),'subset',T.regSide==hem(h)&T.regType==regions(reg));
                set(gca,'xticklabel',fithrf);xlabel('Parameter number');
                title('Fitted parameters');
                
                subplot(2,2,2)
                scatterplot(T.R2,T.Consistency,'identity','label',sn,...
                    'subset',T.regSide==hem(h)&T.regType==regions(reg));
                xlabel('R^2 for 1st-level glm');
                ylabel('R^2 for pattern consistency')
                title('Benefit of fitting: values');
                
                subplot(2,2,3)
                scatterplot(T.Eratio,T.ConsistencyGain,'draworig','label',sn,...
                    'subset',T.regSide==hem(h)&T.regType==regions(reg));hold on
                drawline(1,'dir','horz');drawline(1,'dir','vert');
                xlabel('Error reduction in 1st-level glm');
                ylabel('Gain in consistency')
                title('Benefit of fitting: ratio');
            end
        end
        
        figure('name',[what,'_timeseries'])
        for h=hem
            subplot(1,2,h)
            title(hemName{h});
            for reg = regions
                subset = Ts.regType==reg&Ts.regSide==h;
                traceplot([-pre:post],Ts.y_adj,'errorfcn','stderr',...
                    'split',[Ts.regType],'leg',regname(regions),'subset',subset); % ,
                hold on;
                traceplot([-pre:post],Ts.y_hat,'linestyle','--',...
                    'split',[Ts.regType], 'linewidth',3,'subset',subset); % ,
                %hold off;
            end
            xlabel('TR');
            ylabel('activation');
            
            drawline(0);
        end
        
        % save data
        fitstr = sprintf('%d',fithrf);
        save(fullfile(regDir,sprintf('ROI_timeseries_fit2_%d_%s.mat',glm,fitstr)),'T','Ts');
        
        varargout = {T,Ts};
    case 'ROI_timeseries_comphemi'                  % Compare BOLD suppression for unimanual-ipsilateral trials across hemispheres
        sn = varargin{1};
        glm = 10;
        ROI = 'all';
        regions = [1:8];
        hemis = [1:2];
        pre = 10;
        post = 10;
        dt = 20;
        smooth = 'on';
        
        vararginoptions(varargin(2:end),{'regions','hemis','pre','post'})
        
        types = {[1:6],[7:12],[13:48]}; % uni-left, uni-right, bi
        
        % load data
        T = load(fullfile(regDir,sprintf('hrf_%s_%d.mat',ROI,glm)));
        T = getrow(T,ismember(T.regType,regions));
        T.newtype = T.type;
        for t=1:3
            T.newtype(ismember(T.type,types{t})) = t;
        end
        T = tapply(T,{'SN','newtype','regSide','regType'},{'y_adj','nanmean(x,1)','name','y_adj'});
        
        task = [2:5]+pre+1;
        
        T.y_task = nanmean(T.y_adj(:,task),2);
        T.y_base = nanmean(T.y_adj(:,[1:task(1)-1]),2);
        
        % smooth ?
        switch smooth
            case 'on'
                y_smooth = zeros(size(T.y_adj,1),TR/dt*(size(T.y_adj,2)));
                y_smooth(:,1:TR/dt:end) = T.y_adj;
                t = [-TR:dt:TR]; % msec
                sigma = 1500;
                kernel = exp(-0.5*(t/sigma).^2);
                kernel = kernel/max(kernel);
                for t = 1:size(T.SN,1)
                    T.y_smooth(t,:) = conv(y_smooth(t,:),kernel,'same');
                end
                time = [1:size(y_smooth,2)]*dt-TR*pre;
            case 'off'
                time = [-pre:post];
                T.y_smooth = T.y_adj;
        end
        
        % loop over trial type
        figure('name',sprintf('Timeseries-compare-suppression-across-hemisphere%s',glmName{glm}),...
            'position',[10 10 2.3*sqrt(2)*10 10]);
        
        for h=hemis;
            subset = T.regSide==h & ismember(T.SN,sn) & ismember(T.newtype,h);
            
            % draw filtered BOLD responses for
            % unimanual-ipsilateral trials
            subplot(1,2,h);
            switch smooth
                case 'on'
                    xtick = [find(time==-3*trialTime),find(time==-2*trialTime),find(time==-trialTime),...
                        find(time==0),...
                        find(time==2*TR),find(time==5*TR),...
                        find(time==3*trialTime)];
                    ydata = pivottablerow(T.regType,T.y_smooth,'nanmean(x,1)',...
                        'subset',subset);
                    imagesc_rectangle(ydata,'YDir','reverse','MAP',jet(),...
                        'scale',[-2.5 2.5],'overlap',0.95);caxis([-2.5 2.5]);
                    hold on;
                    colormap(jet);
                    colorbar();
                    drawline([find(time==0)],...
                        'dir','vert','color','k',...
                        'linewidth',1,'linestyle','-');
                    drawline([find(time==2*TR),find(time==5*TR)],...
                        'dir','vert','color','k',...
                        'linewidth',1,'linestyle',':');
                    xlabel('Time (sec)');ylabel('ROI');
                    set(gca,'xtick',sort(xtick),'xticklabel',{'-21','-14','-7','0','5.4','13.5','21'});
                    set(gca,'ytick',[regions],'yticklabel',regname(regions),'tickdir','out');
                case 'off'
                    traceplot(time,T.y_smooth,'errorfcn','stderr',...
                        'split',[T.regType],'subset',subset,...
                        'leg',regname(regions),'leglocation','bestoutside'); % ,
                    hold on;
                    drawline([0],'dir','horz','linestyle','-');
                    drawline([0],'dir','vert','linestyle','-');
                    hold off;
                    xlabel('TR');
                    ylabel('Activation (a.u.)');
                    
            end
            title(hemName{h})
            
        end
        
        % compare across hemispheres and rois
        P = [];Anova = [];
        for reg=[1:5,7,8];
            idxA = T.regType==reg&T.regSide==1&T.newtype==1;
            idxB = T.regType==reg&T.regSide==2&T.newtype==2;
            
            A = T.y_base(idxA)-T.y_task(idxA);
            B = T.y_base(idxB)-T.y_task(idxB);
            
            % ttest
            [~,p] = ttest_mc(A,B,2,'paired');
            P = [P;p];
            
            ano.data = [A;B];
            ano.SN = [T.SN(idxA);T.SN(idxB)];
            ano.hemi = [T.regSide(idxA);T.regSide(idxB)];
            ano.region = [T.regType(idxA);T.regType(idxB)];
            Anova = addstruct(Anova,ano);
        end
        % 2-way rmanova
        anovaMixed(Anova.data,Anova.SN,'within',[Anova.hemi,Anova.region],{'hemi','region'});
        
        % bar graph
        figure('position',[5 5 sqrt(2)*10 10]);
        subplot(1,2,1)
        barplot([T.regType],[T.y_base-T.y_task],...
            'subset',T.regType~=6&T.regSide==1&T.newtype==1,...
            'style_bold','facecolor',[0.7 0.7 0.7]);
        set(gca,'ylim',[-0.05 2],'xticklabel',{});
        set(gca,'View',[90 90],'Xaxislocation','top','YDir','reverse','tickdir','out');
        ylabel('Suppression (a.u.)');title('Left hemisphere')
        
        subplot(1,2,2)
        barplot([T.regType],[T.y_base-T.y_task],...
            'subset',T.regType~=6&T.regSide==2&T.newtype==2,...
            'style_bold','facecolor',[0.7 0.7 0.7]);
        set(gca,'ylim',[-0.05 2],'xticklabel',regname([1:5,7,8]));
        set(gca,'view',[90 90],'tickdir','out');
        ylabel('Suppression (a.u.)');title('Right hemisphere')
        
    case 'ROI_beta_getMean'                         % Extracts ROI data and calc mean beta value
        sn      = varargin{1};
        glm     = 10;
        ROI     = 'all';
        fname 	= 'movement';
        regions = 1:12;
        hemi    = 1:2;
        
        vararginoptions(varargin(2:end),{'glm','ROI','fname','regions','hemi','fig'});
        
        for g=glm
            T = [];
            for s=sn
                fprintf('subj = %d\n',s)
                
                % load SPM.mat
                glmDirSubj = fullfile(baseDir,glmName{g},subj_name{s});
                load(fullfile(glmDirSubj,'SPM.mat'));
                %SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
                cd(glmDirSubj);
                
                % choose functionaly/anatomical ROI
                switch (ROI)
                    case 'func'
                        load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_reg800_%s_%d.mat',fname,glm)]));
                    case 'all'
                        load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_regAll_%s.mat',fname)]));
                end
                
                % Get and tabulate design
                E   = load(fullfile(glmDirSubj,'SPM_info.mat'));
                % calc bimanual contrast for intrinsic,extrinsic,asymmetric
                % movements
                extrinsic = [0:60:300];
                intrinsic = [180,120,60,0,300,240];
                for t=1:length(E.hand)
                    if E.u_or_b(t,1)==1
                        idx = find(extrinsic==E.dirL(t,1));
                        if E.dirR(t,1)==extrinsic(idx)
                            BiCat(t,1) = 1;
                        elseif E.dirR(t,1)==intrinsic(idx)
                            BiCat(t,1) = 2;
                        else
                            BiCat(t,1) = 3;
                        end
                    else
                        BiCat(t,1) = 0;
                    end
                end
                
                % Now get betas from these regions
                NRegions = numel(regions)*numel(hemi);
                nRegions = numel(R);
                for reg = regions
                    for hem = hemi
                        if reg<9
                            roi = numregions_surf*(hem-1)+reg;
                        elseif reg>=9
                            roi = 2*numregions_surf + numregions_BG*(hem-1)+(reg-8);
                        end
                        if (nRegions<NRegions)
                            warning('number of ROIs doesn''t much!');
                            break;
                        end
                        if (~isempty(R{roi}))
                            fprintf('extracting: %s ...',R{roi}.name);
                            
                            % get data
                            data = region_getdata(SPM.Vbeta,R{roi});
                            numVox = size(R{roi}.data,1);
                            % cut out NaN or all-zero voxels
                            idx = sum(data.*data,1)==0|all(isnan(data));
                            if sum(idx)==0
                                fprintf('... done (%d voxels)\n',numVox);
                            else
                                fprintf('... done (%d voxels, %d were discarded since containing no data)\n',numVox,sum(idx));
                            end
                            data = data(:,~idx);
                            numVox = size(data,2);
                            
                            % calc mean beta value (mean activation)
                            %-- uni left
                            S.UL = mean(mean(data(E.hand==0,:)));
                            %-- uni right
                            S.UR = mean(mean(data(E.hand==1,:)));
                            %-- bi symmetric extrinsic
                            S.BEx = mean(mean(data(E.hand==2&BiCat==1,:)));
                            %-- bi symmetric intrinsic
                            S.BIn = mean(mean(data(E.hand==2&BiCat==2,:)));
                            %-- bi asymmetric
                            S.BUn = mean(mean(data(E.hand==2&BiCat==3,:)));
                            
                            S.region    = reg;
                            S.subj      = s;
                            S.hemis     = hem;
                            S.numVox    = numVox;
                            S.glm       = g;
                            
                            T = addstruct(T,S);
                        end;
                    end;
                end;
            end;
            varargout = {T};
            save(fullfile(regDir,sprintf('meanBeta.glm%d.mat',g)),'-struct','T');
        end
    case 'ROI_beta_getPD'                           % Extracts ROI data and calc mean beta value
        sn      = varargin{1};
        glm     = 10;
        ROI     = 'all';
        fname 	= 'movement';
        regions = 1:12;
        hemi    = 1:2;
        
        vararginoptions(varargin(2:end),{'glm','ROI','fname','regions','hemi','fig'});
        
        for g=glm
            T = [];
            for s=sn
                fprintf('subj = %d\n',s)
                
                % load SPM.mat
                glmDirSubj = fullfile(baseDir,glmName{g},subj_name{s});
                load(fullfile(glmDirSubj,'SPM.mat'));
                %SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
                cd(glmDirSubj);
                
                % choose functionaly/anatomical ROI
                switch (ROI)
                    case 'func'
                        load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_reg800_%s_%d.mat',fname,glm)]));
                    case 'all'
                        load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_regAll_%s.mat',fname)]));
                end
                
                % Get and tabulate design
                E   = load(fullfile(glmDirSubj,'SPM_info.mat'));
                % calc bimanual contrast for intrinsic,extrinsic,asymmetric
                % movements
                extrinsic = [0:60:300];
                intrinsic = [180,120,60,0,300,240];
                for t=1:length(E.hand)
                    if E.u_or_b(t,1)==1
                        idx = find(extrinsic==E.dirL(t,1));
                        if E.dirR(t,1)==extrinsic(idx)
                            BiCat(t,1) = 1;
                        elseif E.dirR(t,1)==intrinsic(idx)
                            BiCat(t,1) = 2;
                        else
                            BiCat(t,1) = 3;
                        end
                    else
                        BiCat(t,1) = 0;
                    end
                end
                
                % Now get betas from these regions
                NRegions = numel(regions)*numel(hemi);
                nRegions = numel(R);
                for reg = regions
                    for hem = hemi
                        if reg<9
                            roi = numregions_surf*(hem-1)+reg;
                        elseif reg>=9
                            roi = 2*numregions_surf + numregions_BG*(hem-1)+(reg-8);
                        end
                        if (nRegions<NRegions)
                            warning('number of ROIs doesn''t much!');
                            break;
                        end
                        if (~isempty(R{roi}))
                            fprintf('extracting: %s ...',R{roi}.name);
                            
                            % get data
                            data = region_getdata(SPM.Vbeta,R{roi});
                            data = data(1:length(E.sn),:); % discard intercept
                            numVox = size(R{roi}.data,1);
                            % cut out NaN or all-zero voxels
                            idx = sum(data.*data,1)==0|all(isnan(data));
                            if sum(idx)==0
                                fprintf('... done (%d voxels)\n',numVox);
                            else
                                fprintf('... done (%d voxels, %d were discarded since containing no data)\n',numVox,sum(idx));
                            end
                            data = data(:,~idx);
                            numVox = size(data,2);
                            
                            % average over run (resulting in KxP beta matrix)
                            B = E;
                            B.beta = data; clear data
                            B = tapply(B,{'hand','dirL','dirR'},{'beta','nanmean(x,1)','name','beta'});
                            
                            % get PD
                            UL = getrow(B,B.hand==0);
                            UR = getrow(B,B.hand==1);
                            Bi = getrow(B,B.hand==2);
                            BL = tapply(Bi,{'dirL'},{'beta','nanmean(x,1)','name','beta'});
                            BR = tapply(Bi,{'dirR'},{'beta','nanmean(x,1)','name','beta'});
                            
                            [~,pdul] = max(UL.beta,[],1);
                            [~,pdur] = max(UR.beta,[],1);
                            [~,pdbl] = max(BL.beta,[],1);
                            [~,pdbr] = max(BR.beta,[],1);
                            dpdl = pdul-pdbl;
                            dpdr = pdur-pdbr;
                            dpdl = int8(acosd(cosd(dpdl*60))/60);
                            dpdr = int8(acosd(cosd(dpdr*60))/60);
                            
                            [S.UL,S.bin]   = histcounts(pdul,[1:7]);
                            [S.UR]   = histcounts(pdur,[1:7]);
                            [S.BL]  = histcounts(pdbl,[1:7]);
                            [S.BR]  = histcounts(pdbr,[1:7]);
                            [S.dL]  = histcounts(dpdl,[0:4]);
                            [S.dR,S.bind]  = histcounts(dpdr,[0:4]);
                            
                            S.region    = reg;
                            S.subj      = s;
                            S.hemis     = hem;
                            S.numVox    = numVox;
                            S.glm       = g;
                            
                            T = addstruct(T,S);
                        end;
                    end;
                end;
            end;
            varargout = {T};
            save(fullfile(regDir,sprintf('PD.glm%d.mat',g)),'-struct','T');
        end
    case 'ROI_beta_plotPD'
        glm = 10;
        
        % load PD data
        T = load(fullfile(regDir,sprintf('PD.glm%d.mat',glm)));
        
        subj = unique(T.subj)';
        hemi = unique(T.hemis)';
        regions = unique(T.region)';
        regions = 1:8;
        
        
        % mix hemispheres to assess contra and ipsi representation
        L = getrow(T,T.hemis==1);
        R = getrow(T,T.hemis==2);
        L.dcontra = L.dR;
        L.dipsi = L.dL;
        R.dcontra = R.dL;
        R.dipsi = R.dR;
        All = addstruct(L,R);
        All.dcontra = bsxfun(@rdivide,All.dcontra,All.numVox)*100;
        All.dipsi = bsxfun(@rdivide,All.dipsi,All.numVox)*100;
        
        % plot delta PD
        figure('position',[5 5 25 15]);c=1;
        for reg=regions
            A = getrow(All,All.region==reg);
            
            % compare
            for bin=1:4
                [t,p(bin)] = ttest_mc(A.dcontra(:,bin),A.dipsi(:,bin),2,'paired');
            end
            
            ipsi = vec(A.dipsi');
            contra = vec(A.dcontra');
            bin = vec(A.bind(:,1:4)');
            
            subplot(2,4,reg); %
            if c==1
                xpos = lineplot(bin,[contra,ipsi],'style_shade',...
                    'linecolor',{darkgold,darkpurple},...
                    'shadecolor',{gold,purple},'leg',{'Contra','Ipsi'},...
                    'leglocation','best'); hold on
                title(regname{reg});
                xlabel('|PD_u - PD_b| (degree)');
                ylabel('%voxel');
            else
                xpos = lineplot(bin,[contra,ipsi],'style_shade',...
                    'linecolor',{darkgold,darkpurple},...
                    'shadecolor',{gold,purple}); hold on
                title(regname{reg});
                xlabel('|PD_u - PD_b| (degree)');
                ylabel('%voxel');
            end
            set(gca,'xtick',[xpos],'xticklabel',{'0','60','120','180'},...
                'tickdir','out','ylim',[10 42],'ytick',[10:10:40]);
            
            % plot p-value
            for bin=1:4
                if p(bin)<0.05
                    text(xpos(bin),40,sprintf('%1.4f',p(bin)),'color','r');
                end
            end
            
            c=c+1;
        end
        
    case 'ROI_beta_plotMean'
        
        glm = 10;
        vararginoptions(varargin,{'glm'});
        
        % load data
        T = load(fullfile(regDir,sprintf('meanBeta.glm%d.mat',glm)));
        
        hemi = unique(T.hemis)';
        regions = unique(T.region)';
        sn = unique(T.subj)';
        
        % plot
        figure('position',[5 5 15 10]);
        for h=hemi
            subplot(2,2,h);
            barplot(T.region,[T.UL,T.UR],...
                'subset',T.hemis==h,'style_bold','facecolor',[0.7 0.7 0.7]);
            %set(gca,'xTicklabel',regname(regions));
            ylabel('Mean activation (a.u.)');
        end
        for h=hemi
            subplot(2,2,h+2);
            barplot(T.region,[T.BEx,T.BIn,T.BUn],...
                'subset',T.hemis==h,'style_bold','facecolor',[0.7 0.7 0.7]);
            %set(gca,'xTicklabel',regname(regions));
            ylabel('Mean activation (a.u.)');
        end
    case 'ROI_patternConsistency'                   % Extracts ROI data and calculate pattern consistency using rsa_patternConsistency.m
        T=[];
        sn = varargin{1};
        glms = varargin{2};
        ROI = 'all';
        fname 	= 'movement';       % using metric file which compares any movement against rest
        regions = 1:8;
        hemi = 1:2;
        vararginoptions({varargin{3:end}},{'glm','ROI','fname','regions','hemi'});
        
        Data=[];
        for glm=glms
            for s=sn
                fprintf('subj = %d\n',s)
                
                % load SPM.mat
                glmDirSubj=fullfile(baseDir,glmName{glm},subj_name{s});
                load(fullfile(glmDirSubj,'SPM.mat'));
                
                %SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
                
                % choose ROI
                switch (ROI)
                    case 'func'
                        load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_reg800_%s_%d.mat',fname,glm)]));
                    case 'all'
                        load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_regAll_%s.mat',fname)]));
                end
                
                % Get and tabulate design
                E   = load(fullfile(glmDirSubj,'SPM_info.mat'));
                indx= find(E.run<11);
                E   = getrow(E,indx);
                Condition = E.movType(indx,1);
                Partition = E.run(indx,1);
                
                % Now get prewhitened betas from these regions
                NRegions = numel(regions)*numel(hemi);
                nRegions = numel(R);
                for reg = regions
                    for hem = hemi
                        if reg<9
                            roi = numregions_surf*(hem-1)+reg;
                        elseif reg>=9
                            roi = 2*numregions_surf + numregions_BG*(hem-1)+(reg-8);
                        end
                        if (nRegions<NRegions)
                            warning('number of ROIs doesn''t much!');
                            break;
                        end
                        if (~isempty(R{roi}))
                            fprintf('extracting: %s ...',R{roi}.name);
                            
                            % get data
                            data = region_getdata(SPM.xY.VY,R{roi});
                            numVox = size(R{roi}.data,1);
                            % cut out NaN or all-zero voxels
                            idx = sum(data.*data,1)==0;
                            if sum(idx)==0
                                fprintf('... done (%d voxels)\n',numVox);
                            else
                                fprintf('... done (%d voxels, %d were discarded since containing no data)\n',numVox,sum(idx));
                            end
                            data = data(:,~idx);
                            numVox = size(data,2);
                            
                            % get pre-whitened beta
                            %B = mva_prewhiten_beta(data,SPM);
                            [B]=rsa_noiseNormalizeBeta(data,SPM,'normmode','overall');
                            
                            % calc pattern consistency
                            R2 = rsa_patternConsistency(B,Partition,Condition,'removeMean',1);
                            R2_ = rsa_patternConsistency(B,Partition,Condition,'removeMean',0);
                            
                            S.ConsR2_womean = R2;
                            S.ConsR2_wmean = R2_;
                            S.region = reg;
                            S.subj  = s;
                            S.hemis = hem;
                            S.numVox = numVox;
                            
                            T=addstruct(T,S);
                        end;
                    end;
                end;
            end;
            varargout={T};
            save(fullfile(regDir,sprintf('ROI_patternConsistency.glm%d.mat',glm)),'-struct','T');
        end
    case 'ROI_compare_patternConsistency'
        glm     = varargin{1};
        regions = [1:8];
        hemis   = [1,2];
        
        T = [];        
        glmnames = glmtypes(glm);
        % load consistency data
        for g=glm;
            C = load(fullfile(regDir,sprintf('ROI_patternConsistency.glm%d.mat',g)));
            
            C.glm = repmat(g,size(C.subj));
            
            T = addstruct(T,C);
        end
        
        % compare
        fields = {{'ConsR2_wmean','ConsR2_womean'}};
        newfields = {{'R2withMean','R2withoutMean'}};
        for f=1:numel(fields)
            T.(newfields{f}{1}) = T.(fields{f}{1});
            T.(newfields{f}{2}) = T.(fields{f}{2});
        end
        for f=1:numel(fields)
            for h=hemis
                figure('name',sprintf('Pattern consistency (%s)',hemName{h}),...
                    'position',[5 5 40 20]);
                
                subplot(2,1,1); % consistency calculated with mean pattern
                x=barplot(T.region,T.(newfields{f}{1}),'split',T.glm,...
                    'leg',glmnames,'leglocation','eastoutside',...
                    'subset',T.hemis==h);
                set(gca,'XTick',x(1:numel(glm):end),'Xticklabel',regname(regions));
                %set(gca,'Ylim',[0.08 0.19]);
                title(newfields{f}{1})
                
                subplot(2,1,2); % consistency calculated after mean pattern subtraction
                x=barplot(T.region,T.(newfields{f}{2}),'split',T.glm,...
                    'leg',glmnames,'leglocation','eastoutside',...
                    'subset',T.hemis==h);
                set(gca,'XTick',x(1:numel(glm):end),'Xticklabel',regname(regions));
                %set(gca,'Ylim',[0.08 0.13]);
                title(newfields{f}{2})
            end
        end
        
        varargout = {};
        
    case 'ROI_get_prewhitened_beta'             % Extracts prewhitened beta
        T=[];
        sn = varargin{1};
        glm = varargin{2}; % 10
        ROI = 'all';
        fname 	= 'movement';       % using metric file which compares any movement against rest
        regions = 1:8;
        hemi = 1:2;
        append = 0;
        vararginoptions({varargin{3:end}},{'glm','ROI','fname','regions','hemi','append'});
        
        if append == 0
            Data=[];
            Design=[];
        elseif append == 1
            % load Data and Design
            load(fullfile(regDir,sprintf('ROI_pwhBeta.glm%d.mat',glm)));
        end
        for s=sn
            fprintf('subj = %d\n',s)
            
            % load SPM.mat
            glmDirSubj=fullfile(baseDir,glmName{glm},subj_name{s});
            load(fullfile(glmDirSubj,'SPM.mat'));
            
            SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
            
            
            % choose ROI
            switch (ROI)
                case 'func'
                    load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_reg800_%s_%d.mat',fname,glm)]));
                case 'all'
                    load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_regAll_%s.mat',fname)]));
            end
            
            % Get and tabulate design
            E   = load(fullfile(glmDirSubj,'SPM_info.mat'));
            indx= find(E.run<11);
            E   = getrow(E,indx);
            % Movement
            Movement = E.movType(indx,1);
            
            % Partition
            Partition   = E.run(indx,1);
            
            % Condition (uni left/uni right/bimanual)
            Condition = E.hand(indx,1);
            
            % save design
            design.condition    = Movement';
            design.hand         = Condition';
            design.partition    = Partition';
            design.subj         = s;
            Design              = addstruct(Design,design);
            
            % Get prewhitened betas from these regions
            Nconditions = length(Movement);
            NRegions = numel(regions)*numel(hemi);
            nRegions = numel(R);
            for reg = regions
                for hem = hemi
                    if reg<9
                        roi = numregions_surf*(hem-1)+reg;
                    elseif reg>=9
                        roi = 2*numregions_surf + numregions_BG*(hem-1)+(reg-8);
                    end
                    if (nRegions<NRegions)
                        warning('number of ROIs doesn''t much!');
                        break;
                    end
                    if (~isempty(R{roi}))
                        fprintf('extracting: %s ...',R{roi}.name);
                        
                        % get data
                        data = region_getdata(SPM.xY.VY,R{roi});
                        numVox = size(R{roi}.data,1);
                        % cut out NaN or all-zero voxels
                        idx = sum(data.*data,1)==0;
                        if sum(idx)==0
                            fprintf('... done (%d voxels)\n',numVox);
                        else
                            fprintf('... done (%d voxels, %d were discarded since containing no data)\n',numVox,sum(idx));
                        end
                        data = data(:,~idx);
                        numVox = size(data,2);
                        
                        % Get pre-whitened beta
                        %B = mva_prewhiten_beta(data,SPM);
                        [B,~,Sw] = rsa.spm.noiseNormalizeBeta(data,SPM);
                        
                        % prewhitened beta
                        S.beta              = {B(1:Nconditions,:)};
                        S.beta_nointerest   = {B(Nconditions+1:end,:)};
                        %S.Sw                = {Sw};
                        
                        % other info
                        S.glm = glm;
                        S.region = reg;
                        S.subj  = s;
                        S.hemis = hem;
                        S.numVox = numVox;
                        
                        Data=addstruct(Data,S);
                    end;
                end;
            end;
        end;
        save(fullfile(regDir,sprintf('ROI_pwhBeta.glm%d.mat',glm)),'Data','Design');
        varargout={Data,Design};
        
    case 'ROI_extract_prewhitened_data_group'       % extract ROI data and save prewhitend data. good way of saving time for prewhitening
        %data.(subjname).(roiname) = N x P matrix
        sn = varargin{1};
        chunkset = 'A';
        glm = 1;
        ROI = 'all';
        fname 	= 'movement';       % using metric file which compares any movement against rest
        regions = 1:12;
        hemi = 1:2;
        vararginoptions({varargin{2:end}},{'chunkset','glm','ROI','fname','regions','hemi'});
        
        Data=[];
        for s=sn
            fprintf('subj = %d\n',s)
            
            % load SPM.mat
            glmDirSubj=fullfile(baseDir,glmName{glm},subj_name{s});
            load(fullfile(glmDirSubj,'SPM.mat'));
            %SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
            
            % choose ROI
            switch (ROI)
                case 'func'
                    load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_reg800_%s_%d.mat',fname,glm)]));
                case 'all'
                    load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_regAll_%s_%d.mat',fname,glm)]));
            end
            
            % Get and tabulate design
            E   = load(fullfile(glmDirSubj,'SPM_info.mat'));
            indx= find(E.run<11);
            Design.(subj_name{s}).condition = E.seqType(indx,1);
            Design.(subj_name{s}).partition = E.run(indx,1);
            
            % Now get prewhitened betas from these regions
            NRegions = numel(regions)*numel(hemi);
            nRegions = numel(R);
            for reg = regions
                for hem = hemi
                    if reg<9
                        roi = numregions_surf*(hem-1)+reg;
                    elseif reg>=9
                        roi = 2*numregions_surf + numregions_BG*(hem-1)+(reg-8);
                    end
                    if (nRegions<NRegions)
                        warning('number of ROIs doesn''t much!');
                        break;
                    end
                    if (~isempty(R{roi}))
                        fprintf('extracting: %s ...',R{roi}.name);
                        
                        % get data
                        data = region_getdata(SPM.xY.VY,R{roi});
                        
                        % cut out NaN or all-zero voxels
                        idx = sum(data.*data,1)==0;
                        data = data(:,~idx);
                        
                        % pre-whiten data
                        B = mva_prewhiten_beta(data,SPM);
                        Data.(subj_name{s}).(R{roi}.name) = B(indx,:);
                        [N,P] = size(B);
                        
                        if sum(idx)==0
                            fprintf('... done (%d voxels)\n',P);
                        else
                            fprintf('... done (%d voxels, %d were discarded since containing no data)\n',P,sum(idx));
                        end
                    end;
                end;
            end;
        end;
        varargout={Data,Design};
        save(fullfile(regDir,sprintf('set%s_ROI_prewhiten_data.mat',chunkset)),'Design','Data');
    case 'ROI_extract_prewhitened_data_individ'     % extract ROI data and save prewhitend data. good way of saving time for prewhitening
        sn = varargin{1};
        glm = 1;
        ROI = 'all';
        fname 	= 'movement';       % using metric file which compares any movement against rest
        regions = 1:8;
        hemi = 1:2;
        vararginoptions(varargin(2:end),{'glm','ROI','fname','regions','hemi'});
        
        % loop over subjects
        for s=sn
            Data = [];
            subj = subj_name{s};
            fprintf('subj = %s\n',subj)
            
            % load SPM.mat
            glmDirSubj=fullfile(baseDir,glmName{glm},subj_name{s});
            load(fullfile(glmDirSubj,'SPM.mat'));
            %SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
            
            % choose ROI
            switch (ROI)
                case 'all'
                    load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_regAll_%s.mat',fname)]));
                otherwise
                    load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_reg%s_%s_%d.mat',ROI,fname,glm)]));
                    
            end
            
            % Get and tabulate design
            E   = load(fullfile(glmDirSubj,'SPM_info.mat'));
            indx= find(E.run<11);
            
            % Save into Design
            Design.condition    = E.movType(indx,1)';
            Design.partition    = E.run(indx,1)';
            Design.SN           = s;
            Design.name         = subj;
            
            % Now get prewhitened betas from these regions
            NRegions = numel(regions)*numel(hemi);
            nRegions = numel(R);
            for reg = regions % diogo changed this
                for hem = hemi
                    if reg<9     % ROI cortical
                        roi = numregions_surf*(hem-1)+reg;
                    elseif reg>=9 % remaining ROI's (BG)
                        roi = 2*numregions_surf + numregions_BG*(hem-1)+(reg-8);
                    end
                    if (nRegions<NRegions)
                        warning('number of ROIs doesn''t much!');
                        break;
                    end
                    if (~isempty(R{roi}))
                        fprintf('extracting: %s ...',R{roi}.name);
                        
                        % get data
                        data = region_getdata(SPM.xY.VY,R{roi});
                        
                        % cut out NaN or all-zero voxels
                        idx = sum(data.*data,1)==0;
                        data = data(:,~idx);
                        
                        % pre-whiten data
                        B = mva_prewhiten_beta(data,SPM);
                        % [B] = rsa.spm.noiseNormalizeBeta(data,SPM) % new
                        
                        [N,P] = size(B);
                        
                        % Save into Data
                        Data_.Beta   = {B(indx,:)};
                        Data_.SN     = s;
                        Data_.region = roi;
                        Data_.name   = {R{roi}.name};
                        Data_.nVox   = P;
                        
                        Data         = addstruct(Data,Data_);
                        
                        if sum(idx)==0
                            fprintf('... done (%d voxels)\n',P);
                        else
                            fprintf('... done (%d voxels, %d were discarded since containing no data)\n',P,sum(idx));
                        end
                    end;
                end;
            end;
            save(fullfile(regDir,subj,sprintf('%s_prewhitened_data_%s_%d.mat',subj,ROI,glm)),'Design','Data');
            size(Data.Beta)
        end;
        varargout={};
    case 'ROI_extract_prewhitened_data_individ_rsa' % extract ROI data and save prewhitend data. good way of saving time for prewhitening
        sn = varargin{1};
        glm = 1;
        ROI = 'all';
        fname 	= 'movement';       % using metric file which compares any movement against rest
        regions = 1:8;
        hemi = 1:2;
        vararginoptions(varargin(2:end),{'glm','ROI','fname','regions','hemi'});
        
        % loop over subjects
        for s=sn
            Data = [];
            subj = subj_name{s};
            fprintf('subj = %s\n',subj)
            
            % load SPM.mat
            glmDirSubj=fullfile(baseDir,glmName{glm},subj_name{s});
            load(fullfile(glmDirSubj,'SPM.mat'));
            %SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
            % save KX
            KX = SPM.xX.xKXs.X;
            
            % choose ROI
            switch (ROI)
                case 'all'
                    load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_regAll_%s.mat',fname)]));
                otherwise
                    load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_reg%s_%s_%d.mat',ROI,fname,glm)]));
                    
            end
            
            % Get and tabulate design
            E   = load(fullfile(glmDirSubj,'SPM_info.mat'));
            indx= find(E.run<11);
            
            % Save into Design
            Design.condition    = E.movType(indx,1)';
            Design.partition    = E.run(indx,1)';
            Design.SN           = s;
            Design.name         = subj;
            
            % Now get prewhitened betas from these regions
            NRegions = numel(regions)*numel(hemi);
            nRegions = numel(R);
            for reg = regions % diogo changed this
                for hem = hemi
                    if reg<9     % ROI cortical
                        roi = numregions_surf*(hem-1)+reg;
                    elseif reg>=9 % remaining ROI's (BG)
                        roi = 2*numregions_surf + numregions_BG*(hem-1)+(reg-8);
                    end
                    if (nRegions<NRegions)
                        warning('number of ROIs doesn''t much!');
                        break;
                    end
                    if (~isempty(R{roi}))
                        fprintf('extracting: %s ...',R{roi}.name);
                        
                        % get data
                        data = region_getdata(SPM.xY.VY,R{roi});
                        
                        % cut out NaN or all-zero voxels
                        idx = sum(data.*data,1)==0;
                        data = data(:,~idx);
                        
                        % pre-whiten data
                        % B = mva_prewhiten_beta(data,SPM);
                        B = rsa.spm.noiseNormalizeBeta(data,SPM); % new
                        % this preserves betas for regressors of
                        % nointerests (i.e., run effect, etc.)
                        
                        [N,P] = size(B);
                        
                        % Save into Data
                        Data_.Beta   = {B};%{B(indx,:)};
                        Data_.SN     = s;
                        Data_.region = roi;
                        Data_.name   = {R{roi}.name};
                        Data_.nVox   = P;
                        
                        Data         = addstruct(Data,Data_);
                        
                        if sum(idx)==0
                            fprintf('... done (%d voxels)\n',P);
                        else
                            fprintf('... done (%d voxels, %d were discarded since containing no data)\n',P,sum(idx));
                        end
                    end;
                end;
            end;
            save(fullfile(regDir,subj,sprintf('%s_prewhitened_data_rsa_%s_%d.mat',subj,ROI,glm)),'Design','Data','KX');
            size(Data.Beta)
        end;
        varargout={};
    case 'ROI_concat_prewhitened_data'              % Concatenate prewhitened data into one data frame
        sn = varargin{1};
        ROI = 'all';
        vararginoptions(varargin(2:end),{'ROI'});
        
        % Loop over subjects
        allDesign = [];
        allData = [];
        for s = sn
            subj = subj_name{s}
            load(fullfile(regDir,subj,sprintf('%s_prewhitened_data_%s.mat',subj,ROI)));
            
            Design.condition = {Design.condition};
            Design.partition = {Design.partition};
            
            allDesign = addstruct(allDesign,Design);
            allData   = addstruct(allData,Data);
        end
        % Save data
        save(fullfile(regDir,sprintf('Allsub_prewhitened_data_%s.mat',ROI)),'allDesign','allData');
        varargout = {allDesign,allData};
    case 'ROI_get_RDMs_toshare'                     % calc RDM for each ROI and save in rsatoolbox format
        sn 		= varargin{1};
        Nsimu   = 1;
        SVD     = 0;
        demean  = 0;
        Nvox    = 160;
        ROI     = 'all';
        Pexclude= []; % partition to be excluded
        nDirL = 8;
        nDirR = 8;
        options = {'Nsimu','SVD','ROI','demean','nDirL','nDirR','Pexclude'};
        vararginoptions({varargin{2:end}},options);
        
        nConditions = nDirL+nDirR+nDirL*nDirR;
        % record options
        for op=1:numel(options)
            Options.(options{op}) = eval(options{op});
        end
        Options.sn = sn;
        
        % get condition no.
        nBi = nDirL*nDirR;
        nConditons = nDirL+nDirR+nBi;
        Conditions = [1:nConditions];
        
        lG 		= nConditions^2; % length of G vector
        ldist 	= nConditions*(nConditions-1)/2; % length of d vector
        
        % Run ROI analysis
        for s = sn
            RDMs = [];
            subj = subj_name{s};
            fprintf('subj = %s\n',subj)
            load(fullfile(regDir,subj,sprintf('%s_prewhitened_data_%s.mat',subj,ROI)))
            
            Partition   = Design.partition;
            Condition   = Design.condition;
            Con     = indicatorMatrix('allpairs',unique(Condition));
            
            
            % loop over rois
            for roi = 1:numel(Data.name)
                % get roi data
                roiname     = Data.name{roi};
                roiname     = strsplit(roiname,'_');
                roiname     = roiname(2:end)
                data_roi    = Data.Beta{roi};
                [N,P]       = size(data_roi);
                
                % loop over multiple iterations (random sub-space approach)
                for i = 1:Nsimu
                    switch (Nsimu)
                        case 1
                            y = data_roi;
                            Nvox = P;
                        otherwise
                            % get random sub-space of ROI (non-nan)
                            voxind  = sample_wor([1:P],Nvox,1);
                            y       = data_roi(:,voxind);
                    end
                    
                    % Exclude partition if necessary
                    if ~isempty(Pexclude)
                        idx = ismember(Partition, Pexclude);
                        Partition   = Partition(~idx);
                        Condition   = Condition(~idx);
                        y           = y(~idx,:);
                    end
                    
                    % Estimate cross-validated distance with noise estimate
                    Out = bmw1_calcDistVar_raw(y, [], Partition, Condition,...
                        'SVD',SVD,'demean',demean,'nDirL',nDirL,'nDirR',nDirR,...
                        'summary',0);
                    % [d,Sig]=distanceLDCweight(U,partition,conditionVec,X)
                    % store in var
                    dist_mva(:,i) = Out(1:ldist);
                    sigma_mva(:,i) = Out(ldist+1:lG+ldist);
                    
                    % Estimate cross-validated distance with noise estimate
                    % Univariate analysis ver.
                    Out = bmw1_calcDistVar_raw(mean(y,2), [], Partition, Condition,...
                        'SVD',SVD,'demean',demean,'nDirL',nDirL,'nDirR',nDirR,...
                        'summary',0);
                    % [d,Sig]=distanceLDCweight(U,partition,conditionVec,X)
                    % store in var
                    dist_uva(:,i) = Out(1:ldist);
                    sigma_uva(:,i) = Out(ldist+1:lG+ldist);
                end
                
                rdm.RDM(roi,:)      = nanmean(dist_mva,2);
                rdm.Sigma(roi,:)    = nanmean(sigma_mva,2);
                rdm.uRDM(roi,:)     = nanmean(dist_uva,2);
                rdm.uSigma(roi,:)   = nanmean(sigma_uva,2);
                rdm.name{roi,:}     = strjoin(roiname,'-');
                rdm.SN(roi,:)       = s;
                rdm.region(roi,:)   = Data.region(roi);
                rdm.nVox(roi,:)     = Nvox;
                for c = Conditions
                    rdm.meanAct(roi,c)  = mean(nanmean(data_roi(find(ismember(Condition,Conditions(c))),:)));
                end
                rdm.regType(roi,:)  = regType(Data.region(roi));
                rdm.regSide(roi,:)  = regSide(Data.region(roi));
            end
            
            RDMs = addstruct(RDMs,rdm);
            % save individual RDMs
            save(fullfile(regDir,subj,sprintf('%s_RDM_ROI_%s.mat',subj,ROI)),'-struct','RDMs');
            
            % show distance
            for roi = 1:numel(Data.name)
                subplot(5,ceil(numel(Data.name)/5),roi);
                histplot(RDMs.RDM(RDMs.region==roi,:)'); hold on
                title(RDMs.name{RDMs.region==roi})
            end
        end
        
        varargout = {RDMs};
    case 'ROI_get_RDMs_toshare_rsa'                 % calc RDM for each ROI and save in rsatoolbox format
        sn 		= varargin{1};
        fig     = 1;
        Nsimu   = 1;
        SVD     = 0;
        demean  = 0;
        Nvox    = 160;
        ROI     = 'all';
        Pexclude= []; % partition to be excluded
        nDirL = 8;
        nDirR = 8;
        options = {'Nsimu','SVD','ROI','demean','nDirL','nDirR','Pexclude'};
        vararginoptions({varargin{2:end}},options);
        
        nConditions = nDirL+nDirR+nDirL*nDirR;
        % record options
        for op=1:numel(options)
            Options.(options{op}) = eval(options{op});
        end
        Options.sn = sn;
        
        % get condition no.
        nBi = nDirL*nDirR;
        nConditons = nDirL+nDirR+nBi;
        Conditions = [1:nConditions];
        
        %% Run ROI analysis
        for s = sn
            RDMs = [];
            subj = subj_name{s};
            fprintf('subj = %s\n',subj)
            load(fullfile(regDir,subj,sprintf('%s_prewhitened_data_rsa_%s.mat',subj,ROI)))
            
            Partition   = Design.partition;
            Condition   = Design.condition;
            Con     = indicatorMatrix('allpairs',unique(Condition));
            
            
            % loop over rois
            for roi = 1:numel(Data.name)
                % get roi data
                roiname     = Data.name{roi};
                roiname     = strsplit(roiname,'_');
                roiname     = roiname(2:end)
                data_roi    = Data.Beta{roi};
                [N,P]       = size(data_roi);
                
                % loop over multiple iterations (random sub-space approach)
                for i = 1:Nsimu
                    switch (Nsimu)
                        case 1
                            y = data_roi;
                            Nvox = P;
                        otherwise
                            % get random sub-space of ROI (non-nan)
                            voxind  = sample_wor([1:P],Nvox,1);
                            y       = data_roi(:,voxind);
                    end
                    
                    % Exclude partition if necessary
                    if ~isempty(Pexclude)
                        idx = ismember(Partition, Pexclude);
                        Partition   = Partition(~idx);
                        Condition   = Condition(~idx);
                        y           = y(~idx,:);
                    end
                    
                    % Estimate cross-validated distance with noise estimate
                    %Out = sh1_calcDistVar_raw(y, [], Partition, Condition,...
                    %    'SVD',SVD,'demean',demean);
                    dist_mva(:,i) = rsa.distanceLDCweight(y,Partition',Condition',KX); % nVox already adjusted!
                    sigma_mva(:,i) = NaN; % not yet implemented
                    
                    % calc univariate distance
                    dist_uva(:,i) = rsa.distanceLDCweight(mean(y,2),Partition',Condition',KX); % nVox already adjusted!
                    sigma_uva(:,i) = NaN; % not yet implemented
                    
                    % store in var
                    %dist(:,i) = Out(1:ldist);
                    %sigma(:,i) = Out(ldist+1:lG+ldist);
                end
                
                rdm.RDM(roi,:)      = nanmean(dist_mva,2);
                rdm.Sigma(roi,:)    = nanmean(sigma_mva,2);
                rdm.uRDM(roi,:)     = nanmean(dist_uva,2);
                rdm.uSigma(roi,:)   = nanmean(sigma_uva,2);
                rdm.name{roi,:}     = strjoin(roiname,'-');
                rdm.SN(roi,:)       = s;
                rdm.region(roi,:)   = Data.region(roi);
                rdm.nVox(roi,:)     = Nvox;
                for c = Conditions
                    rdm.meanAct(roi,c)  = mean(nanmean(data_roi(find(ismember(Condition,Conditions(c))),:)));
                end
                rdm.regType(roi,:)  = regType(Data.region(roi));
                rdm.regSide(roi,:)  = regSide(Data.region(roi));
            end
            
            RDMs = addstruct(RDMs,rdm);
            % save individual RDMs
            save(fullfile(regDir,subj,sprintf('%s_RDM_ROI_rsa_%s.mat',subj,ROI)),'-struct','RDMs');
            
            % show distance
            if fig
                figure;
                for roi = 1:numel(Data.name)
                    subplot(5,ceil(numel(Data.name)/5),roi);
                    histplot(RDMs.RDM(RDMs.region==roi,:)'); hold on
                    title(RDMs.name{RDMs.region==roi})
                end
            end
        end
        
        varargout = {RDMs};
    case 'ROI_concat_RDMs'
        sn = varargin{1};
        ROI = 'all';
        RDMs = [];
        for s = sn
            subj = subj_name{s};
            rdms = load(fullfile(regDir,subj,sprintf('%s_RDM_ROI_%s.mat',subj,ROI)));
            
            RDMs = addstruct(RDMs,rdms);
        end
        % save
        save(fullfile(regDir,sprintf('allRDMs_ROI_%s.mat',ROI)),'-struct','RDMs');
        
        varargout = {RDMs};
    case 'ROI_distraw_effVox'                       % Extracts ROI data and calculates distances and effective number of voxels
        T=[];
        sn = varargin{1};
        glm = 8;
        ROI = 'all';
        fname 	= 'movement';       % using metric file which compares any movement against rest
        regions = 1:12;
        hemi = 1:2;
        vararginoptions({varargin{2:end}},{'glm','ROI','fname','regions','hemi'});
        
        Data=[];
        for s=sn
            fprintf('subj = %d\n',s)
            
            % load SPM.mat
            glmDirSubj=fullfile(baseDir,glmName{glm},subj_name{s});
            load(fullfile(glmDirSubj,'SPM.mat'));
            
            %SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
            
            % choose ROI
            switch (ROI)
                case 'func'
                    load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_reg800_%s_%d.mat',fname,glm)]));
                case 'all'
                    load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_regAll_%s.mat',fname)]));
            end
            
            % Get and tabulate design
            E   = load(fullfile(glmDirSubj,'SPM_info.mat'));
            indx= find(E.run<11);
            E   = getrow(E,indx);
            
            % Now get prewhitened betas from these regions
            NRegions = numel(regions)*numel(hemi);
            nRegions = numel(R);
            for reg = regions
                for hem = hemi
                    if reg<9
                        roi = numregions_surf*(hem-1)+reg;
                    elseif reg>=9
                        roi = 2*numregions_surf + numregions_BG*(hem-1)+(reg-8);
                    end
                    if (nRegions<NRegions)
                        warning('number of ROIs doesn''t much!');
                        break;
                    end
                    if (~isempty(R{roi}))
                        fprintf('extracting: %s ...',R{roi}.name);
                        
                        % get data
                        data = region_getdata(SPM.xY.VY,R{roi});
                        numVox = size(R{roi}.data,1);
                        % cut out NaN or all-zero voxels
                        idx = sum(data.*data,1)==0;
                        if sum(idx)==0
                            fprintf('... done (%d voxels)\n',numVox);
                        else
                            fprintf('... done (%d voxels, %d were discarded since containing no data)\n',numVox,sum(idx));
                        end
                        data = data(:,~idx);
                        numVox = size(data,2);
                        
                        % get the distances
                        [S.RDM,Sw,S.effVox,S.trSS]=rsa_distanceLDCsepPerm(data,SPM,E.movType);
                        S.Sigma  = Sw(:)';
                        S.region = reg;
                        S.subj  = s;
                        S.hemis = hem;
                        S.numVox = numVox;
                        
                        T=addstruct(T,S);
                    end;
                end;
            end;
        end;
        varargout={T};
        save(fullfile(regDir,sprintf('distances_sepPerm.glm%d.mat',glm)),'-struct','T');
    case 'ROI_distraw_changeBF'                     % Extracts ROI data and distances with adjusted BF
        T=[];
        sn = varargin{1};
        glm = 8;
        normmode = 'overall';
        ROI = 'all';
        fname 	= 'movement';       % using metric file which compares any movement against rest
        regions = 1:12;
        hemi = 1:2;
        %Phrf   = [3.389  16  1   1    6  1.0140    0.4990]'; % unit:seconds
        Phrf    = [3.389  16  1   1    6  2.8101    0.1727]';
        duration = 1;
        onsetshift = -5000/2720;
        vararginoptions(varargin(2:end),{'glm','ROI','fname','regions','hemi','Phrf'});
        
        for s=sn
            fprintf('subj = %d\n',s)
            
            % load SPM.mat
            glmDirSubj=fullfile(baseDir,glmName{glm},subj_name{s});
            load(fullfile(glmDirSubj,'SPM.mat'));
            %SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
            
            % Adjust hrf parameter and change design matrix
            %--------------------------------------------------------------
            % 1. clearn default onset and duration
            for r = 1:length(SPM.nscan)
                for u=1:length(SPM.Sess(r).U)
                    SPM.Sess(r).U(u).dur = ones(size(SPM.Sess(r).U(u).dur))*duration; % 1
                    SPM.Sess(r).U(u).ons = SPM.Sess(r).U(u).ons+onsetshift; % return to TR at announce trial
                end;
                SPM.Sess(r).U=spm_get_ons(SPM,r);
            end;
            % 2. change hrf and design matrix
            %figure(33);plot(SPM.xBF.bf);hold on;
            SPM.xBF.bf = spmj_hrf(SPM.xBF.dt,Phrf(1:7));
            %figure(33);plot(SPM.xBF.bf);hold on;
            SPM = spmj_fMRI_design_changeBF(SPM);
            fprintf('hrf parameter changed:\n');
            fprintf('P = [');fprintf('%d ',Phrf);fprintf(']\n');
            
            % choose functionaly/anatomical ROI
            switch (ROI)
                case 'func'
                    load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_reg800_%s_%d.mat',fname,glm)]));
                case 'all'
                    load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_regAll_%s.mat',fname)]));
            end
            
            % Get and tabulate design
            E   = load(fullfile(glmDirSubj,'SPM_info.mat'));
            indx= find(E.run<11);
            E   = getrow(E,indx);
            
            % Now get prewhitened betas from these regions
            NRegions = numel(regions)*numel(hemi);
            nRegions = numel(R);
            for reg = regions
                for hem = hemi
                    if reg<9
                        roi = numregions_surf*(hem-1)+reg;
                    elseif reg>=9
                        roi = 2*numregions_surf + numregions_BG*(hem-1)+(reg-8);
                    end
                    if (nRegions<NRegions)
                        warning('number of ROIs doesn''t much!');
                        break;
                    end
                    if (~isempty(R{roi}))
                        fprintf('extracting: %s ...',R{roi}.name);
                        
                        % get data
                        data = region_getdata(SPM.xY.VY,R{roi});
                        numVox = size(R{roi}.data,1);
                        % cut out NaN or all-zero voxels
                        idx = sum(data.*data,1)==0;
                        if sum(idx)==0
                            fprintf('... done (%d voxels)\n',numVox);
                        else
                            fprintf('... done (%d voxels, %d were discarded since containing no data)\n',numVox,sum(idx));
                        end
                        data = data(:,~idx);
                        numVox = size(data,2);
                        
                        % get the distances
                        [S.RDM,Sw,S.reliability,S.RDM_A,S.RDM_B] = rsa.spm.distanceLDCraw(data,SPM,E.movType,'normmode',normmode);
                        
                        S.Sigma  = Sw(:)';
                        S.region = reg;
                        S.subj  = s;
                        S.hemis = hem;
                        S.numVox = numVox;
                        
                        T=addstruct(T,S);
                    end;
                end;
            end;
        end;
        varargout={T};
        save(fullfile(regDir,sprintf('distances_adjhrf.glm%d.mat',glm)),'-struct','T');
    case 'ROI_Graw_changeBF'                        % Extracts ROI data and G with adjusted BF
        T=[];
        sn = varargin{1};
        glm = 8;
        normmode = 'overall';
        ROI = 'all';
        fname 	= 'movement';       % using metric file which compares any movement against rest
        regions = 1:12;
        hemi = 1:2;
        %Phrf   = [3.389  16  1   1    6  1.0140    0.4990]'; % unit:seconds
        Phrf    = [3.389  16  1   1    6  2.8101    0.1727]';
        duration = 1;
        onsetshift = -5000/2720;
        vararginoptions(varargin(2:end),{'glm','ROI','fname','regions','hemi','Phrf'});
        
        for s=sn
            fprintf('subj = %d\n',s)
            
            % load SPM.mat
            glmDirSubj=fullfile(baseDir,glmName{glm},subj_name{s});
            load(fullfile(glmDirSubj,'SPM.mat'));
            %SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
            
            % Adjust hrf parameter and change design matrix
            %--------------------------------------------------------------
            % 1. clearn default onset and duration
            for r = 1:length(SPM.nscan)
                for u=1:length(SPM.Sess(r).U)
                    SPM.Sess(r).U(u).dur = ones(size(SPM.Sess(r).U(u).dur))*duration; % 1
                    SPM.Sess(r).U(u).ons = SPM.Sess(r).U(u).ons+onsetshift; % return to TR at announce trial
                end;
                SPM.Sess(r).U=spm_get_ons(SPM,r);
            end;
            % 2. change hrf and design matrix
            %figure(33);plot(SPM.xBF.bf);hold on;
            SPM.xBF.bf = spmj_hrf(SPM.xBF.dt,Phrf(1:7));
            %figure(33);plot(SPM.xBF.bf);hold on;
            SPM = spmj_fMRI_design_changeBF(SPM);
            fprintf('hrf parameter changed:\n');
            fprintf('P = [');fprintf('%d ',Phrf);fprintf(']\n');
            
            % choose functionaly/anatomical ROI
            switch (ROI)
                case 'func'
                    load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_reg800_%s_%d.mat',fname,glm)]));
                case 'all'
                    load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_regAll_%s.mat',fname)]));
            end
            
            % Get and tabulate design
            E   = load(fullfile(glmDirSubj,'SPM_info.mat'));
            indx= find(E.run<11);
            E   = getrow(E,indx);
            
            % Now get prewhitened betas from these regions
            NRegions = numel(regions)*numel(hemi);
            nRegions = numel(R);
            for reg = regions
                for hem = hemi
                    if reg<9
                        roi = numregions_surf*(hem-1)+reg;
                    elseif reg>=9
                        roi = 2*numregions_surf + numregions_BG*(hem-1)+(reg-8);
                    end
                    if (nRegions<NRegions)
                        warning('number of ROIs doesn''t much!');
                        break;
                    end
                    if (~isempty(R{roi}))
                        fprintf('extracting: %s ...',R{roi}.name);
                        
                        % get data
                        data = region_getdata(SPM.xY.VY,R{roi});
                        numVox = size(R{roi}.data,1);
                        % cut out NaN or all-zero voxels
                        idx = sum(data.*data,1)==0;
                        if sum(idx)==0
                            fprintf('... done (%d voxels)\n',numVox);
                        else
                            fprintf('... done (%d voxels, %d were discarded since containing no data)\n',numVox,sum(idx));
                        end
                        data = data(:,~idx);
                        numVox = size(data,2);
                        
                        % get G
                        [~,Sw,G] = rsa.spm.distanceLDCraw(data,SPM,E.movType,'normmode',normmode,'G',1,'Sigma',1);
                        
                        S.G = G(:)';
                        S.Sigma  = Sw(:)';
                        S.region = reg;
                        S.subj  = s;
                        S.hemis = hem;
                        S.numVox = numVox;
                        
                        T=addstruct(T,S);
                    end;
                end;
            end;
        end;
        varargout={T};
        save(fullfile(regDir,sprintf('G_adjhrf.glm%d.mat',glm)),'-struct','T');
    case 'ROI_Graw_glm6direction'                         % Extracts ROI data and G with setting of glm 9-12
        sn = varargin{1};
        glms = 2;%fast no hpf
        normmode = 'runwise';
        ROI = 'all';
        fname 	= 'movement';       % using metric file which compares any movement against rest
        regions = 1:8;
        hemi = 1:2;
        
        vararginoptions(varargin(2:end),{'glm','ROI','fname','regions','hemi'});
        
        % make centering matrix
        Mu  = ones(nDirections,1)*ones(1,nDirections)/nDirections;
        Mb  = ones(nDirections^2,1)*ones(1,nDirections^2)/nDirections^2;
        M   = blockdiag(Mu,Mu,Mb);
        C   = eye(nConditions)-M;
        
        for glm=glms
            T = [];
            for s=sn
                fprintf('subj = %d\n',s)
                
                % load SPM.mat
                glmDirSubj=fullfile(baseDir,glmName{glm},subj_name{s});
                load(fullfile(glmDirSubj,'SPM.mat'));
                %SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
                
                % choose functionaly/anatomical ROI
                switch (ROI)
                    case 'func'
                        load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_reg800_%s_%d.mat',fname,glm)]));
                    case 'all'
                        load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_regAll_%s.mat',fname)]));
                end
                
                % Get and tabulate design
                E   = load(fullfile(glmDirSubj,'SPM_info.mat'));
                indx= find(E.run<11);
                E   = getrow(E,indx);
                
                % Now get prewhitened betas from these regions
                NRegions = numel(regions)*numel(hemi);
                nRegions = numel(R);
                for reg = regions
                    for hem = hemi
                        if reg<9
                            roi = numregions_surf*(hem-1)+reg;
                        elseif reg>=9
                            roi = 2*numregions_surf + numregions_BG*(hem-1)+(reg-8);
                        end
                        if (nRegions<NRegions)
                            warning('number of ROIs doesn''t much!');
                            break;
                        end
                        if (~isempty(R{roi}))
                            fprintf('extracting: %s ...',R{roi}.name);
                            
                            % get data
                            data = region_getdata(SPM.xY.VY,R{roi});
                            numVox = size(R{roi}.data,1);
                            % cut out NaN or all-zero voxels
                            idx = sum(data.*data,1)==0;
                            if sum(idx)==0
                                fprintf('... done (%d voxels)\n',numVox);
                            else
                                fprintf('... done (%d voxels, %d were discarded since containing no data)\n',numVox,sum(idx));
                            end
                            data = data(:,~idx);
                            numVox = size(data,2);
                            
                            % get G
                            [~,Sw] = rsa.spm.distanceLDCraw(data,SPM,E.movType,...
                                'normmode',normmode);
                            G = rsa.spm.crossvalIPMraw(data,SPM,E.movType,...
                                'normmode',normmode);
                            %[Beta] = rsa.spm.noiseNormalizeBeta(data,SPM,'normmode',normmode);
                            %[G,Sw] = crossval_estG(Beta(1:nConditions,:),indicatorMatrix('identity_p',1:nConditions),E.run);
                            
                            % decompose G
                            Ga = C*G*C'; % pattern around mean
                            Gb = M*G*M'; % mean pattern
                            Gab= C*G*M';
                            Gba= M*G*C';
                            
                            S.G         = G(:)';
                            S.Ga        = Ga(:)';
                            S.Gb        = Gb(:)';
                            S.Gab       = Gab(:)';
                            S.Gba       = Gba(:)';
                            S.Sigma     = Sw(:)';
                            S.region    = reg;
                            S.subj      = s;
                            S.hemis     = hem;
                            S.numVox    = numVox;
                            S.glm       = glm;
                            
                            T = addstruct(T,S);
                        end;
                    end;
                end;
            end;
            varargout = {T};
            save(fullfile(regDir,sprintf('G_raw.glm%d.mat',glm)),'-struct','T');
        end
    case 'ROI_sepGraw_glm9-12'                      % Extracts ROI data and calculate IPM and RDM for both pattern and mean patterns
        sn      = varargin{1};
        glms    = 10;
        fig     = 0;
        normmode = 'runwise';
        ROI     = 'all';
        fname 	= 'movement';       % using metric file which compares any movement against rest
        regions = 1:12;
        hemi    = 1:2;
        
        vararginoptions(varargin(2:end),{'glms','ROI','fname','regions','hemi','fig'});
        
        for glm=glms
            T = [];
            for s=sn
                fprintf('subj = %d\n',s)
                
                % load SPM.mat
                glmDirSubj=fullfile(baseDir,glmName{glm},subj_name{s});
                load(fullfile(glmDirSubj,'SPM.mat'));
                %SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
                % get 1st-level design matrix for optimal extimatin of U
                X = SPM.xX.xKXs.X;
                
                % choose functionaly/anatomical ROI
                switch (ROI)
                    case 'func'
                        load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_reg800_%s_%d.mat',fname,glm)]));
                    case 'all'
                        load(fullfile(regDir,subj_name{s},[subj_name{s} sprintf('_regAll_%s.mat',fname)]));
                end
                
                % Get and tabulate design
                E   = load(fullfile(glmDirSubj,'SPM_info.mat'));
                partition   = E.run;
                partitions  = unique(E.run);
                partition   = [partition;partitions]; % add intercept
                Npartition  = length(partitions);
                Z = indicatorMatrix('identity_p',repmat(1:(nConditions),1,10)');
                Z = blockdiag(Z,eye(Npartition)); % add intercept
                C = indicatorMatrix('allpairs',[1:nConditions*2]);
                
                % Now get prewhitened betas from these regions
                NRegions = numel(regions)*numel(hemi);
                nRegions = numel(R);
                for reg = regions
                    for hem = hemi
                        if reg<9
                            roi = numregions_surf*(hem-1)+reg;
                        elseif reg>=9
                            roi = 2*numregions_surf + numregions_BG*(hem-1)+(reg-8);
                        end
                        if (nRegions<NRegions)
                            warning('number of ROIs doesn''t much!');
                            break;
                        end
                        if (~isempty(R{roi}))
                            fprintf('extracting: %s ...',R{roi}.name);
                            
                            % get data
                            data = region_getdata(SPM.xY.VY,R{roi});
                            numVox = size(R{roi}.data,1);
                            % cut out NaN or all-zero voxels
                            idx = sum(data.*data,1)==0;
                            if sum(idx)==0
                                fprintf('... done (%d voxels)\n',numVox);
                            else
                                fprintf('... done (%d voxels, %d were discarded since containing no data)\n',numVox,sum(idx));
                            end
                            data = data(:,~idx);
                            numVox = size(data,2);
                            
                            % get prewhitened beta
                            [Beta] = rsa.spm.noiseNormalizeBeta(data,SPM,'normmode',normmode);
                            
                            % estimate patterns
                            Ua = zeros(nConditions,numVox,Npartition);
                            Ub = Ua;
                            for r = 1:Npartition
                                idx     = partition==partitions(r);
                                rBeta   = Beta(idx,:);
                                rZ      = Z(idx,any(Z,1));
                                rZ      = rZ(:,any(rZ,1));
                                rX      = X(:,idx);
                                U       = pinv(rX*rZ)*rX*rBeta; % optimal extimate taking 1st-level design matrix into account
                                U       = real(U(1:nConditions,:)); % discard resting baseline (run wise intercept)
                                
                                % subtract mean per each movement categoly
                                % (uni-left, uni-right, and bimanual)
                                Ua(1:6,:,r)     = bsxfun(@minus,U(1:6,:),nanmean(U(1:6,:),1));
                                Ua(7:12,:,r)    = bsxfun(@minus,U(7:12,:),nanmean(U(7:12,:),1));
                                Ua(13:nConditions,:,r)   = bsxfun(@minus,U(13:nConditions,:),nanmean(U(13:nConditions,:),1));
                                
                                % get mean pattern
                                Ub(:,:,r) = U-Ua(:,:,r);
                            end
                            
                            % get crossvalidated inner product for both Ua
                            % and Ub (U=[Ua;Ub])
                            Npair = Npartition*(Npartition-1);
                            IPM = zeros(2*nConditions,2*nConditions,Npair);
                            RDM = zeros(2*nConditions,2*nConditions);
                            count = 1;
                            for i=1:Npartition
                                A = [Ua(:,:,i);Ub(:,:,i)];
                                for j=1:Npartition
                                    if i~=j
                                        B = [Ua(:,:,j);Ub(:,:,j)];
                                        IPM(:,:,count) = A*B';
                                        count = count+1;
                                    end
                                end
                            end
                            IPM = sum(IPM,3)/numVox/Npair;
                            IPM = 0.5*(IPM+IPM');
                            RDM = squareform(diag(C*IPM*C'));
                            
                            if fig
                                subplot(2,1,1)
                                imagesc(IPM);colorbar();title('IPM')
                                subplot(2,1,2)
                                imagesc(RDM);colorbar();title('RDM')
                            end
                            
                            S.IPM = rsa_vectorizeIPM(IPM);
                            S.RDM = rsa_vectorizeRDM(RDM);
                            S.region    = reg;
                            S.subj      = s;
                            S.hemis     = hem;
                            S.numVox    = numVox;
                            S.glm       = glm;
                            
                            T = addstruct(T,S);
                        end;
                    end;
                end;
            end;
            varargout = {T};
            save(fullfile(regDir,sprintf('IPMRDMsep.glm%d.mat',glm)),'-struct','T');
        end
    case 'ROI_compare_method' % Check if 'ROI_sepGraw_glm9-12' and 'ROI_Graw_glm9-12' gives the same value
        glm = 10;
        
        % load data
        T1 = load(fullfile(regDir,sprintf('G_raw.glm%d.mat',glm)));
        T2 = load(fullfile(regDir,sprintf('IPMRDMsep.glm%d.mat',glm)));
        
        for hemis = [1,2];
            for reg=1:8;
                t1 = getrow(T1,T1.hemis==hemis&T1.region==reg);
                t2 = getrow(T2,T2.hemis==hemis&T2.region==reg);
                
                % get IPMs
                IPMa    = nanmean(t1.Ga,1);
                IPMb    = nanmean(t1.Gb,1);
                IPMab   = nanmean(t1.Gab,1);
                IPMba   = nanmean(t1.Gba,1);
                IPM2    = nanmean(t2.IPM,1);
                
                Ga  = reshape(IPMa,48,48);
                Gb  = reshape(IPMb,48,48);
                Gab = reshape(IPMab,48,48);
                Gba = reshape(IPMba,48,48);
                
                sum(sum((t1.G-t1.Ga-t1.Gb-t1.Gab-t1.Gba).*(t1.G-t1.Ga-t1.Gb-t1.Gab-t1.Gba)))
                
                G1  = [Ga,Gab;Gba,Gb];
                G2  = rsa_squareIPM(IPM2);
                
                UL1 = G1(1:6,1:6);
                UL2 = G2(1:6,1:6);
                
                D   = G2 - G1;
                
                corrcoef(UL1(:),UL2(:))
                corrcoef(G1(:),G2(:))
                scale = sum(UL2(:))/sum(UL1(:))
                SSR = sum(sum(D.*D))
            end
        end
        
    case 'ROI_show_RDMs'
        type = varargin{1}; % 'individual','average','both'
        glm = 8;
        regions = [1:8];
        hemi = [2];
        subspace = 0;
        cond = 'All';
        squareTransform = 1;
        vararginoptions(varargin(2:end),{'glm','regions','hemi','subspace','cond'})
        
        % load data set
        if subspace==0
            %             T = load(fullfile(regDir,sprintf('distances_sepPerm.glm%d.mat',glm)));
            T = load(fullfile(regDir,sprintf('distances_adjhrf.glm%d.mat',glm)));
        elseif  subspace==1
            T = load(fullfile(regDir,sprintf('distances_sepPerm.glm%d_subspace.mat',glm)));
            T = load(fullfile(regDir,sprintf('distances_sepPerm.glm%d.iter1000_subspace.mat',glm)));
        end
        if squareTransform
            T.RDM = ssqrt(T.RDM);
        end
        
        rdm = [];
        for reg = regions;%unique(T.region)'
            for h = hemi
                subset  = T.region==reg&T.hemis==h;
                D       = getrow(T,subset);
                
                switch (type)
                    case 'average'
                        R.RDM = nanmean(D.RDM,1);
                        R.region = reg;
                        R.hemis = h;
                        R.name = {sprintf('%s.%s-average',regname{reg},hemName{h})};%{''};%{sprintf('%s_%s_set%s',regname{reg},hemName_s{h},D.chunkset(1))};
                        rdm = addstruct(rdm,R);
                    case 'individual'
                        D.name = repmat({sprintf('%s.%s',regname{reg},hemName{h})},length(D.subj),1);%{''};%{sprintf('%s_%s_set%s',regname{reg},hemName_s{h},D.chunkset(1))};
                        rdm = addstruct(rdm,D);
                    case 'both'
                        % average
                        R.RDM = nanmean(D.RDM,1);
                        R.region = reg;
                        R.hemis = h;
                        R.name = {sprintf('%s.%s-average',regname{reg},hemName{h})};%{''};%{sprintf('%s_%s_set%s',regname{reg},hemName_s{h},D.chunkset(1))};
                        rdm = addstruct(rdm,R);
                        
                        % individual
                        D.name = repmat({sprintf('%s.%s',regname{reg},hemName{h})},length(D.subj),1);%{''};%{sprintf('%s_%s_set%s',regname{reg},hemName_s{h},D.chunkset(1))};
                        rdm = addstruct(rdm,D);
                end
                
            end
        end
        RDM = rdm;
        RDM.RDM = splitUB(RDM.RDM,cond);
        %figure('position',[50 50 900/6 900]);
        rsa.fig.imageRDMs(RDM,'aspect',1/5,'clims',[-1 1]);
        %end
    case 'ROI_reliability'                          % Calc within- and between-subject reliability of RDM
        glm = 8;
        sep = 1; % if we calculate reliability separately for unimanual and bimanual
        vararginoptions(varargin(:),{'sep'});
        
        % load data
        T = load(fullfile(regDir,sprintf('distances_adjhrf.glm%d.mat',glm)));
        
        regions = unique(T.region);
        hemis = unique(T.hemis);
        cases = {'Uni','UL','UR','Bi'};
        Reliability = [];
        for h=hemis'
            for reg=regions'
                index = find(T.hemis==h&T.region==reg);
                rdm = T.RDM(index,:);
                rdma = T.RDM_A(index,:);
                rdmb = T.RDM_B(index,:);
                
                if sep == 0
                    % within
                    Rwithin = mean(T.reliability(index,1));
                    SDwithin = std(T.reliability(index,1));
                    T.reliability_within = T.reliability;
                    
                    % between
                    Ctmp = corr(rdm');
                    Rbtw = mean(Ctmp(tril(~logical(eye(size(Ctmp))))));
                    SDbtw = std(Ctmp(tril(~logical(eye(size(Ctmp))))));
                    T.reliability_btw(index,:) = Ctmp;
                else
                    for c=1:numel(cases)
                        % split RDM into relevant parts
                        rdma_ = splitUB(rdma,cases{c});
                        rdmb_ = splitUB(rdmb,cases{c});
                        rdm_ = splitUB(rdm,cases{c});
                        
                        % within
                        Rwithin(1,c) = mean(diag(corr(rdma_',rdmb_')));
                        SDwithin(1,c) = std(diag(corr(rdma_',rdmb_')));
                        T.reliability_within(index,c,:) = diag(corr(rdma_',rdmb_'));
                        
                        % between
                        Ctmp = corr(rdm_');
                        Rbtw(1,c) = mean(Ctmp(tril(~logical(eye(size(Ctmp))))));
                        SDbtw(1,c) = std(Ctmp(tril(~logical(eye(size(Ctmp))))));
                        T.reliability_btw(index,c,:) = Ctmp;
                        
                        % same between-subject correlation using 0.5(rdma+rdmb)
                        rdm_split = 0.5*(rdma_+rdmb_);
                        Ctmp2 = corr(rdm_split');
                        Rbtw2(1,c) = mean(Ctmp2(tril(~logical(eye(size(Ctmp2))))));
                        SDbtw2(1,c) = std(Ctmp2(tril(~logical(eye(size(Ctmp2))))));
                        T.reliability_btw2(index,c,:) = Ctmp2;
                        
                        % correlation between normally calculated RDM and
                        % mean of split halves
                        Cspno = corr(rdm_',rdm_split');
                        Rspno(1,c) = mean(diag(Cspno));
                        T.reliability_split_normal(index,c,:) = diag(Cspno);
                    end
                end
                Rel.Rwithin     = Rwithin;
                Rel.Rbtw        = Rbtw;
                Rel.Rbtw2       = Rbtw2;
                Rel.Rspno       = Rspno; % split vs normal
                Rel.SDwithin    = SDwithin;
                Rel.SDbtw       = SDbtw;
                Rel.SDbtw2      = SDbtw2;
                Rel.region      = reg;
                Rel.hemis       = h;
                Reliability     = addstruct(Reliability,Rel);
            end
        end
        % display summary
        disp('Within-subject reliability for each ROI');
        disp('ROIs      | Reliability');
        pivottablerow([Reliability.hemis,Reliability.region],...
            Reliability.Rwithin,'(x)');
        
        disp('Between-subject reliability for each ROI');
        disp('ROIs      | Reliability');
        pivottablerow([Reliability.hemis,Reliability.region],...
            Reliability.Rbtw,'(x)');
        
        disp('Between-subject reliability for each ROI (RDMA+RDMB)/2');
        disp('ROIs      | Reliability');
        pivottablerow([Reliability.hemis,Reliability.region],...
            Reliability.Rbtw2,'(x)');
        
        varargout = {Reliability,T};
    case 'ROI_plot'
        sn      = varargin{1};
        ROIdef  = varargin{2};
        reg     = [1 2 3 4 5 7 8];
        nDirL   = 8;
        nDirR   = 8;
        mva     = 1;
        vararginoptions(varargin(3:end),{'reg','nDirL','nDirR','mva'})
        
        for s = sn
            for r = 1:numel(ROIdef)
                ROI = ROIdef{r};
                bmw1_imana('ROI_show_RDMs',s,ROI,'nDirL',nDirL,'nDirR',nDirR,'reg',reg,'mva',mva);
                bmw1_imana('ROI_get_summary_RDMs',s,ROI,'nDirL',nDirL,'nDirR',nDirR,'reg',reg,'mva',mva);
                pause();
            end
        end
    case 'ROI_show_individRDMs'                     % Show RDM for each ROI
        sn      = varargin{1};
        ROI     = varargin{2};
        reg     = [1 2 3 4 5 7 8];
        nDirL   = 8;
        nDirR   = 8;
        mva     = 1;
        vararginoptions(varargin(3:end),{'reg','nDirL','nDirR','mva'})
        nFig = 3;
        
        nBi = nDirL*nDirR;
        nAll = nDirL+nDirR+nBi;
        
        % load RDM
        RDMs = load(fullfile(regDir,subj_name{sn},sprintf('%s_RDM_ROI_%s.mat',subj_name{sn},ROI)));
        
        % mva option
        if ~mva
            if isfield(RDMs,'uRDM')
                RDMs.RDM    = RDMs.uRDM;
                RDMs.Sigma  = RDMs.uSigma;
            else
                warning('Missing uRDM!')
                return;
            end
        end
        % normalise RDM by nVox if necessary
        if isempty(strfind(ROI,'rsa'))
            RDMs.RDM = RDMs.RDM./repmat(RDMs.nVox,1,size(RDMs.RDM,2));
        end
        
        % plot RDMs
        offset = gcf;
        for hemi = 1:2
            figure(offset+hemi);
            
            tmp_rdm = getrow(RDMs,RDMs.regSide==hemi&ismember(RDMs.regType,reg));
            showRDMs(rsa_foldRDMs(tmp_rdm),hemi,0); hold on
            
            for i=1:numel(reg);
                subplot(nFig, nFig, i);
                drawline([nDirL+0.5, (nDirL+nDirR)+0.5],'dir','horiz','color',[1 1 1]);
                drawline([nDirL+0.5, (nDirL+nDirR)+0.5],'dir','vert','color',[1 1 1]);
                drawline([(nDirL+nDirR+0.5):nDirR:(nAll+0.5-nDirR)],'dir','vert','lim',[(nDirL+nDirR)+0.5, nAll+0.5],'color',[1 1 1]);
                drawline([(nDirL+nDirR+0.5):nDirR:(nAll+0.5-nDirR)],'dir','horiz','lim',[(nDirL+nDirR)+0.5, nAll+0.5],'color',[1 1 1]);
                set(gca,'fontsize',12);
            end
            set(gcf,'name',sprintf('RDMs_%s_%s',subj_name{sn},ROI),...
                'position',[50 50 800 800],'color','w');
        end
        
        % plot histgram
        figure('name','Histgram of distances');
        for roi = 1:numel(unique(RDMs.region))
            subplot(5,ceil(numel(unique(RDMs.name))/5),roi);
            histplot(RDMs.RDM(RDMs.region==roi,:)'); hold on
            title(RDMs.name{RDMs.region==roi})
        end
    case 'ROI_get_summary_RDMs'                     %
        sn      = varargin{1};
        ROI     = varargin{2};
        reg     = [1 2 3 4 5 7 8];
        nDirL   = 8;
        nDirR   = 8;
        mva     = 1;
        vararginoptions(varargin(3:end),{'reg','nDirL','nDirR','mva'})
        
        regnam    = {'S1','M1','PMd','PMv','SMA','SPLa','SPLp'};
        nConditions = nDirL+nDirR+nDirL*nDirR;
        nBi = nDirL*nDirR;
        if nDirL<8
            ylim = [-0.05 0.1];
        else
            ylim = [-0.1 0.2];
        end
        
        % load RDM
        RDMs = load(fullfile(regDir,subj_name{sn},sprintf('%s_RDM_ROI_%s.mat',subj_name{sn},ROI)));
        
        % mva option
        if ~mva
            if isfield(RDMs,'uRDM')
                RDMs.RDM    = RDMs.uRDM;
                RDMs.Sigma  = RDMs.uSigma;
            else
                warning('Missing uRDM!')
                return;
            end
        end
        % normalise RDM by nVox if necessary
        if isempty(strfind(ROI,'rsa'))
            RDMs.RDM = RDMs.RDM./repmat(RDMs.nVox,1,size(RDMs.RDM,2));
        end
        
        % get summary distance
        S           = getrow(RDMs,ismember(RDMs.regType,reg));
        [Mean, SD]  = bmw1_processRDM('summary',S.RDM,nDirL,nDirR);
        f           = fieldnames(Mean);
        for i=1:numel(f)
            S.(['mean_',f{i}]) = Mean.(f{i});
            S.(['sd_',f{i}]) = SD.(f{i});
        end
        labels = S.name;
        
        % plot results
        ys      = {'mean_ul','mean_cross_u','mean_ur','mean_cross_bl','mean_bi','mean_cross_br'};
        %errs    = {'se_ul','se_cross_u','se_ur','se_cross_bl','se_bi','se_cross_br'};
        errs    = {'sd_ul','sd_cross_u','sd_ur','sd_cross_bl','sd_bi','sd_cross_br'};
        titles  = {'Unimanual left','Uni L vs Uni R','Unimanual right','Uni L vs Bi','Bimanual','Uni R vs Bi'};
        figure('name',sprintf('SummaryRDMs_%s_%s',subj_name{sn},ROI));
        for type = 1:6
            y   = S.(ys{type});
            err = S.(errs{type});
            subplot(3,2,type);
            xloc = barplot([S.regSide,S.regType],y,'errorfcn',err',...
                'barwidth',0.8,'style_bold');
            title(titles{type})
            xlabel({'LH <--  ROIs --> RH'});
            ylabel('mean distance')
            set(gca,'Xticklabel',labels);
            %         set(gca,'ylim',ylim);
            
            % show summary
            %fprintf('\n|--- mean(distance(%s)) ---| \n',titles{type})
            %pivottable(S.regSide,S.regType,y,'nanmean');
            %fprintf('\n|--- se(distance(%s)) ---| \n',titles{type})
            %pivottable(S.regSide,S.regType,err,'nanmean');
        end
        varargout = {S};
        
    case 'ROI_MDS'                                  % Show MDS plot with different contrasts
        sn  = varargin{1};
        reg = varargin{3};
        hem = varargin{4};
        glm = varargin{2};
        
        vararginoptions(varargin(5:end),{'ROI'});
        
        % Load data
        T = load(fullfile(regDir,sprintf('G_raw.glm%d.mat',glm)));
        T = getrow(T,ismember(T.subj,sn));
        
        % Make contrasts (thetaL,thetaR,UL/UR/Bi)
        thetaL  = [kron([1 0]',[1:6]');zeros(36,1)];
        thetaR  = [kron([0 1]',[1:6]');zeros(36,1)];
        comLR   = [zeros(12,1);[1:36]'];
        thetaLB = [zeros(12,1);kron([1:6]',ones(6,1))];
        thetaRB = [zeros(12,1);kron(ones(6,1),[1:6]')];
        movType = [ones(6,1);2*ones(6,1);3*ones(36,1)];
        
        CL = indicatorMatrix('identity_p',thetaL);
        CR = indicatorMatrix('identity_p',thetaR);
        Cm = indicatorMatrix('identity',movType);
        Cb = indicatorMatrix('identity_p',comLR);
        CLb= indicatorMatrix('identity_p',thetaLB);
        CRb= indicatorMatrix('identity_p',thetaRB);
        CL = bsxfun(@minus,CL,mean(CL,2));
        CR = bsxfun(@minus,CR,mean(CR,2));
        Cm = bsxfun(@minus,Cm,mean(Cm,2));
        Cb = bsxfun(@minus,Cb,mean(Cb,2));
        CLb = bsxfun(@minus,CLb,mean(CLb,2));
        CRb = bsxfun(@minus,CRb,mean(CRb,2));
        Ca = eye(nConditions)-ones(nConditions)/nConditions;
        
        % label
        label = kron([60;60],[1:nDirections]');
        for i=1:nDirections
            for j=1:nDirections
                label = [label; str2double(sprintf('%d%d',60*i,60*j))];
            end
        end
        
        % Get MDS per roi
        for h=hem
            for roi=reg
                D = getrow(T,T.region==roi&T.hemis==h);
                
                IPM = mean(D.G,1);
                IPM = rsa_vectorizeIPM(reshape(IPM,nConditions,nConditions));
                
                Y{1} = rsa_classicalMDS(IPM,'mode','IPM'); % data driven
                Y{2} = rsa_classicalMDS(IPM,'mode','IPM','contrast',Ca); % all conditions
                Y{3} = rsa_classicalMDS(IPM,'mode','IPM','contrast',CL); % left wrist direction (uni)
                Y{4} = rsa_classicalMDS(IPM,'mode','IPM','contrast',CR); % right wrist direction (uni)
                Y{5} = rsa_classicalMDS(IPM,'mode','IPM','contrast',Cb); % all bimanual combinations
                Y{6} = rsa_classicalMDS(IPM,'mode','IPM','contrast',Cm); % uni-left/uni-right/bi
                Y{7} = rsa_classicalMDS(IPM,'mode','IPM','contrast',CLb); % left wrist direction (bi)
                Y{8} = rsa_classicalMDS(IPM,'mode','IPM','contrast',CRb); % right wrist direction (bi)
                
                
                %clf;
                % h=axes('position',[0 0 1 1],'Visible','off');
                figure('Name',sprintf('%s-%s (glm=%s)',regname{roi},hemName{h},glmName{glm}),'NumberTitle','off')
                
                subplot(2,3,1);
                scatterplotMDS(Y{1}(:,1:3),movType,label);
                title('raw')
                
                %                 subplot(2,3,2);
                %                 scatterplotMDS(Y{2}(:,1:3),movType,label);
                %                 title('all conditions')
                
                subplot(2,3,2);
                scatterplotMDS(Y{3}(:,1:3),movType,label);
                title('left directions (uni)')
                
                subplot(2,3,3);
                scatterplotMDS(Y{4}(:,1:3),movType,label);
                title('right directions (uni)')
                
                %                 subplot(2,3,5);
                %                 scatterplotMDS(Y{5}(:,1:3),movType,label);
                %                 title('all bimanual combinations')
                
                subplot(2,3,4);
                scatterplotMDS(Y{6}(:,1:3),movType,label);
                title('uni-left/uni-right/bi')
                
                subplot(2,3,5);
                scatterplotMDS(Y{7}(:,1:3),movType,label);
                title('left directions (bi)')
                
                subplot(2,3,6);
                scatterplotMDS(Y{8}(:,1:3),movType,label);
                title('right directions (bi)')
            end
        end
        
    case 'ROI_show_IPM_RDM'
        
        glm = 10;
        hemi = 1:2;
        regions = [2 3 6 8];
        data = 'IPM';
        
        vararginoptions(varargin,{'hemi','regions','data'})
        
        % load data
        T = load(fullfile(regDir,sprintf('IPMRDMsep.glm%d.mat',glm)));
        
        D = tapply(T,{'hemis','region'},{'IPM','nanmean(x,1)','name','IPM'},...
            {'RDM','nanmean(x,1)','name','RDM'});
        
        % plot
        nhem = length(hemi);
        nreg = length(regions);
        figure('name',data);
        rc=0;
        for reg=regions
            rc=rc+1;
            for h=hemi
                D_ = getrow(D,D.hemis==h&D.region==reg);
                
                tmp = rsa_squareIPM(D_.IPM);
                dtmp = diag(tmp);
                for i=1:96
                    for j=1:96
                        tmp(i,j) = tmp(i,j)/ssqrt(dtmp(i)*dtmp(j));
                    end
                end
                
                subplot(nreg,nhem,nhem*(rc-1)+h)
                map = eval(sprintf('rsa_square%s(D_.%s)',data,data));
                
                imagesc(map);colorbar();axis square
                caxis([-0.01 0.08]);
                title(sprintf('%s-%s',hemName{h},regname{reg}));
            end
        end
        
    case 'ROI_RDM_meandist'                         % Show mean distance for within-condition RDMs
        hemi = [1:2];
        regions = [1:8];
        glm = 10;
        sqrtTransform = 1;
        ex180=1;
        subjects = [1:7];
        plotfcn = 'barplot';
        
        vararginoptions(varargin,{'hemi','regions','sqrtTransform','ex180','subjects','plotfcn','glm'});
        
        % load data
        T = load(fullfile(regDir,sprintf('IPMRDMsep.glm%d.mat',glm)));
        T = getrow(T,ismember(T.subj,subjects));
        R = [];
        
        % exclude pair with 180 apart movement
        C = true(nDirections);
        for i=1:nDirections
            C(i,i) = false;
        end
        Cb_ = true(36);
        for j=1:36
            Cb_(j,j) = false;
        end
        C_  = C;
        C   = circshift(C,[0,3]);%,2);
        Cl  = kron(C,true(nDirections));
        Cr  = kron(true(nDirections),C);
        C   = C&C_;
        Cb  = Cl&Cr&Cb_;
        
        if ex180==0
            C = C_;
            Cb= Cb_;
        end
        
        % make index for bimanual contra and ipsi distance
        Cc = logical(eye(nDirections));
        Ccr= kron(Cc,true(nDirections))&Cb; % was Cb_
        Ccl= kron(true(nDirections),Cc)&Cb; % was Cb_
        
        % another way of calculating bimanual contra&ipsi
        Zcontra = kron(eye(nDirections),ones(nDirections,1));
        Zipsi = kron(ones(nDirections,1),eye(nDirections));
        Z = [zeros(2*nDirections),[Zcontra';Zipsi']];
        Z = bsxfun(@rdivide,Z,sum(Z,2));
        Con = indicatorMatrix('allpairs',[1:2*nDirections]);
        
        % loop over subjects,hemispheres,rois
        %subjects = unique(T.subj);
        for s=subjects
            for h=hemi
                for reg = regions
                    T_ = getrow(T,T.hemis==h&T.region==reg&T.subj==s);
                    
                    RDM = rsa_squareRDM(T_.RDM);
                    IPM = rsa_squareIPM(T_.IPM);
                    IPM = IPM(1:nConditions,1:nConditions);
                    
                    % pick corresponding parts of RDM
                    UL = RDM(1:nDirections,1:nDirections);
                    UR = RDM(nDirections+1:2*nDirections,nDirections+1:2*nDirections);
                    Bi = RDM(2*nDirections+1:nConditions,2*nDirections+1:nConditions);
                    
                    % another way of calculating bimanual contra&ipsi
                    Bi2 = rsa_squareRDM(diag(Con*Z*IPM*Z'*Con')');
                    Bi2L = Bi2(1:nDirections,1:nDirections);
                    Bi2R = Bi2(nDirections+1:2*nDirections,nDirections+1:2*nDirections);
                    
                    % ssqrt
                    if sqrtTransform
                        UL = ssqrt(UL);
                        UR = ssqrt(UR);
                        Bi = ssqrt(Bi);
                        Bi2L = ssqrt(Bi2L);
                        Bi2R = ssqrt(Bi2R);
                        ylim = [-0.01 0.08];
                        ylim1 = [-0.01 0.3];
                    else
                        ylim = [-0.001 0.07];
                        ylim1 = [-0.001 0.09];
                    end
                    
                    % overall distance
                    S.UL = nanmean(UL(C));
                    S.UR = nanmean(UR(C));
                    S.Bi = nanmean(Bi(Cb));
                    S.BiL = nanmean(Bi(Ccl));
                    S.BiR = nanmean(Bi(Ccr));
                    
                    S.Bi2L = nanmean(Bi2L(C));
                    S.Bi2R = nanmean(Bi2R(C));
                    
                    S.subj = s;
                    S.hemis = h;
                    S.region = reg;
                    
                    R = addstruct(R,S);
                end
            end
        end
        
        % plot
        fillcolor   = {purple,gold,darkpurple,darkgold};
        variables   = {'UL','UR','Bi2L','Bi2R'};
        labels      = {'UniL','UniR','BiL','BiR'};
        bmw1_imana('ROI_RDM_plot_testzero2',R,variables,...
            'fillcolor',fillcolor,'labels',labels,...
            'yrange',[-0.05 0.3]);
        
        %         variables = {'BiL','BiR','Bi2L','Bi2R'};
        %         bmw1_imana('ROI_RDM_plot_testzero',R,variables,plotfcn,'Bimanual ipsi- vs contra-');
        
        bmw1_imana('ROI_RDM_plot_conip',R);
        
        varargout = {R};
    case 'ROI_RDM_contrast'                         % Compare extrinsic, intrinsic, unrelated pairs of movements
        hemi = [1:2];
        regions = [1:8];
        glm = 10;
        sqrtTransform = 1;
        outback = 1;
        subjects = [1:7];
        
        % todo separate out-out/back-back and out-back/back-out
        
        vararginoptions(varargin,{'hemi','regions','sqrtTransform','outback','subjects'});
        
        % load data
        T = load(fullfile(regDir,sprintf('IPMRDMsep.glm%d.mat',glm)));
        R = [];
        
        % define contrast for unimanual movements
        Coverall = ones(nDirections);
        C = Coverall-1;
        extrinsic = [1:nDirections];
        intrinsic = [4 3 2 1 6 5];
        for i=1:nDirections
            for j=1:nDirections
                if j==extrinsic(i)
                    CU(i,j) = 1;
                end
                if j==intrinsic(i)
                    CU(i,j) = 2;
                end
            end
        end
        
        % deal with back phase of movement?
        CUback = circshift(CU,[0,3]); % = circshift(CU,3,2)
        CUback(CUback==1) = 3;
        CUback(CUback==2) = 6;
        if outback==1;CU = CU+CUback;end
        
        % extend into bimanual movements
        Coverall_Bi = kron(Coverall,ones(1,nDirections));
        CBL = kron(CU,ones(1,nDirections));
        CBR = kron(ones(1,nDirections),CU);
        
        % define within-bimanual contrast
        CBw = 3*ones(nDirections^2,nDirections^2);
        thetaR = kron(ones(nDirections,1),[1:nDirections]');
        thetaL = kron([1:nDirections]',ones(nDirections,1));
        
        for i=1:length(thetaL)
            for j=1:length(thetaR)
                if (thetaL(i)==extrinsic(thetaR(i))&(thetaL(j)==extrinsic(thetaR(j))))
                    CBw(i,j) = 1;
                end
                if (thetaL(i)==intrinsic(thetaR(i))&(thetaL(j)==intrinsic(thetaR(j))))
                    CBw(i,j) = 2;
                end
            end
        end
        CBw(logical(eye(length(thetaL)))) = 0;
        
        % check contrasts
        %         figure('position',[5 5 15,30]);
        %         subplot(1,3,1); imagesc(CU); xlabel('R'); ylabel('L'); axis equal
        %         subplot(1,3,2); imagesc(CBL); xlabel('Bi-L'); ylabel('R'); axis equal
        %         subplot(1,3,3); imagesc(CBR); xlabel('Bi-R'); ylabel('L'); axis equal
        
        % apply contrasts: loop over subjects,hemispheres,rois
        %subjects = unique(T.subj);
        for s=subjects
            for h=hemi
                for reg = regions
                    T_ = getrow(T,T.hemis==h&T.region==reg&T.subj==s);
                    
                    RDM = rsa_squareRDM(T_.RDM);
                    
                    % pick corresponding parts of RDM
                    U   = RDM(1:nDirections,nDirections+1:2*nDirections); % left uni vs right uni part
                    ULB = RDM(1:nDirections,2*nDirections+1:nConditions); % left uni vs bi part
                    URB = RDM(nDirections+1:2*nDirections,2*nDirections+1:nConditions); % right uni vs bi part
                    B   = RDM(2*nDirections+1:nConditions,2*nDirections+1:nConditions); % bimanual
                    
                    % ssqrt
                    if sqrtTransform
                        U = ssqrt(U);
                        ULB = ssqrt(ULB);
                        URB = ssqrt(URB);
                        ylim = [-0.01 0.08];
                        ylim1 = [-0.01 0.3];
                    else
                        ylim = [-0.001 0.07];
                        ylim1 = [-0.001 0.09];
                    end
                    
                    % overall distance
                    S.Uall = mean(mean(U.*Coverall));
                    S.ULBall = mean(mean(ULB.*Coverall_Bi));
                    S.URBall = mean(mean(URB.*Coverall_Bi));
                    
                    % uni contra vs ipsi
                    S.Uex = mean(mean(U(CU==1))); % parallel in extrinsic space (out-out, back-bacl)
                    S.Uin = mean(mean(U(CU==2))); % parallel in intrinsic space
                    S.Uex2 = mean(mean(U(CU==3))); % parallel in extrinsic space (out-back,back-out)
                    S.Uin2 = mean(mean(U(CU==6))); % parallel in intrinsic space
                    S.Uun = mean(mean(U(CU==0))); % unrelated movements
                    
                    % uni vs bi-ipsi
                    S.ULBex = mean(mean(ULB(CBR==1)));
                    S.ULBin = mean(mean(ULB(CBR==2)));
                    S.ULBex2 = mean(mean(ULB(CBR==3)));
                    S.ULBin2 = mean(mean(ULB(CBR==6)));
                    S.ULBun = mean(mean(ULB(CBR==0)));
                    
                    S.URBex = mean(mean(URB(CBL==1)));
                    S.URBin = mean(mean(URB(CBL==2)));
                    S.URBex2 = mean(mean(URB(CBL==3)));
                    S.URBin2 = mean(mean(URB(CBL==6)));
                    S.URBun = mean(mean(URB(CBL==0)));
                    
                    % uni vs bi-contra (too tribial?)
                    S.ULBcon = mean(mean(ULB(CBL==1)));
                    S.URBcon = mean(mean(URB(CBR==1)));
                    
                    % within-bi? (bimanual-parallel-extrinsic, bimanual-parallel-intrinsic)
                    S.BwEx = mean(mean(B(CBw==1))); % extrinsically parallel bimanual combinations
                    S.BwIn = mean(mean(B(CBw==2))); % intrinsically parallel bimanual combinations
                    S.BwUn = mean(mean(B(CBw==3))); % unrelated bimanual combinations
                    
                    S.subj = s;
                    S.hemis = h;
                    S.region = reg;
                    
                    R = addstruct(R,S);
                end
            end
        end
        CAT.fillcolor   = {[0.3 0.3 1],[1,0.3,0.3],'w',[0.3 0.3 1],[1 0.3 0.3]};%mat2cell(jet(length(postfix)),ones(length(postfix),1),3);
        CAT.linecolor   = repmat({'k'},size(CAT.fillcolor,1),1);
        CAT.linewidth   = repmat({1},size(CAT.linecolor,1),1);
        CAT.mediancolor = repmat({'k'},size(CAT.linecolor,1),1);
        CAT.medianwidth = repmat({1},size(CAT.linecolor,1),1);
        CAT.markersize  = repmat({5},size(CAT.medianwidth));
        CAT.markertype  = repmat({'o'},size(CAT.medianwidth));
        CAT.markerfill  = CAT.fillcolor;
        CAT.markercolor = CAT.fillcolor;
        CAT.facecolor   = CAT.fillcolor;
        CAT.linewidth   = {1.5};
        CAT.errorwidth  = {1.5};
        CAT.errorcolor  = {[0 0 0]};
        
        % scatter plot uex and uex2, etc
        if (0)
            for h=hemi
                figure('position',[5 5 25,25*sqrt(2)],'name','Corr');
                for reg = regions
                    R_ = getrow(R,R.hemis==h&R.region==reg);
                    
                    Xdata = [R_.Uex;R_.Uin];
                    Ydata = [R_.Uex2;R_.Uin2];
                    split = [ones(size(R_.Uex));2*ones(size(R_.Uin))];
                    
                    subplot(length(regions)/2,2,reg); warning off;
                    [r2,~,~,p]=scatterplot(Xdata,Ydata,'split',split,'CAT',CAT,...
                        'intercept',1,'regression','linear','printcorr'); axis square
                    xlabel('out-out/back-back pairs');
                    ylabel('out-back/back-out pairs');
                    
                    title({sprintf('%s-%s',hemName{h},regname{reg}),...
                        sprintf('p=%1.4f',p(end,:))});
                end
            end
        end
        R.extrinsic = R.Uun-R.Uex;
        R.intrinsic = R.Uun-R.Uin;
        variables   = {'extrinsic','intrinsic'};
        labels      = {'Ex','In'};
        fillcolor   = {red,blue};
        bmw1_imana('ROI_RDM_plot_testzero2',R,variables,...
            'fillcolor',fillcolor,'labels',labels,...
            'yrange',[-0.05 0.2],'yname','\Deltad (a.u.)',...
            'split','on');
        
        
        % see group statistics
        %         variables = {'U'}; %'ULB', 'URB'
        %         postfix = {'ex','in','un'};%,'ex2','in2'};
        %         fillcolor = {red,blue,'w'};
        %         labels = {'Ex','In','Un'};
        %         bmw1_imana('ROI_RDM_plot_compare',R,variables,postfix,3,...
        %             'fillcolor',fillcolor,'labels',labels);
        %
        % within-bimanual
        %         variables = {'Bw'};
        %         postfix = {'Ex','In','Un'};
        %         fillcolor = {darkred,darkblue,'w'};
        %         bmw1_imana('ROI_RDM_plot_compare',R,variables,postfix,3,...
        %             'fillcolor',fillcolor);
        
        varargout = {R};
    case 'ROI_RDM_plot_testzero'
        R = varargin{1};
        variables = varargin{2};
        plotfcn = varargin{3};
        figName = varargin{4};
        fillcolor = {'w'};
        labels = variables;
        yname = 'Average distance (a.u.)';
        yrange = 'auto';
        vararginoptions(varargin(5:end),{'fillcolor','labels','yname','yrange'});
        
        hemi = unique(R.hemis)';
        regions = unique(R.region)';
        
        
        CAT.fillcolor   = fillcolor;%{'w'};%mat2cell(gray(length(variables)),ones(length(variables),1),3);
        CAT.linecolor   = repmat({'k'},size(CAT.fillcolor,1),1);
        CAT.linewidth   = repmat({1},size(CAT.linecolor,1),1);
        CAT.mediancolor = repmat({'k'},size(CAT.linecolor,1),1);
        CAT.medianwidth = repmat({1},size(CAT.linecolor,1),1);
        CAT.markersize  = repmat({2},size(CAT.medianwidth));
        CAT.markertype  = repmat({'o'},size(CAT.medianwidth));
        CAT.markerfill  = CAT.fillcolor;
        CAT.markercolor = CAT.fillcolor;
        CAT.facecolor   = CAT.fillcolor;
        CAT.gapwidth    = repmat({0.1},length(variables),1);
        CAT.barwidth    = repmat({0.8},length(variables),1);
        CAT.linewidth   = {1.5};
        CAT.errorwidth  = {1.5};
        CAT.errorcolor  = {[0 0 0]};
        
        plotvals = [];
        for v=1:length(variables)
            plotvals = [plotvals,R.(variables{v})];
        end
        R.plotvals = plotvals;
        
        PFDR = [];
        for h=hemi
            c = 1;
            figure('position',[5 5 15,30],'name',figName);
            for reg = regions
                R_ = getrow(R,R.hemis==h&R.region==reg);
                
                subplot(length(regions)/2,2,reg);
                
                warning off
                switch (plotfcn)
                    case 'boxplot'
                        myboxplot([],R_.plotvals,'CAT',CAT);
                        drawline(0,'dir','horz');
                        set(gca,'xticklabel',labels);
                    case 'barplot'
                        [xpos,meanvals,ebars] = barplot([],R_.plotvals,...
                            'CAT',CAT,'barwidth',0.8);%hold off
                        
                        switch ischar(yrange)
                            case 1
                                ylim = [min(meanvals),max(meanvals)]+...
                                    [-ebars(meanvals==min(meanvals)),ebars(meanvals==max(meanvals))];
                                yrange_ = 0.1*(ylim(1)-ylim(1));
                                ylim = [ylim(1)-yrange_,ylim(2)+yrange_];
                            otherwise
                                ylim = yrange;
                        end
                        
                        for v=1:numel(variables);
                            [t(v),p(v)] = ttest_mc(R_.(variables{v}),0,1,'onesample');
                            drawline(0,'dir','horz');
                            %p(v) = signrank(R_.(variables{v}));
                            ypos(v) = 1.3*(meanvals(v)+ebars(v));
                            if p(v)<0.05
                                text(xpos(v),ypos(v),sprintf('p=%1.4f',p(v)),...
                                    'horizontalalignment','center',...
                                    'color','r');
                            end
                        end
                        set(gca,'xtick',xpos,'xticklabel',labels);
                        set(gca,'ylim',ylim);
                        PFDR = [PFDR,p];
                end
                warning on
                title(sprintf('%s-%s',hemName_s{h},regname{reg}));
                set(gca,'tickdir','out');
                if c==1;
                    ylabel(yname);
                else
                    ylabel('');
                end
                c = c+1;
            end
        end
        % p-FDR
        [pFDR] = FDR(PFDR,0.05);
        fprintf('p(FDR=0.05)=%1.4f\n',pFDR);
    case 'ROI_RDM_plot_testzero2' % horizontal bar plot
        R = varargin{1};
        variables = varargin{2};
        fillcolor = {'w'};
        labels = variables;
        yname = 'Average distance (a.u.)';
        yrange = [];
        split = 'off';
        vararginoptions(varargin(3:end),{'fillcolor','labels','yname','yrange','split'});
        
        if ~isempty(yrange); ylim=yrange; end
        hemi = unique(R.hemis)';
        regions = unique(R.region)';
        
        CAT.fillcolor   = fillcolor;%{'w'};%mat2cell(gray(length(variables)),ones(length(variables),1),3);
        CAT.linecolor   = repmat({'k'},size(CAT.fillcolor,1),1);
        CAT.linewidth   = repmat({1},size(CAT.linecolor,1),1);
        CAT.mediancolor = repmat({'k'},size(CAT.linecolor,1),1);
        CAT.medianwidth = repmat({1},size(CAT.linecolor,1),1);
        CAT.markersize  = repmat({2},size(CAT.medianwidth));
        CAT.markertype  = repmat({'o'},size(CAT.medianwidth));
        CAT.markerfill  = CAT.fillcolor;
        CAT.markercolor = CAT.fillcolor;
        CAT.facecolor   = CAT.fillcolor;
        CAT.gapwidth    = repmat({0.1},length(variables),1);
        CAT.barwidth    = repmat({0.8},length(variables),1);
        CAT.linewidth   = {1.5};
        CAT.errorwidth  = {1.5};
        CAT.errorcolor  = {[0 0 0]};
        
        % loop over variables
        PFDR = [];
        for v=1:numel(variables)
            plotval = R.(variables{v});
            % bar graph
            figure('position',[5 5 sqrt(2)*10 10],'name',labels{v});
            subplot(1,2,1)
            [xl,ml] = barplot([R.region],[plotval],...
                'subset',R.hemis==1,...
                'CAT',CAT,'facecolor',CAT.fillcolor{v});
            set(gca,'ylim',ylim,'xticklabel',{});
            set(gca,'View',[90 90],'Xaxislocation','top','YDir','reverse',...
                'tickdir','out','ticklength',[0.03,0.01]);
            ylabel(yname);title('Left hem')
            
            subplot(1,2,2)
            [xr,mr] = barplot([R.region],[R.(variables{v})],...
                'subset',R.hemis==2,...
                'CAT',CAT,'facecolor',CAT.fillcolor{v});
            set(gca,'ylim',ylim,'xticklabel',regname(regions));
            set(gca,'view',[90 90],'tickdir','out','ticklength',[0.03,0.01]);
            ylabel(yname);title('Right hem')
            
            xvars = {xl,xr};
            mvars={ml,mr};
            for h=hemi
                subplot(1,2,h);
                p=[];
                for reg=regions
                    X = plotval(R.hemis==h&R.region==reg);
                    [~,p(end+1)] = ttest_mc(X,0,1,'onesample');
                    %p(end+1) = signrank(X);
                    xpos = xvars{h}(reg);
                    ypos = mvars{h}(reg)*1.2;
                    if p(end)<0.05
                        text(xpos,ypos,sprintf('p=%1.4f',p(end)),...
                            'horizontalalignment','center',...
                            'color','r');
                    end
                end
                PFDR = [PFDR,p];
            end
        end
        % p-FDR
        [pFDR] = FDR(PFDR,0.05);
        fprintf('p(FDR=0.05)=%1.4f\n',pFDR);
    case 'ROI_RDM_plot_compare'
        R = varargin{1};
        variables = varargin{2};
        postfix = varargin{3};
        com = varargin{4};
        fillcolor = {'w'};
        labels = [];
        vararginoptions(varargin(5:end),{'fillcolor','labels'});
        if isempty(labels); labels=postfix; end
        hemi = unique(R.hemis)';
        regions = unique(R.region)';
        
        CAT.fillcolor   = fillcolor;%{[0.3 0.3 1],[1,0.3,0.3],'w',[0.3 0.3 1],[1 0.3 0.3]};%mat2cell(jet(length(postfix)),ones(length(postfix),1),3);
        CAT.linecolor   = repmat({'k'},size(CAT.fillcolor,1),1);
        CAT.linewidth   = repmat({1},size(CAT.linecolor,1),1);
        CAT.mediancolor = repmat({'k'},size(CAT.linecolor,1),1);
        CAT.medianwidth = repmat({1},size(CAT.linecolor,1),1);
        CAT.markersize  = repmat({5},size(CAT.medianwidth));
        CAT.markertype  = repmat({'o'},size(CAT.medianwidth));
        CAT.markerfill  = CAT.fillcolor;
        CAT.markercolor = CAT.fillcolor;
        CAT.facecolor   = CAT.fillcolor;
        CAT.linewidth   = {1.5};
        CAT.errorwidth  = {1.5};
        CAT.errorcolor  = {[0 0 0]};
        
        PFDR = [];
        for v=1:length(variables)
            for h=hemi
                figure('position',[5 5 15,30],'name',variables{v});
                c = 1;
                for reg = regions
                    fprintf('%s-%s-%s\n',hemName{h},regname{reg},variables{v});
                    R_ = getrow(R,R.hemis==h&R.region==reg);
                    
                    plotdataY = [];plotdataX = [];
                    for p=1:length(postfix)
                        plotdataY = [plotdataY;R_.(sprintf('%s%s',variables{v},postfix{p}))];
                        plotdataX = [plotdataX;repmat(p,length(R_.subj),1)];
                    end
                    
                    subplot(length(regions)/2,2,reg); warning off;
                    [xpos,meandata,errors] = barplot([],plotdataY,'CAT',CAT,...
                        'split',plotdataX,'gapwidth',[0.1 0.1 0.1],'barwidth',0.8);hold on; warning on;
                    hold on;
                    drawline(0,'dir','horz');
                    %lineplot(plotdataX,plotdataY,'CAT',CAT);
                    
                    set(gca,'xtick',xpos,'xticklabel',labels);
                    set(gca,'tickdir','out');
                    title(sprintf('%s-%s',hemName{h},regname{reg}));
                    pval=[];
                    for p=setdiff([1:length(postfix)],com)
                        A = R_.(sprintf('%s%s',variables{v},postfix{com}));
                        B = R_.(sprintf('%s%s',variables{v},postfix{p}));
                        
                        [t,pval(end+1)] = ttest_mc(A,B,1,'paired');
                        ypos = 1.5*meandata(com)+errors(com);
                        
                        if pval(end)<0.05
                            text(0.5*(xpos(p)+xpos(com)),ypos*(1-0.08*p),sprintf('%1.4f',pval(end)),...
                                'color','r','horizontalalignment','center');
                            drawline(ypos*(1-0.08*p),'dir','horz','lim',[xpos(p),xpos(com)]);
                        end
                        fprintf('%s-p=%2.5f\n',postfix{p},pval(end));
                    end
                    PFDR = [PFDR,pval];
                    set(gca,'ylim',[min([0,min(meandata)-errors(meandata==min(meandata))]),ypos]);
                    if c==1;
                        ylabel('Average distance (a.u.)');
                    else
                        ylabel('');
                    end
                    c = c+1;
                end
            end
        end
        % p-FDR
        [pFDR] = FDR(PFDR,0.05);
        fprintf('p(FDR=0.05)=%1.4f\n',pFDR);
    case 'ROI_RDM_plot_conip' % re-organize in terms of contra and ipsi
        R = varargin{1};
        
        regions = unique(R.region)';
        
        var = {'ConU','IpsU','ConB','IpsB'};
        T = [];P = [];
        for reg=regions
            lh = R.hemis==1&R.region==reg;
            rh = R.hemis==2&R.region==reg;
            
            % Contra
            t.ConU = mean([R.UL(rh),R.UR(lh)],2);
            t.ConB = mean([R.Bi2L(rh),R.Bi2R(lh)],2);
            
            % Ipsi
            t.IpsU = mean([R.UR(rh),R.UL(lh)],2);
            t.IpsB = mean([R.Bi2R(rh),R.Bi2L(lh)],2);
            
            % t-test
            for v=1:numel(var)
                [~,p.p(v)] = ttest_mc(t.(var{v}),0,1,'onesample');
            end
            
            t.region = R.region(lh);
            t.subj = R.subj(lh);
            p.region = reg;
            
            P = addstruct(P,p);
            T = addstruct(T,t);
        end
        % plot
        figure('position',[5 5 20 20])
        for reg=regions
            subplot(2,4,reg)
            % barplot
            [x,m,e] = barplot(T.region,[T.ConU,T.IpsU,T.ConB,T.IpsB],'gapwidth',[0.4,0.1],...
                'barwidth',0.8,'facecolor',{red blue darkred darkblue},'subset',T.region==reg);
            title(regname{reg});
            for v=1:numel(var)
                if P.p(P.region==reg,v)<0.05
                    text(x(v),[m(v)+e(v)+mean(e)],sprintf('%1.4f',P.p(P.region==reg,v)),...
                        'color','r','horizontalalignment','center');
                end
            end
            
        end
        fprintf('p(FDR==0.05)=%1.4f',FDR(P.p(:),0.05));
        pivottablerow(P.region,P.p,'double(x<0.05)+double(x<0.01)+double(x<0.001)');
        varargout = {};
        
    case 'ROI_IPM_compare_uniLR'                    % Compare unimanual left and right (c.f. Diedrichsen et al.,  Cereb Cortex 2013)
        hemi = [1:2];
        regions = [1:8];
        glm = 10;
        subjects = [1:7];
        
        vararginoptions(varargin,{'hemi','regions','subjects','glm'});
        
        % load data
        T = load(fullfile(regDir,sprintf('IPMRDMsep.glm%d.mat',glm)));
        T = getrow(T,ismember(T.subj,subjects));
        R = [];
        D = [];
        
        % loop over subjects,hemispheres,rois
        %subjects = unique(T.subj);
        for s=subjects
            for h=hemi
                for reg = regions
                    T_ = getrow(T,T.hemis==h&T.region==reg&T.subj==s);
                    
                    IPM = rsa_squareIPM(T_.IPM);
                    IPM = IPM(1:nConditions,1:nConditions);
                    
                    % pick IPM for unimanual-LR pairs
                    UL  = IPM(1:nDirections,1:nDirections);
                    UR  = IPM(nDirections+1:2*nDirections,nDirections+1:2*nDirections);
                    ULR = IPM(1:nDirections,nDirections+1:2*nDirections);
                    
                    % get phi
                    [S.ul,S.dtheta]   = RDM2Tuning(UL,'uni',1);
                    S.ur            = RDM2Tuning(UR,'uni',1);
                    S.ulr           = RDM2Tuning(ULR,'uni',1);
                    
                    S.subj = repmat(s,length(S.dtheta),1);
                    S.hemis = repmat(h,length(S.dtheta),1);
                    S.region = repmat(reg,length(S.dtheta),1);
                    R = addstruct(R,S);
                    
                    % index for evaluatoin: distance from Phi(0)
                    remove180 = ~ismember(S.dtheta,[-180,180]);
                    zero      = S.dtheta==0;
                    remove1800= ~ismember(S.dtheta,[-180,0,180]);
                    N       = sum(~ismember(S.dtheta,[-180,0,180]));
                    
                    ul      = S.ul-S.ul(zero);
                    ur      = S.ur-S.ur(zero);
                    ulr     = S.ulr-S.ulr(zero);
                    xlu     = -ul(remove1800);
                    xru     = -ur(remove1800);
                    xlru    = -ulr(remove1800);
                    ul      = sum(xlu)/N;
                    ur      = sum(xru)/N;
                    ulr     = sum(xlru)/N;
                    
                    % correlation between Phi(x)-Phi(0)
                    D_.corrLR = (corrN(ssqrt(xlu).*ssqrt(xru),xlru));
                    D_.corrLR2= ulr/(ssqrt(ul)*ssqrt(ur));
                    
                    D_.subj = s;
                    D_.hemis = h;
                    D_.region= reg;
                    D = addstruct(D,D_);
                end
            end
        end
        
        % plot Phi
        for i=1:2
            for h=hemi
                figure('position',[10 10 25 25],'name',['Phi(uniLR)']);
                y = [R.ul,R.ur,R.ulr];
                leg = {'PhiL(uL)','PhiL(uR)','PhiL(uLR)'};
                
                for reg = regions
                    subplot(length(regions)/2,2,reg);
                    lineplot(R.dtheta,y,'subset',R.hemis==h&R.region==reg,...
                        'style_shade','leg',leg);
                    xlabel('\Delta\theta (across condition)');
                    title(sprintf('%s-%s',regname{reg},hemName{h}))
                end
            end
        end
        
        % scatter plot between cross-phi and sqrt(phi(uni).*phi(bi)) (1)
        %         R1 = [];
        %         idx = ismember(R.dtheta,[0]);
        %         for i=1:2
        %         for h=hemi
        %             figure('position',[10 10 15 25],'name',['Phi(unixbi)']);
        %
        %             switch i
        %                 case 1 % left
        %                     x = R.l;
        %                     y = mean(R.ubl,2);
        %                 case 2 % right
        %                     x = R.r;
        %                     y = mean(R.ubr,2);
        %             end
        %
        %             for reg = regions
        %                 subplot(length(regions)/2,2,reg);
        %                 [r1.R2,r1.B,r1.t,r1.p] = ...
        %                     scatterplot(x,y,'subset',R.hemis==h&R.region==reg&idx,...
        %                     'identity','regression','linear','intercept',0,...
        %                     'printslope','printpval','draworig');
        %                 xlabel('sqrt(Phi_u*Phi_b)');
        %                 ylabel('Phi_{ub}')
        %                 title(sprintf('%s-%s',regname{reg},hemName{h}))
        %
        %                 r1.hand = i;
        %                 r1.hemis = h;
        %                 r1.region = reg;
        %
        %                 R1=addstruct(R1,r1);
        %             end
        %         end
        %         end
        
        % scatter plot between cross-phi and sqrt(phi(uni).*phi(bi)) (2)
        %         R3 = [];
        %         for i=1:2
        %         for h=hemi
        %             figure('position',[10 10 15 25],'name',['Phi(unixbi)']);
        %
        %             switch i
        %                 case 1 % left
        %                     x = ssqrt(D.ul).*ssqrt(D.bml);
        %                     y = D.ubl;
        %                 case 2 % right
        %                     x = ssqrt(D.ur).*ssqrt(D.bmr);
        %                     y = D.ubr;
        %             end
        %
        %             for reg = regions
        %                 subplot(length(regions)/2,2,reg);
        %                 [r2.R2,r2.B,r2.t,r2.p] =...
        %                     scatterplot(x,y,'subset',D.hemis==h&D.region==reg,...
        %                     'identity','regression','linear','intercept',0,...
        %                     'printcorr','draworig');
        %                 xlabel('sqrt(d_u*d_b)');
        %                 ylabel('d_{ub}')
        %                 title(sprintf('%s-%s',regname{reg},hemName{h}))
        %
        %                 r2.hand = i;
        %                 r2.hemis = h;
        %                 r2region = reg;
        %
        %                 R3=addstruct(R3,r2);
        %             end
        %         end
        %         end
        
        % plot correlation (unixbi)
        bmw1_imana('ROI_RDM_plot_testzero2',D,{'corrLR'},...
            'yname','Shared variance','yrange',[-0.5 1],'fillcolor',{'w','w'})
        
        varargout = {R};
    case 'ROI_IPM_compare_unibi'                    % Compare unimanual and bimanual (cf. Nozaki et al., 2006)
        % Compare overall similarity (e.g., correlation) between unimanual
        % and bimanual movements.
        % We can assess both 1) overall dissimilarity or 2) dissimilarities
        % for different opposite hand movement directions (here opposite
        % hand means the hand which was not used in unimanual condition)
        % 1) is in relation to Nozaki et al., NN 2006, and 2) is in relation
        % to Yokoi et al., sfn 2009
        %
        % Two ways of assessing correlation:
        % 1. use only diagonal of IPM (= Phi at dtheta=0)
        % 2. use modulation of Phi as a substitute to cancell the effect of common
        %    patterns
        
        hemi = [1:2];
        regions = [1:8];
        glm = 10;
        sqrtTransform = 1;
        ex180=1;
        subjects = [1:7];
        
        vararginoptions(varargin,{'hemi','regions','sqrtTransform','ex180','subjects','glm'});
        
        % load data (see 'ROI_sepGraw_glm9-12')
        T = load(fullfile(regDir,sprintf('IPMRDMsep.glm%d.mat',glm))); % mean pattern was subtracted for each condition separately
        T = getrow(T,ismember(T.subj,subjects));
        R = [];
        D = [];
        C = [];
        
        % make index for bimanual contra and ipsi distance
        C1 = [1:nDirections];
        C2 = ones(nDirections);
        CLB = kron(C2,C1);
        CRB = kron(C1,C2);
        
        % another way of calculating bimanual contra&ipsi
        Zcontra = kron(eye(nDirections),ones(nDirections,1));
        Zipsi = kron(ones(nDirections,1),eye(nDirections));
        Z = [zeros(2*nDirections),[Zcontra';Zipsi']];
        Z = [Z;[eye(2*nDirections),zeros(12,36)]];
        Z = bsxfun(@rdivide,Z,sum(Z,2));
        
        % loop over subjects,hemispheres,rois
        %subjects = unique(T.subj);
        for s=subjects
            for h=hemi
                for reg = regions
                    T_ = getrow(T,T.hemis==h&T.region==reg&T.subj==s);
                    
                    IPM = rsa_squareIPM(T_.IPM);
                    IPM = IPM(1:nConditions,1:nConditions); % choose ipm of interest (see 'ROI_sepGraw_glm9-12')
                    Bi2 = Z*IPM*Z';
                    
                    % pick IPM for unimanual-bimanual pairs
                    ULxB = IPM(1:nDirections,2*nDirections+1:end);
                    URxB = IPM(nDirections+1:2*nDirections,2*nDirections+1:end);
                    
                    % reorganize IPM in terms of uni directions under the
                    % fixed direction of opposite hand
                    for d=1:nDirections
                        UBL(:,:,d) = reshape(ULxB(CLB==d),nDirections,nDirections);
                        UBR(:,:,d) = reshape(URxB(CRB==d),nDirections,nDirections);
                    end
                    
                    % get phi (looking at the effect of opposite hand
                    % direction on uni and bi)
                    for i=1:nDirections
                        [S.ubl(:,i)] = RDM2Tuning(UBL(:,:,i),'uni',1);
                        [S.ubr(:,i),S.dtheta] = RDM2Tuning(UBR(:,:,i),'uni',1);
                        S.power_ubl(:,i) = repmat(trace(UBL(:,:,i)),size(S.dtheta));
                        S.power_ubr(:,i) = repmat(trace(UBR(:,:,i)),size(S.dtheta));
                    end
                    
                    % adjust dtheta with respect to opposite wrist
                    % direction?
                    
                    % get unimanual phi
                    UL = Bi2(2*nDirections+1:3*nDirections,2*nDirections+1:3*nDirections);
                    UR = Bi2(3*nDirections+1:4*nDirections,3*nDirections+1:4*nDirections);
                    [S.ul] = RDM2Tuning(UL,'uni',1);
                    [S.ur] = RDM2Tuning(UR,'uni',1);
                    
                    % get marginal bimanual phi
                    BmL = Bi2(1:nDirections,1:nDirections);
                    BmR = Bi2(nDirections+1:2*nDirections,nDirections+1:2*nDirections);
                    [S.bml] = RDM2Tuning(BmL,'uni',1);
                    [S.bmr] = RDM2Tuning(BmR,'uni',1);
                    
                    % calc sqrt(var(uni)*var(bi) -> maybe not correct
                    S.l = ssqrt(S.ul).*ssqrt(S.bml);
                    S.r = ssqrt(S.ur).*ssqrt(S.bmr);
                    
                    S.subj = repmat(s,length(S.dtheta),1);
                    S.hemis = repmat(h,length(S.dtheta),1);
                    S.region = repmat(reg,length(S.dtheta),1);
                    
                    % compare pattern strength between uni and bi
                    S.power_u = repmat(trace(Bi2(2*nDirections+1:end,2*nDirections+1:end)),...
                        size(S.subj));
                    S.power_b = repmat(trace(Bi2(1:2*nDirections,1:2*nDirections)),...
                        size(S.subj));
                    
                    R = addstruct(R,S);
                    
                    
                    % index for evaluatoin: distance from Phi(0)
                    remove180 = ~ismember(S.dtheta,[-180,180]);
                    zero      = S.dtheta==0;
                    remove1800= ~ismember(S.dtheta,[-180,0,180]);
                    N       = sum(~ismember(S.dtheta,[-180,0,180]));
                    
                    for i=1:nDirections
                        dubl    = S.ubl(:,i);
                        dubr    = S.ubr(:,i);
                        dubl    = dubl-dubl(zero);
                        dubr    = dubr-dubr(zero);
                        yl      = -dubl(remove1800);
                        yr      = -dubr(remove1800);
                        dubl    = sum(yl)/N;
                        dubr    = sum(yr)/N;
                        ul      = S.ul-S.ul(zero);
                        ur      = S.ur-S.ur(zero);
                        xlu     = -ul(remove1800);
                        xru     = -ur(remove1800);
                        ul      = sum(xlu)/N;
                        ur      = sum(xru)/N;
                        bml     = S.bml-S.bml(S.dtheta==0);
                        bmr     = S.bmr-S.bmr(S.dtheta==0);
                        xlb     = -bml(remove1800);
                        xrb     = -bmr(remove1800);
                        bml     = sum(xlb)/N;
                        bmr     = sum(xlb)/N;
                        
                        dubl_(i) = dubl;
                        dubr_(i) = dubr;
                        % correlation between Phi(x)-Phi(0)
                        c.corrL(i,1) = (corrN(ssqrt(xlb).*ssqrt(xlu),yl)); % corr(uniL,biL)
                        c.corrR(i,1) = (corrN(ssqrt(xrb).*ssqrt(xru),yr)); % corr(uniR,biR)
                        c.direction(i,1) = Directions(i);
                        c.yl(i,:) = yl';
                        c.yr(i,:) = yr';
                        c.xl(i,:) = [ssqrt(xlb).*ssqrt(xlu)]';
                        c.xr(i,:) = [ssqrt(xrb).*ssqrt(xru)]';
                    end
                    c.subj = repmat(s,6,1);
                    c.hemis = repmat(h,6,1);
                    c.region = repmat(reg,6,1);
                    C = addstruct(C,c);
                    
                    % distance calculated from Phi
                    D_.ul = 2*ul;
                    D_.ur = 2*ur;
                    D_.bml= 2*bml;
                    D_.bmr= 2*bmr;
                    D_.ubl= 2*dubl_;
                    D_.ubr= 2*dubr_;
                    
                    D_.subj = s;
                    D_.hemis = h;
                    D_.region= reg;
                    D = addstruct(D,D_);
                end
            end
        end
        
        % plot Phi
        %         for i=1:2
        %         for h=hemi
        %             figure('position',[10 10 25 25],'name',['Phi(unixbi)']);
        %
        %             switch i
        %                 case 1 % left
        %                     y = [R.ul,mean(R.ubl,2),R.bml];
        %                     leg = {'PhiL(uxu)','PhiL(uxb)','PhiL(bxb)'};
        %                 case 2 % right
        %                     y = [R.ur,mean(R.ubr,2),R.bmr];
        %                     leg = {'PhiR(uxu)','PhiR(uxb)','PhiR(bxb)'};
        %             end
        %
        %             for reg = regions
        %                 subplot(length(regions)/2,2,reg);
        %                 lineplot(R.dtheta,y,'subset',R.hemis==h&R.region==reg,...
        %                     'style_shade','leg',leg);
        %                 xlabel('\Delta\theta (across condition)');
        %                 title(sprintf('%s-%s',regname{reg},hemName{h}))
        %             end
        %         end
        %         end
        
        % scatter plot between cross-phi and sqrt(phi(uni).*phi(bi)) (1)
        R1 = [];
        for i=1:2
            for h=hemi
                figure('position',[10 10 15 25],'name',['Phi(unixbi)']);
                
                switch i
                    case 1 % left
                        x = mean(C.yl,2);
                        y = mean(C.xl,2);
                    case 2 % right
                        x = mean(C.yr,2);
                        y = mean(C.xr,2);
                end
                
                for reg = regions
                    subplot(length(regions)/2,2,reg);
                    [r1.R2,r1.B,r1.t,r1.p] = ...
                        scatterplot(x,y,'subset',C.hemis==h&C.region==reg,...
                        'identity','regression','linear','intercept',0,...
                        'printslope','printpval','draworig');
                    xlabel('sqrt(Phi_u*Phi_b)');
                    ylabel('Phi_{ub}')
                    title(sprintf('%s-%s',regname{reg},hemName{h}))
                    
                    r1.hand = i;
                    r1.hemis = h;
                    r1.region = reg;
                    
                    R1=addstruct(R1,r1);
                end
            end
        end
        
        % plot corrN on each opposite hand direction (collapse hemi and show ipsi and contra)
        C1 = getrow(C,C.hemis==1);
        C2 = getrow(C,C.hemis==2);
        
        C1.corr = [C1.corrL+C2.corrR]/2; % corr(uni,bi) for ipsilateral representation
        C1.ipsi = ones(size(C1.corr)); % identifier
        
        C2.corr = [C2.corrL+C1.corrR]/2; % corr(uni,bi) for contralateral representation
        C2.ipsi = zeros(size(C2.corr)); % identifier
        
        C3 = addstruct(C1,C2);
        
        label={'contra','ipsi'};
        figure('position',[5 5 50 50/sqrt(2)/2]);
        for reg=regions
            subplot(2,4,reg);
            if reg==1
                lineplot(C3.direction,C3.corr,'subset',C3.region==reg,...
                    'split',C3.ipsi,'style_shade','linecolor',{darkred darkblue},...
                    'shadecolor',{[0.9 0.9 0.9]},'leg',label,'transp',1);
                xlabel('Opposite wrist direction');
                ylabel('Correlation without intercept');
            else
                lineplot(C3.direction,C3.corr,'subset',C3.region==reg,...
                    'split',C3.ipsi,'style_shade','linecolor',{darkred darkblue},...
                    'shadecolor',{[0.9 0.9 0.9]},'transp',1);
                xlabel('');ylabel('');
            end
            drawline(0,'dir','horz');
            set(gca,'xlim',[-3,303],'ylim',[-0.5,1],'xtick',[0:60:300],...
                'tickdir','out','ticklength',[0.03 0.01]);
            
            title(regname{reg});
        end
        
        % plot summary of corrN averaged over direction of opposite hand
        C4 = tapply(C3,{'subj','hemis','region','ipsi'},{'corr','nanmean','name','corrUB'});
        figure('position',[5 5 20 20/sqrt(2)/2]);
        subplot(1,2,1);
        [x,y,e]=barplot([C4.hemis,C4.region],C4.corrUB,'subset',C4.ipsi==0,...
            'edgecolor',[0 0 0],...lineColor,...
            'facecolor',[0.5 0.5 0.5],...lineColor,...
            'barwidth',0.8);
        title('Contralateral representation');
        ylabel('corr(U,B)');
        set(gca,'ylim',[-0.2,1.2],'xticklabel',{});
        regs = regname(regions);
        for r=1:length(x)
            text(x(r),-0.2,regs{r},'horizontalalignment','right',...
                'rotation',90)
        end
        
        % calc t-statistics and p for one sample t-test
        t       = y./e;
        df      = length(subjects)-1;
        pval    = 1-tcdf(t,df);
        
        drawasterisk(pval,x,y+e+0.1);
        
        subplot(1,2,2);
        [x,y,e]=barplot([C4.hemis,C4.region],C4.corrUB,'subset',C4.ipsi==1,...
            'edgecolor',[0 0 0],...lineColor,...
            'facecolor',[0.5 0.5 0.5],...lineColor,...
            'barwidth',0.8);
        title('Ipsilateral representation');
        ylabel('corr(U,B)');
        set(gca,'ylim',[-0.2,1.2],'xticklabel',{});
        for r=1:length(x)
            text(x(r),-0.2,regs{r},'horizontalalignment','right',...
                'rotation',90)
        end
        
        % calc t-statistics and p for one sample t-test
        t       = y./e;
        df      = length(subjects)-1;
        pval    = 1-tcdf(t,df);
        
        drawasterisk(pval,x,y+e+0.1);
        
        
        varargout = {R,R1};
        
    case 'ROI_tuning_get_G'                         % Sort G by delta theta to get tuning function and display them
        direction = [0:60:300]';
        
        glm = 10;
        center = 1;
        %vararginoptions(varargin,{'glm','center'});
        
        % load data
        %fnameG = sprintf('G_adjhrf.glm%d.mat',glm);
        fnameG = sprintf('G_raw.glm%d.mat',glm);
        T = load(fullfile(regDir,fnameG));
        
        regions     = [1:8]';%unique(T.region);
        hemis       = unique(T.hemis);
        subjects    = unique(T.subj);
        
        vararginoptions(varargin,{'regions','hemis','subjects','center'})
        
        
        Tuning      = []; Tu = [];
        for s = subjects'
            for h=hemis'
                for reg=regions'
                    D=getrow(T,T.region==reg&T.hemis==h&T.subj==s);
                    
                    G = reshape(D.G,nConditions,nConditions);
                    
                    % symmetrize G
                    G = 0.5*(G+G');
                    
                    % split RDM into unimanual and bimnaual parts (UL,UR,Bi)
                    UL = G(1:nDirections,1:nDirections);
                    UR = G(nDirections+1:2*nDirections,nDirections+1:2*nDirections);
                    Bi = G(2*nDirections+1:end,2*nDirections+1:end);
                    
                    % center
                    if center
                        [UL,UL_] = doubleCenter(UL);
                        [UR,UR_] = doubleCenter(UR);
                        [Bi,Bi_] = doubleCenter(Bi);
                    else
                        UL_ = zeros(size(UL));
                        UR_ = zeros(size(UR));
                        Bi_ = zeros(size(Bi));
                    end
                    
                    % sort IPM by dtheta
                    [ul,diffs]  = RDM2Tuning(UL,'uni',1);
                    ur          = RDM2Tuning(UR,'uni',1);
                    bi          = RDM2Tuning(Bi,'bi',1);
                    ul_         = RDM2Tuning(UL_,'uni',1);
                    ur_         = RDM2Tuning(UR_,'uni',1);
                    bi_         = RDM2Tuning(Bi_,'bi',1);
                    
                    % get mean dist
                    Cuni = indicatorMatrix('allpairs',[1:6]);
                    Cbim = indicatorMatrix('allpairs',[1:36]);
                    mdL = nanmean(diag(Cuni*UL*Cuni'));
                    mdR = nanmean(diag(Cuni*UR*Cuni'));
                    mdB = nanmean(diag(Cbim*Bi*Cbim'));
                    tuning.meandist = repmat([mdL,mdR,mdB],length(diffs),1);
                    
                    % tuning data for plot
                    tuning.ul       = ul;
                    tuning.ur       = ur;
                    tuning.bi       = bi;
                    tuning.ul_      = ul_;
                    tuning.ur_      = ur_;
                    tuning.bi_      = bi_;
                    
                    tuning.dtheta   = diffs; ndiff = length(diffs);
                    tuning.subj     = repmat(s,ndiff,1);
                    tuning.region   = repmat(reg,ndiff,1);
                    tuning.hemi     = repmat(h,ndiff,1);
                    tuning.glm      = repmat(glm,ndiff,1);
                    tuning.center   = repmat(center,ndiff,1);
                    Tuning          = addstruct(Tuning,tuning);
                end
            end
        end
        % save result
        switch center
            case 0
                fname = fullfile(analyzeDir,'bmw1_Gtuning_forplot_uncentered.mat');
            case 1
                fname = fullfile(analyzeDir,'bmw1_Gtuning_forplot_centered.mat');
        end
        
        save(fname,'-struct','Tuning');
        
        % plot
        %bmw1_imana('ROI_tuning_plot',Tuning,'Ylim_uni',[-2*10e-4 5*10e-4],...
        %    'Ylim_bi',[-2*10e-4 4*10e-4],'arrange','axes','mode','IPM');
        bmw1_imana('ROI_tuning_plot',Tuning,'arrange','axes','mode','IPM',...
            'adjustcenter',1,'Ylim_uni','auto','Ylim_bi','auto',...
            'linecolor',{red blue darkred darkblue});
        
        varargout = {Tuning,Tu,T};
    case 'ROI_tuning_plot'
        Tuning = varargin{1};
        Ylim_uni = [-0.001,0.01];
        Ylim_bi = [0 0.2];
        arrange = 'subplot'; % 'subplot', 'axes'
        posAngle = [135 180 225 270 90 315 45 0];
        axsize = [0.2 0.14];
        sqrtTransform = 0;
        adjustcenter = 1;
        mode = 'IPM';
        linecolor = {purple gold darkpurple darkgold};
        
        if isempty(Tuning)
            Tuning = load(fullfile(analyzeDir,'bmw1_Gtuning_forplot_centered.mat'));%'bmw1_RDMtuning_forplot.mat'));
        end
        % get bimanual marginal
        Tuning = bmw1_imana('ROI_tuning_get_marginal',Tuning);
        
        
        regions     = unique(Tuning.region);
        hemis       = unique(Tuning.hemi);
        subjects    = unique(Tuning.subj);
        
        vararginoptions(varargin(2:end),{'regions','hemis','subjects','linecolor',...
            'Ylim_uni','Ylim_bi','arrange','sqrtTransform','mode','adjustcenter'});
        
        if sqrtTransform==1
            Ylim_uni = ssqrt(Ylim_uni);
            Ylim_bi = ssqrt(Ylim_bi);
        end
        
        switch mode
            case 'RDM'
                Tuning.ul = Tuning.rdmUL;
                Tuning.ur = Tuning.rdmUR;
                Tuning.bi = Tuning.rdmBi;
        end
        
        % plot summary
        Tuning = getrow(Tuning,ismember(Tuning.subj,subjects));
        
        if sqrtTransform==1
            Tuning.ul = ssqrt(Tuning.ul);
            Tuning.ur = ssqrt(Tuning.ur);
            Tuning.bi = ssqrt(Tuning.bi);
        end
        
        % adjust center
        if adjustcenter==1
            Tuning.ul = Tuning.ul-kron(Tuning.ul(Tuning.dtheta==0,:),ones(length(unique(Tuning.dtheta)),1));
            Tuning.ur = Tuning.ur-kron(Tuning.ur(Tuning.dtheta==0,:),ones(length(unique(Tuning.dtheta)),1));
            Tuning.Bi_ml = Tuning.Bi_ml-kron(Tuning.Bi_ml(Tuning.dtheta==0,:),ones(length(unique(Tuning.dtheta)),1));
            Tuning.Bi_mr = Tuning.Bi_mr-kron(Tuning.Bi_mr(Tuning.dtheta==0,:),ones(length(unique(Tuning.dtheta)),1));
        end
        
        % Unimanual
        switch (arrange)
            case 'subplot'
                figure('name','Unimanual trials',...
                    'position',[10 10 25*sqrt(2),25]);
            otherwise
        end
        c=1;
        for h=hemis'
            switch (arrange)
                case 'axes'
                    figure('name',sprintf('Unimanual trials-%s',hemName{h}),...
                        'position',[10 10 22*sqrt(2),22]);
                otherwise
            end
            
            for reg=regions'
                subset = Tuning.hemi==h&Tuning.region==reg;
                
                switch (arrange)
                    case 'subplot'
                        subplot(numel(hemis),numel(regions),c);
                    case 'axes'
                        R = 0.5*(1-axsize(1))*0.9;
                        orig = [R/0.9,R/0.9];[0.7-axsize(1),0.7-axsize(2)];
                        pos = [orig(1)+R*cosd(posAngle(reg)),orig(2)+R*sind(posAngle(reg)),...
                            axsize(1),axsize(2)];
                        axes('position',pos);
                end
                
                if c==1;
                    lineplot(Tuning.dtheta,[Tuning.ul,Tuning.ur,Tuning.Bi_ml,Tuning.Bi_mr],'subset',subset,...
                        'style_shade','linecolor',linecolor,...
                        'linestyle',{'-','-','--','--'},'errorbars',{''},...
                        'leg',{'UL','UR','BL','BR'},'leglocation','best','transp',1);
                else
                    lineplot(Tuning.dtheta,[Tuning.ul,Tuning.ur,Tuning.Bi_ml,Tuning.Bi_mr],'subset',subset,...
                        'style_shade','linecolor',linecolor,...
                        'linestyle',{'-','-','--','--'},'errorbars',{''},'transp',1);
                end
                drawline(0,'dir','horz')
                if c==1;
                    xlabel('\Delta\theta (degree)');
                    ylabel('Value (a.u.)');
                else
                    xlabel('');
                    ylabel('');
                end
                
                if ischar(Ylim_uni)
                    set(gca,'Xlim',[min(Tuning.dtheta), max(Tuning.dtheta)]);
                elseif reg==6|reg==7|reg==8;
                    set(gca,'Xlim',[min(Tuning.dtheta), max(Tuning.dtheta)]);
                else
                    set(gca,'Xlim',[min(Tuning.dtheta), max(Tuning.dtheta)],'Ylim',Ylim_uni);
                end
                title(sprintf('%s-%s',regname{reg},hemName_s{h}));
                set(gca,'tickdir','out','ticklength',[0.03 0.01]);
                
                c=c+1;
            end
        end
        
        % Bimanual
        for h=hemis'
            figure('name','Bimanual trials',...
                'position',[10 10 25*sqrt(2),25]);
            %set(gcf,'DefaultPatchLineStyle','non');
            c=1;
            for reg=regions'
                
                subset = Tuning.hemi==h&Tuning.region==reg;
                
                Dist = pivottablerow(Tuning.dtheta,Tuning.bi,'nanmean(x,1)','subset',subset);
                pivottablerow(Tuning.dtheta,Tuning.bi,'length','subset',subset);
                dtheta = unique(Tuning.dtheta);
                
                switch (arrange)
                    case 'subplot'
                        subplot(ceil(length(regions)/2),2,c);
                    case 'axes'
                        R = 0.5*(1-axsize(1))*0.9;
                        orig = [R/0.9,R/0.9];[0.7-axsize(1),0.7-axsize(2)];
                        pos = [orig(1)+R*cosd(posAngle(reg)),orig(2)+R*sind(posAngle(reg)),...
                            axsize(1),axsize(2)];
                        axes('position',pos);
                end
                
                colormap(jet);
                if ischar(Ylim_bi)
                    Ylim = [min(Dist(:)),max(Dist(:))];
                    mycontourf(dtheta,dtheta,Dist,'scale',Ylim); axis square
                    caxis(Ylim);
                elseif reg==6|reg==7|reg==8
                    Ylim = [min(Dist(:)),max(Dist(:))];
                    mycontourf(dtheta,dtheta,Dist,'scale',Ylim); axis square
                    caxis(Ylim);
                else
                    mycontourf(dtheta,dtheta,Dist,'scale',Ylim_bi); axis square
                    caxis(Ylim_bi);
                end
                %colorbar('location','eastoutside')
                xlabel('\Delta\theta_R');
                ylabel('\Delta\theta_L');
                set(gca,'Xlim',[min(Tuning.dtheta), max(Tuning.dtheta)],...
                    'Ylim',[min(Tuning.dtheta), max(Tuning.dtheta)])
                %set(gca,'XTick',[unique(Tuning.dtheta)],...
                %   'YTick',[unique(Tuning.dtheta)])
                set(gca,'XTick',[-180 -60 0 60 180],...
                    'YTick',[-180 -60 0 60 180])
                
                title(sprintf('%s-%s',regname{reg},hemName_s{h}));
                
                c=c+1;
            end
        end
    case 'ROI_tuning_plot_all1'
        % load data
        %T = load(fullfile(analyzeDir,'bmw1_RDMtuning_forplot.mat'));
        T = load(fullfile(analyzeDir,'bmw1_Gtuning_forplot_centered.mat'));
        
        subjects = unique(T.subj)';
        hemis = unique(T.hemi)';
        regions = unique(T.region)';
        
        vararginoptions(varargin,{'is180','subjects','hemis','regions'});
        
        Tuning = T;
        
        % get marginal tuning (just only taking center)
        Tuning = bmw1_imana('ROI_tuning_get_marginal',Tuning);
        
        % plot summary
        for h=hemis
            for reg=regions
                figname = sprintf('%s-%s',regname{reg},hemName{h});
                figure('name',figname,...
                    'position',[10 10 20 20]);
                subset = Tuning.hemi==h&Tuning.region==reg;
                
                % Bimanual
                %--------------------------------------------------------------------------
                Bi = pivottablerow(Tuning.dtheta,Tuning.bi,'nanmean(x,1)','subset',subset);
                dtheta = unique(Tuning.dtheta);
                position = [0.55 0.55 0.4 0.4];
                axes('position',position);
                mycontourf(dtheta,dtheta,Bi,'scale',[0 0.2]); %axis square
                %colorbar('location','eastoutside')
                xlabel('\Delta\theta_R (degree)');
                ylabel('\Delta\theta_L (degree)');
                
                % Bimanual marginal Left
                %--------------------------------------------------------------------------
                position = [0.275 0.55 0.15 0.4];
                axes('position',position);
                lineplot(Tuning.dtheta,Tuning.Bi_ml,'subset',subset,...
                    'style_shade','linecolor',[0.1 1 0.1],...
                    'shadecolor',[0.5 1 0.5]);
                set(gca,'ylim',[-0.05 0.2],'box','off');
                set(gca,'ytick',[0:0.1:0.2]);
                set(gca,'xlim',[-120 180]);
                set(gca,'xtick',[-120:60:180],'xticklabel',{'','','','',''})
                set(gca,'box','off','view',[-90 90]);
                set(gca,'xaxislocation','top','yaxislocation','right');
                xlabel('');ylabel('');
                
                % Unimanual Left
                %---------------------------------------------------------------------------
                position = [0.075 0.55 0.15 0.4];
                hold on;%axes('position',position);
                lineplot(Tuning.dtheta,Tuning.ul,'subset',subset,...
                    'style_shade','linecolor',[1 0.1 0.1],...
                    'shadecolor',[1 0.5 0.5]);
                set(gca,'ylim',[-0.05 0.2]);
                set(gca,'ytick',[0:0.1:0.2]);
                set(gca,'xlim',[-120 180]);
                set(gca,'xtick',[-120:60:180],'xticklabel',{'','','','',''})
                set(gca,'box','off','view',[-90 90]);
                set(gca,'xaxislocation','top','yaxislocation','right');
                xlabel('');ylabel('');
                
                % Bimanual marginal Right
                %--------------------------------------------------------------------------
                position = [0.55 0.275 0.4 0.15];
                axes('position',position);
                lineplot(Tuning.dtheta,Tuning.Bi_mr,'subset',subset,...
                    'style_shade','linecolor',[0.1 1 0.1],...
                    'shadecolor',[0.5 1 0.5]);
                set(gca,'ylim',[-0.05 0.2]);
                set(gca,'ytick',[0:0.1:0.2]);
                set(gca,'xlim',[-120 180]);
                set(gca,'xtick',[-120:60:180],'xticklabel',{'','','','',''})
                set(gca,'ydir','reverse','xaxislocation','top','box','off')
                set(gca,'xaxislocation','top','yaxislocation','right');
                xlabel('');ylabel('');
                
                % Unimanual Right
                %---------------------------------------------------------------------------
                position = [0.55 0.075 0.4 0.15];
                hold on;%axes('position',position);
                lineplot(Tuning.dtheta,Tuning.ur,'subset',subset,...
                    'style_shade','linecolor',[1 0.1 0.1],...
                    'shadecolor',[1 0.5 0.5]);
                set(gca,'ylim',[-0.05 0.2]);
                set(gca,'ytick',[0:0.1:0.2]);
                set(gca,'xlim',[-120 180]);
                set(gca,'xtick',[-120:60:180],'xticklabel',{'','','','',''})
                set(gca,'ydir','reverse','xaxislocation','top','box','off')
                set(gca,'xaxislocation','top','yaxislocation','right');
                xlabel('');ylabel('');
                
                % save fig
                if (0)
                    if ~exist(figDir,'dir');
                        mkdir(figDir);
                    end
                    fname = sprintf('Tuning_group_%s.fig',figname);
                    saveas(gcf,fullfile(figDir,fname),'fig');
                end
            end
        end
        
        vararegout = {Tuning};
    case 'ROI_tuning_plot_all2'
        % load data
        %T = load(fullfile(analyzeDir,'bmw1_RDMtuning_forplot.mat'));
        T = load(fullfile(analyzeDir,'bmw1_Gtuning_forplot_centered.mat'));
        
        subjects = unique(T.subj)';
        hemis = unique(T.hemi)';
        regions = unique(T.region)';
        
        vararginoptions(varargin,{'is180','subjects','hemis','regions'});
        
        Tuning = T;
        
        % Define color
        nDtheta = unique(Tuning.dtheta);
        lineColor = jet(numel(nDtheta));
        shadeColor= min(1,lineColor + 0.3);
        lineColor = mat2cell(lineColor,ones(numel(nDtheta),1),3);
        shadeColor= mat2cell(shadeColor,ones(numel(nDtheta),1),3);
        
        CAT.linecolor = lineColor;
        CAT.shadecolor= shadeColor;
        
        % plot summary
        for h=hemis
            for reg=regions
                figname = sprintf('%s-%s',regname{reg},hemName{h});
                figure('name',figname,...
                    'position',[10 10 20 20]);
                subset = Tuning.hemi==h&Tuning.region==reg;
                
                % Bimanual
                %--------------------------------------------------------------------------
                Bi = pivottablerow(Tuning.dtheta,Tuning.bi,'nanmean(x,1)','subset',subset);
                dtheta = unique(Tuning.dtheta);
                position = [0.55 0.55 0.4 0.4];
                axes('position',position);
                mycontourf(dtheta,dtheta,Bi,'scale',[0 0.2]); %axis square
                %colorbar('location','eastoutside')
                xlabel('\Delta\theta_R (degree)');
                ylabel('\Delta\theta_L (degree)');
                
                % Bimanual marginal Left
                %--------------------------------------------------------------------------
                position = [0.275 0.55 0.15 0.4];
                axes('position',position);
                lineplot(Tuning.dtheta,Tuning.bi,'subset',subset,...
                    'style_shade','CAT',CAT);
                %set(gca,'ylim',[-0.05 0.2],'box','off');
                set(gca,'ytick',[0:0.1:0.2]);
                set(gca,'xlim',[-120 180]);
                set(gca,'xtick',[-120:60:180],'xticklabel',{'','','','',''})
                set(gca,'box','off','view',[-90 90]);
                set(gca,'xaxislocation','top','yaxislocation','right');
                xlabel('');ylabel('');
                
                % Unimanual Left
                %---------------------------------------------------------------------------
                position = [0.075 0.55 0.15 0.4];
                hold on;%axes('position',position);
                lineplot(Tuning.dtheta,Tuning.ul,'subset',subset,...
                    'style_shade','linecolor',[1 0.1 0.1],...
                    'shadecolor',[1 0.5 0.5]);
                set(gca,'ylim',[-0.05 0.2]);
                set(gca,'ytick',[0:0.1:0.2]);
                set(gca,'xlim',[-120 180]);
                set(gca,'xtick',[-120:60:180],'xticklabel',{'','','','',''})
                set(gca,'box','off','view',[-90 90]);
                set(gca,'xaxislocation','top','yaxislocation','right');
                xlabel('');ylabel('');
                
                % Bimanual marginal Right
                %--------------------------------------------------------------------------
                position = [0.55 0.275 0.4 0.15];
                axes('position',position);
                lineplot(Tuning.dtheta,Tuning.Bi_mr,'subset',subset,...
                    'style_shade','linecolor',[0.1 1 0.1],...
                    'shadecolor',[0.5 1 0.5]);
                set(gca,'ylim',[-0.05 0.2]);
                set(gca,'ytick',[0:0.1:0.2]);
                set(gca,'xlim',[-120 180]);
                set(gca,'xtick',[-120:60:180],'xticklabel',{'','','','',''})
                set(gca,'ydir','reverse','xaxislocation','top','box','off')
                set(gca,'xaxislocation','top','yaxislocation','right');
                xlabel('');ylabel('');
                
                % Unimanual Right
                %---------------------------------------------------------------------------
                position = [0.55 0.075 0.4 0.15];
                hold on;%axes('position',position);
                lineplot(Tuning.dtheta,Tuning.ur,'subset',subset,...
                    'style_shade','linecolor',[1 0.1 0.1],...
                    'shadecolor',[1 0.5 0.5]);
                set(gca,'ylim',[-0.05 0.2]);
                set(gca,'ytick',[0:0.1:0.2]);
                set(gca,'xlim',[-120 180]);
                set(gca,'xtick',[-120:60:180],'xticklabel',{'','','','',''})
                set(gca,'ydir','reverse','xaxislocation','top','box','off')
                set(gca,'xaxislocation','top','yaxislocation','right');
                xlabel('');ylabel('');
                
                % save fig
                if (0)
                    if ~exist(figDir,'dir');
                        mkdir(figDir);
                    end
                    fname = sprintf('Tuning_group_%s.fig',figname);
                    saveas(gcf,fullfile(figDir,fname),'fig');
                end
            end
        end
        
        vararegout = {Tuning};
        
    case 'ROI_tuning_get_marginal'                  % Calculate marginalized distance tuning for bimanual movements
        Tuning = varargin{1};
        type = 'center';
        dtheta = unique(Tuning.dtheta);
        for s=unique(Tuning.subj)'
            for reg=unique(Tuning.region)'
                for h=unique(Tuning.hemi)'
                    idx = find(Tuning.subj==s&Tuning.region==reg&Tuning.hemi==h);
                    Bi = Tuning.bi(idx,:);
                    
                    switch (type)
                        case 'mean'
                            mleft = nanmean(Bi,2);
                            mright= nanmean(Bi,1);
                        case 'center'
                            mleft = Bi(:,dtheta==0);
                            mright = Bi(dtheta==0,:);
                    end
                    Tuning.Bi_ml(idx,:) = mleft;
                    Tuning.Bi_mr(idx,:) = mright';
                end
            end
        end
        
        varargout = {Tuning};
    case 'ROI_tuning_reliability_between'           % Calculate between-subject reliabitliy from tuning functino of distance
        includeZero=1;
        % load data
        T = load(fullfile(analyzeDir,'bmw1_RDMtuning_foranalysis.mat'));
        
        % loop over rois
        regions     = unique(T.region);
        hemis       = unique(T.hemi);
        subjects    = unique(T.subj);
        Tuning      = []; Tu = [];
        
        Ana=[];
        for h=hemis'
            for reg=regions'
                D = getrow(T,T.region==reg&T.hemi==h);
                
                switch includeZero
                    case 0
                        subset_u = ~all(D.ul==0);
                        subset_b = ~all(D.bi==0);
                    case 1
                        subset_u = ~all(isnan(D.ul));
                        subset_b = ~all(isnan(D.bi));
                end
                
                % uni left
                CL = corr(D.ul(:,subset_u)');
                ana.Rbtw_ul = nanmean(CL(tril(~logical(eye(size(D.ul,1))))));
                
                % uni right
                CR = corr(D.ur(:,subset_u)');
                ana.Rbtw_ur = nanmean(CR(tril(~logical(eye(size(D.ur,1))))));
                
                % bi
                CB = corr(D.bi(:,subset_b)');
                ana.Rbtw_bi = nanmean(CB(tril(~logical(eye(size(D.bi,1))))));
                
                ana.region = reg;
                ana.hemi = h;
                Ana = addstruct(Ana,ana);
            end
        end
        % show summary
        pivottablerow([Ana.hemi,Ana.region],[Ana.Rbtw_ul,Ana.Rbtw_ur,Ana.Rbtw_bi],'nanmean(x,1)');
        
        varargout = {Ana};
    case 'ROI_tuning_fit'                           % Estimate tuning parameters under assumption of Gaussian tuning
        is180 = 1;
        center = 1;
        fig = 0;
        vararginoptions(varargin,{'is180','center','fig'});
        % load data
        switch center
            case 0
                fname = fullfile(analyzeDir,'bmw1_Gtuning_forplot_uncentered.mat');
            case 1
                fname = fullfile(analyzeDir,'bmw1_Gtuning_forplot_centered.mat');
        end
        T1 = load(fname);
        
        subjects = unique(T1.subj)';
        hemis = unique(T1.hemi)';
        regions = unique(T1.region)';
        
        vararginoptions(varargin,{'is180','subjects','hemis','regions','center','fig'});
        
        Tu = [];
        for s=subjects
            for h=hemis
                for reg=regions
                    fprintf('%s-%s-%s\n',subj_name{s},hemName{h},regname{reg});
                    
                    % 0. get data
                    tmp     = getrow(T1,T1.subj==s&T1.hemi==h&T1.region==reg);
                    dtheta  = tmp.dtheta;
                    ul      = tmp.ul; T.UL = ul';
                    ur      = tmp.ur; T.UR = ur';
                    bi      = tmp.bi; T.Bi = vec(bi)';
                    range   = 1.3*[min([ul;ur;bi(:)]),max([ul;ur;bi(:)])];
                    
                    % 1. assess modulation (fit Gaussian)
                    % fit options
                    switch (is180)
                        case 1 % two gaussians
                            opt = fitoptions('Method','NonlinearLeastSquares',...
                                'Lower', [0 1.2*range(1) 20 0],...
                                'Upper', [1 1.2*range(2) 360 1],...
                                'StartPoint',[0.02 -0.01 20 0.02],...
                                'Algorithm', 'Levenberg-Marquardt',... % 'Gauss-Newton', 'Trust-Region'
                                'MaxIter', 10e9,...
                                'TolFun', 1e-9,...
                                'TolX', 1e-4...
                                );
                        case 0 % single gaussian
                            opt = fitoptions('Method','NonlinearLeastSquares',...
                                'Lower', [0 1.2*range(1) 20],...
                                'Upper', [1 1.2*range(2) 360],...
                                'StartPoint',[0.02 -0.01 20],...
                                'Algorithm', 'Levenberg-Marquardt',... % 'Gauss-Newton', 'Trust-Region'
                                'MaxIter', 10e9,...
                                'TolFun', 1e-9,...
                                'TolX', 1e-4...
                                );
                    end
                    
                    % fit function
                    switch (is180)
                        case 1 % two gaussians
                            fitstr = sprintf('(A*exp(-0.5*x*x/C/C)) + (D*exp(-0.5*(x-180)^2/C/C)) +(D*exp(-0.5*(x+180)^2/C/C)) + B');
                        case 0 % single gaussian
                            fitstr = sprintf('(A*exp(-0.5*x*x/C/C)+B)');
                    end
                    fitfun = fittype(fitstr,'options',opt);
                    
                    % set x and y
                    if is180==1 % use all data
                        x = dtheta;
                        yul = ul;
                        yur = ur;
                        ybi = bi;
                    else
                        x = dtheta(dtheta~=180&dtheta~=-180);
                        yul = ul(dtheta~=180&dtheta~=-180);
                        yur = ur(dtheta~=180&dtheta~=-180);
                        ybi = bi(dtheta~=180&dtheta~=-180,dtheta~=180&dtheta~=-180);
                    end
                    yblc = ybi(:,dtheta==0);
                    ybrc = ybi(dtheta==0,:)';
                    yblc = nanmean(ybi,2);
                    ybrc = nanmean(ybi,1)';
                    
                    
                    [para_ul, gof_ul] = fit(x,yul,fitfun);
                    [para_ur, gof_ur] = fit(x,yur,fitfun);
                    [para_brc, gof_brc] = fit(x,ybrc,fitfun); % tuning for marginal
                    [para_blc, gof_blc] = fit(x,yblc,fitfun); % tuning for marginal
                    
                    
                    T.R2 = [gof_ul.rsquare,gof_ur.rsquare,gof_blc.rsquare,gof_brc.rsquare];
                    T.A = [para_ul.A,para_ur.A,para_blc.A,para_brc.A]; % tuning parameters
                    T.B = [para_ul.B,para_ur.B,para_blc.B,para_brc.B];
                    T.C = [para_ul.C,para_ur.C,para_blc.C,para_brc.C];
                    switch (is180)
                        case 1 % two gaussians
                            T.D = [para_ul.D,para_ur.D,para_brc.D,para_brc.D];
                        case 0 % single gaussian
                            T.D = [NaN];
                    end
                    
                    if fig==1
                        figure('name',sprintf('%s-%s-%s\n',subj_name{s},hemName{h},regname{reg}));
                        subplot(2,2,1);plot(x,yul,'or');hold on;plot(para_ul);
                        title('uni-left');set(gca,'ylim',range);
                        subplot(2,2,2);plot(x,yur,'or');hold on;plot(para_ur);
                        title('uni-right');set(gca,'ylim',range);
                        subplot(2,2,3);plot(x,yblc,'or');hold on;plot(para_blc);
                        title('bi-left-center');set(gca,'ylim',range);
                        subplot(2,2,4);plot(x,ybrc,'or');hold on;plot(para_brc);
                        title('bi-right-center');set(gca,'ylim',range);
                        pause();close()
                    end
                    
                    % Calc mean distanace from phi [mean(d^2)=2*phi(zero)-sum(phi(nonzero))]
                    
                    T.subj = s;
                    T.hemis = h;
                    T.region = reg;
                    T.center = center;
                    T.is180 = is180;
                    T.meandist = nanmean(tmp.meandist,1);
                    
                    Tu = addstruct(Tu,T);
                end
            end
        end
        % save result
        switch center
            case 0
                fname = fullfile(analyzeDir,'bmw1_tuning_summary_uncentered.mat');
            case 1
                fname = fullfile(analyzeDir,'bmw1_tuning_summary_centered.mat');
        end
        save(fname,'-struct','Tu');
        
        varargout = {Tu};
    case 'ROI_tuning_fit_group'                     % Estimate tuning parameters under assumption of Gaussian tuning
        is180 = 1;
        center = 1;
        fig = 0;
        vararginoptions(varargin,{'is180','center','fig'});
        % load data
        switch center
            case 0
                fname = fullfile(analyzeDir,'bmw1_Gtuning_forplot_uncentered.mat');
            case 1
                fname = fullfile(analyzeDir,'bmw1_Gtuning_forplot_centered.mat');
        end
        T1 = load(fname);
        
        subjects = unique(T1.subj)';
        hemis = unique(T1.hemi)';
        regions = unique(T1.region)';
        
        vararginoptions(varargin,{'is180','subjects','hemis','regions','center','fig'});
        
        T1 = tapply(T1,{'region','hemi','dtheta'},...
            {'ul','nanmean','name','ul'},...
            {'ur','nanmean','name','ur'},...
            {'bi','nanmean','name','bi'});
        
        Tu = [];
        
        for h=hemis
            for reg=regions
                %fprintf('%s-%s-%s\n',subj_name{s},hemName{h},regname{reg});
                
                % 0. get data
                tmp     = getrow(T1,T1.hemi==h&T1.region==reg);
                dtheta  = tmp.dtheta;
                ul      = tmp.ul; T.UL = ul';
                ur      = tmp.ur; T.UR = ur';
                bi      = tmp.bi; T.Bi = vec(bi)';
                
                range = 1.3*[min([ul;ur;bi(:)]),max([ul;ur;bi(:)])];
                
                % 1. assess modulation (fit Gaussian)
                % fit options
                switch (is180)
                    case 1 % two gaussians
                        opt = fitoptions('Method','NonlinearLeastSquares',...
                            'Lower', [0 1.2*range(1) 20 0],...
                            'Upper', [1 1.2*range(2) 360 1],...
                            'StartPoint',[0.02 -0.01 20 0.02],...
                            'Algorithm', 'Levenberg-Marquardt',... % 'Gauss-Newton', 'Trust-Region'
                            'MaxIter', 10e9,...
                            'TolFun', 1e-9,...
                            'TolX', 1e-4...
                            );
                    case 0 % single gaussian
                        opt = fitoptions('Method','NonlinearLeastSquares',...
                            'Lower', [0 1.2*range(1) 20],...
                            'Upper', [1 1.2*range(2) 360],...
                            'StartPoint',[0.02 -0.01 20],...
                            'Algorithm', 'Levenberg-Marquardt',... % 'Gauss-Newton', 'Trust-Region'
                            'MaxIter', 10e9,...
                            'TolFun', 1e-9,...
                            'TolX', 1e-4...
                            );
                end
                
                % fit function
                switch (is180)
                    case 1 % two gaussians
                        fitstr = sprintf('(A*exp(-0.5*x*x/C/C)) + (D*exp(-0.5*(x-180)^2/C/C)) +(D*exp(-0.5*(x+180)^2/C/C)) + B');
                    case 0 % single gaussian
                        fitstr = sprintf('(A*exp(-0.5*x*x/C/C)+B)');
                end
                fitfun = fittype(fitstr,'options',opt);
                
                % set x and y
                if is180==1 % use all data
                    x = dtheta;
                    yul = ul;
                    yur = ur;
                    ybi = bi;
                else
                    x = dtheta(dtheta~=180&dtheta~=-180);
                    yul = ul(dtheta~=180&dtheta~=-180);
                    yur = ur(dtheta~=180&dtheta~=-180);
                    ybi = bi(dtheta~=180&dtheta~=-180,dtheta~=180&dtheta~=-180);
                end
                yblc = ybi(:,dtheta==0);
                ybrc = ybi(dtheta==0,:)';
                yblc = nanmean(ybi,2);
                ybrc = nanmean(ybi,1)';
                
                
                [para_ul, gof_ul] = fit(x,yul,fitfun);
                [para_ur, gof_ur] = fit(x,yur,fitfun);
                [para_brc, gof_brc] = fit(x,ybrc,fitfun); % tuning for marginal
                [para_blc, gof_blc] = fit(x,yblc,fitfun); % tuning for marginal
                
                
                T.R2 = [gof_ul.rsquare,gof_ur.rsquare,gof_blc.rsquare,gof_brc.rsquare];
                T.A = [para_ul.A,para_ur.A,para_blc.A,para_brc.A]; % tuning parameters
                T.B = [para_ul.B,para_ur.B,para_blc.B,para_brc.B];
                T.C = [para_ul.C,para_ur.C,para_blc.C,para_brc.C];
                switch (is180)
                    case 1 % two gaussians
                        T.D = [para_ul.D,para_ur.D,para_brc.D,para_brc.D];
                    case 0 % single gaussian
                        T.D = [NaN];
                end
                
                if fig==1
                    figure('name',sprintf('%s-%s\n',hemName{h},regname{reg}),...
                        'position',[15,15,20*sqrt(2),20]);
                    subplot(2,2,1);plot(x,yul,'or');hold on;plot(para_ul);
                    title('uni-left');set(gca,'ylim',range);
                    subplot(2,2,2);plot(x,yur,'or');hold on;plot(para_ur);
                    title('uni-right');set(gca,'ylim',range);
                    subplot(2,2,3);plot(x,yblc,'or');hold on;plot(para_blc);
                    title('bi-left-center');set(gca,'ylim',range);
                    subplot(2,2,4);plot(x,ybrc,'or');hold on;plot(para_brc);
                    title('bi-right-center');set(gca,'ylim',range);
                    %pause();close()
                end
                
                
                T.hemis = h;
                T.region = reg;
                T.center = center;
                T.is180 = is180;
                
                
                Tu = addstruct(Tu,T);
            end
        end
        
        % save result
        switch center
            case 0
                fname = fullfile(analyzeDir,'bmw1_tuning_summary_uncentered.mat');
            case 1
                fname = fullfile(analyzeDir,'bmw1_tuning_summary_centered.mat');
        end
        save(fname,'-struct','Tu');
        
        varargout = {Tu};
    case 'ROI_tuning_fit_plot'
        center = 1;
        vararginoptions(varargin,{'center'});
        
        % load data
        switch center
            case 0
                fname = fullfile(analyzeDir,'bmw1_tuning_summary_uncentered.mat');
            case 1
                fname = fullfile(analyzeDir,'bmw1_tuning_summary_centered.mat');
        end
        T = load(fname); %T=rmfield(T,'D');
        hemis = unique(T.hemis);
        regions = unique(T.region);
        
        % plot overall R-squared (V12 excluded)
        figure('name','Average R^2 (V12 excluded)','position',[10 10 20 20*sqrt(2)]);
        subplot(1,2,1)
        barplot(T.subj,(T.R2),'subset',T.region~=6&T.hemis==1,'leg',{'UL','UR','BL','BR'});
        xlabel('SN');ylabel('R^2');title('Left hem')
        
        subplot(1,2,2)
        barplot(T.subj,(T.R2),'subset',T.region~=6&T.hemis==2,'leg',{'UL','UR','BL','BR'});
        xlabel('SN');ylabel('R^2');title('Right hem')
        
        % plot each estimated parameters
        type = {'UL','UR','BL','BR'};
        fields = {'R2','A','B','C','D'};
        fieldsname = {'R^2','Amplitude (a.u.)','Offset (a.u.)','Tuning width (degree)','Amplitude for back movement (a.u.)'};
        for t=1:length(type)
            figure('name',type{t},'position',[10 10 20 20*sqrt(2)]);
            for f=1:length(fields)
                subplot(2,3,f);
                
                myboxplot([T.hemis,T.region],T.(fields{f})(:,t),'xtickoff');
                set(gca,'XTicklabel',regname(regions));
                title(fieldsname{f});
            end
        end
        
        % scatter plot for (amp1xamp2), (amp1xwidth)
        CAT.markercolor = {[1,0.3,0.3],[0.3,1,0.3]};
        CAT.markerfill = CAT.markercolor;
        CAT.markersize = 5;
        CAT.leglocation = 'best';
        for reg=regions'
            figure('name',sprintf('Amp1xAmp2xWidth-%s',regname{reg}),'position',[10 10 30*sqrt(2) 30]);
            
            for h=hemis'
                
                idx = T.hemis==h&T.region==reg;
                T_ = getrow(T,idx);
                typenum = repmat([1:4],size(T_.region,1),1);
                subj = repmat(T_.subj,1,4);
                x = reshape(T_.A,numel(T_.A),1);
                y = reshape(T_.D,numel(T_.D),1);
                z = reshape(T_.C,numel(T_.C),1);
                typenum = reshape(typenum,numel(typenum),1);
                subj = reshape(subj,numel(subj),1);
                
                %--- amp1 x amp2
                % ul and ur
                subplot(4,3,1+3*(h-1));
                scatterplot(x,y,'split',typenum,'subset',ismember(typenum,[1,2]),...
                    'label',subj,'identity','leg',{'UL','UR'},'CAT',CAT);
                xlabel('Amplitude (1)');
                ylabel({hemName{h},'Amplitude (2)'});title('UL vs UR')
                
                % ul and bi-l
                subplot(4,3,2+3*(h-1));
                scatterplot(x,y,'split',typenum,'subset',ismember(typenum,[1,3]),...
                    'label',subj,'identity','leg',{'UL','BiL'},'CAT',CAT);
                xlabel('Amplitude (1)');
                ylabel('Amplitude (2)');title('UL vs BiL')
                
                % ul and bi-l
                subplot(4,3,3+3*(h-1));
                scatterplot(x,y,'split',typenum,'subset',ismember(typenum,[2,4]),...
                    'label',subj,'identity','leg',{'UR','BiR'},'CAT',CAT);
                xlabel('Amplitude (1)');
                ylabel('Amplitude (2)');title('UR vs BiR')
                
                %--- amp1 x width
                % ul and ur
                subplot(4,3,7+3*(h-1));
                scatterplot(x,z,'split',typenum,'subset',ismember(typenum,[1,2]),...
                    'label',subj,'identity','leg',{'UL','UR'},'CAT',CAT);
                xlabel('Amplitude (1)');
                ylabel({hemName{h},'Width'});title('UL vs UR')
                
                
                % ul and bi-l
                subplot(4,3,8+3*(h-1));
                scatterplot(x,z,'split',typenum,'subset',ismember(typenum,[1,3]),...
                    'label',subj,'identity','leg',{'UL','BiL'},'CAT',CAT);
                xlabel('Amplitude (1)');
                ylabel('Width');title('UL vs BiL')
                
                % ul and bi-l
                subplot(4,3,9+3*(h-1));
                scatterplot(x,z,'split',typenum,'subset',ismember(typenum,[2,4]),...
                    'label',subj,'identity','leg',{'UR','BiR'},'CAT',CAT);
                xlabel('Amplitude (1)');
                ylabel('Width');title('UR vs BiR')
            end
        end
        
        varargout = {T};
    case 'ROI_tuning_summary'
        %         Check for encoding of movement directions for each ROI, separated for unimanual and bimanual trial types (if failed, exclude that ROI). RDMs are sorted (re-plotted) by dtheta as usual tuning (generalization) analysis.
        %
        %         1. group statistic (t-test or signtest against zero) of mean distance
        %
        %         2. R2 for fitting gaussian tuning function
        %         - small amplitude and large tuning width -> less likely to have modulation
        %         - R2<threshold -> no modulation (tuning)
        %         - somehow arbitrary selection (I admit)
        %
        %         3. test multiplicative/additive bimanual encoding (for ROIs that are tuned to both wrist movement)
        %         - get model predictions from fitted curve (fitting is done for only center of data, i.e., d(dtheta_L,0) and d(0,dtheta_R)).
        %         - get R2 gain for multiplicative and additive model predictions from null-model (only intercept)
        %         - test R2 gain against zero.
        %         - compare R2 gains between multiplicative and additive models.
        center = 1;
        vararginoptions(varargin,{'center'});
        
        % load data
        switch center
            case 0
                fname = fullfile(analyzeDir,'bmw1_tuning_summary_uncentered.mat');
            case 1
                fname = fullfile(analyzeDir,'bmw1_tuning_summary_centered.mat');
        end
        T = load(fname); %T=rmfield(T,'D');
        hemis = unique(T.hemis);
        regions = unique(T.region);
        
        E = [];
        for h=hemis'
            for reg=regions'
                D=getrow(T,T.hemis==h&T.region==reg);
                % 1. group statistic (t-test or signtest agrainst zero) of mean distance
                for i=1:3
                    [S.t_meandist(i),S.p_meandist(i)] = ttest_mc(ssqrt(D.meandist(:,i)),0,1,'onesample');
                end
                
                % 2. R2 for fitting gaussian tuning function
                waratio     = D.C./D.A; % width-amplitude ratio
                threshold   = 0.05e6;%200; % arbitrary value iffy
                encoding    = ones(size(D.R2));
                encoding(D.R2<0) = 0; % if negative R2, deem no encoding
                encoding(waratio>threshold) = 0; % small amplitude and large width is less likely to have encoding
                % get binomial probability
                S.binomP_encoding = 1-binocdf(sum(encoding,1),size(encoding,1),0.5);
                
                % 3. test multiplicative/additive bimanual encoding (for ROIs that are tuned to both wrist movement)
                alpha = 0.05;
                if (0)%(S.binomP_encoding(3)<alpha&&S.binomP_encoding(4)<alpha)
                    % better than null model?
                    for i=1:2
                        [S.t_R2gain(i),S.p_R2gain(i)] = ttest_mc(D.R2gain(:,i),0,1,'onesample');
                    end
                    
                    % which model is better?
                    [S.t_R2gain(3),S.p_R2gain(3)] = ttest_mc(D.R2gain(:,1),D.R2gain(:,2),1,'paired');
                else
                    S.t_R2gain = [NaN,NaN,NaN];
                    S.p_R2gain = S.t_R2gain;
                end
                S.hemis = h;
                S.region = reg;
                S.regname = {sprintf('%s-%s',regname{reg},hemName{h})};
                
                E = addstruct(E,S);
            end
        end
        % save result
        switch center
            case 0
                fname = fullfile(analyzeDir,'bmw1_ROI_encoding_summary_uncentered.mat');
            case 1
                fname = fullfile(analyzeDir,'bmw1_ROI_encoding_summary_centered.mat');
        end
        save(fname,'-struct','E');
        
        varargout = {E};
    case 'ROI_tuning_noiseceiling'                  % Get noise ceiling for phi (using correlation without subtracting intercept)
        center = 1;
        
        % load data
        switch center
            case 0
                fname = fullfile(analyzeDir,'bmw1_Gtuning_forplot_uncentered.mat');
            case 1
                fname = fullfile(analyzeDir,'bmw1_Gtuning_forplot_centered.mat');
        end
        T = load(fname);
        
        dtheta = unique(T.dtheta);
        subjects = unique(T.subj);
        hemis = unique(T.hemi);
        regions = unique(T.region);
        
        C = []; % C for ceiling
        % loop over hemisphere/region/subjects
        for h=hemis';
            for reg=regions';
                D = getrow(T,T.hemi==h&T.region==reg);
                
                % get overall mean of data
                ul_all = pivottablerow(D.dtheta,D.ul,'nanmean');
                ur_all = pivottablerow(D.dtheta,D.ur,'nanmean');
                bi_all = pivottablerow(D.dtheta,D.bi,'nanmean(x,1)');
                
                for s=subjects';
                    subset = D.subj~=s;
                    
                    % get leave-one-subject-out mean of data
                    ul_loo = pivottablerow(D.dtheta,D.ul,'nanmean','subset',subset); % leave one out
                    ur_loo = pivottablerow(D.dtheta,D.ur,'nanmean','subset',subset);
                    bi_loo = pivottablerow(D.dtheta,D.bi,'nanmean(x,1)','subset',subset);
                    
                    % left-out-data
                    ul = D.ul(~subset,:);
                    ur = D.ur(~subset,:);
                    bi = D.bi(~subset,:);
                    
                    % calc lower bound
                    c.ul(1,1) = corrN(ul,ul_loo);
                    c.ur(1,1) = corrN(ur,ur_loo);
                    c.bi(1,1) = corrN(bi(:),bi_loo(:));
                    
                    % calc upper bound
                    c.ul(1,2) = corrN(ul,ul_all);
                    c.ur(1,2) = corrN(ur,ur_all);
                    c.bi(1,2) = corrN(bi(:),bi_all(:));
                    
                    % other values
                    c.subj = s;
                    c.hemi = h;
                    c.region = reg;
                    
                    C = addstruct(C,c);
                end
            end
        end
        varargout = {C};
    case 'ROI_tuning_compare_model_baseline'     % Compare additive vs multiplicative models with different level of baseline subtraction
        arrange='axes';
        posAngle = [135 180 225 270 90 315 45 0];
        axsize = [0.2 0.14];
        fig = 0;
        center = 1;
        ex180 = 1;
        vararginoptions(varargin,{'fig','center','ex180'});
        
        % load data
        switch center
            case 0
                fname = fullfile(analyzeDir,'bmw1_Gtuning_forplot_uncentered.mat');
            case 1
                fname = fullfile(analyzeDir,'bmw1_Gtuning_forplot_centered.mat');
        end
        T = load(fname);
        
        subjects = unique(T.subj);
        hemis = unique(T.hemi);
        regions = unique(T.region);
        A = [1 0.3 0.1 0.03 0.01 0.003 0.001];
        
        
        R = [];
        for s=subjects'
            for h=hemis'
                for reg=regions'
                    T_ = getrow(T,T.subj==s&T.hemi==h&T.region==reg);
                    Phi = T_.bi;
                    dtheta = T_.dtheta;
                    Phi_ = mean(mean(T_.bi_));
                    
                    if ex180
                        Phi = Phi(dtheta~=180&dtheta~=-180,dtheta~=180&dtheta~=-180);
                        dtheta = dtheta(dtheta~=180&dtheta~=-180);
                    end
                    
                    % normalize Phi
                    Phi = Phi/Phi_;
                    
                    % TSS
                    Phi_tss = Phi;
                    tmp = vec(Phi_tss(dtheta~=0,dtheta~=0));
                    TSS = sum(tmp.*tmp);
                    
                    % Additive prediction (residual is not affected by baseline)
                    Phi_d0 = Phi(:,dtheta==0);
                    Phi_0d = Phi(dtheta==0,:);
                    Phi_00 = Phi(dtheta==0,dtheta==0);
                    Phi_true= (Phi);
                    Phi_pred = bsxfun(@plus,Phi_d0,Phi_0d)-Phi_00;
                    Phi_add = Phi_pred;
                    res_add = Phi_true-Phi_pred;
                    res = res_add(dtheta~=0,dtheta~=0);
                    
                    RSS_add = (sum(sum(res.*res)));
                    
                    R2_add = 1-RSS_add/TSS;
                    
                    % Multiplicative prediction (affected by baseline)
                    R2_multi = [];R_multi = [];RSS_multi=[];
                    % Loop over different levels of baseline subtraction
                    for a=A;
                        Phi_d0 = Phi(:,dtheta==0)+a;
                        Phi_0d = Phi(dtheta==0,:)+a;
                        Phi_00 = Phi(dtheta==0,dtheta==0)+a;
                        Phi_pred = bsxfun(@times,Phi_d0,Phi_0d)/Phi_00;
                        Phi_true = (Phi+a);
                        Phi_multi = Phi_pred-a; % adjust for plot
                        res_multi = Phi_true-Phi_pred;
                        res = res_multi(dtheta~=0,dtheta~=0);
                        
                        RSS_multi(end+1) = (sum(sum(res.*res)));
                        
                        R2_multi(end+1) = 1-RSS_multi(end)/TSS;
                        
                        % plot
                        if (fig)
                            figure('name',sprintf('%s-%s-%s-%d',subj_name{s},hemName{h},regname{reg},a));
                            
                            subplot(2,3,1);
                            mycontourf(dtheta,dtheta,Phi);axis square; title('true \Phi');
                            colorbar();
                            
                            subplot(2,3,2);
                            mycontourf(dtheta,dtheta,Phi_multi);axis square; title('predicted \Phi (multi)');
                            colorbar();
                            
                            subplot(2,3,3);
                            mycontourf(dtheta,dtheta,Phi_add);axis square; title('predicted \Phi (add)');
                            colorbar();
                            
                            subplot(2,3,4); % scatterplot
                            plot(vec(Phi_true(dtheta~=0,dtheta~=0)),vec(Phi_add(dtheta~=0,dtheta~=0)+a*Phi_),'b+');hold on
                            plot(vec(Phi_true(dtheta~=0,dtheta~=0)),vec(Phi_pred(dtheta~=0,dtheta~=0)),'r+');hold on
                            
                            subplot(2,3,5);
                            mycontourf(dtheta,dtheta,res_multi);axis square; title('residual (multi)');
                            colorbar();
                            
                            subplot(2,3,6);
                            mycontourf(dtheta,dtheta,res_add);axis square; title('residual (add)');
                            colorbar();
                            
                            pause();close();
                        end
                    end
                    
                    r.subj=repmat(s,length(A),1);
                    r.hemi=repmat(h,length(A),1);
                    r.region=repmat(reg,length(A),1);
                    r.A = A(end:-1:1)';
                    r.RSS_multi=RSS_multi';
                    r.RSS_add=repmat(RSS_add,length(A),1);
                    r.R2_multi=R2_multi';
                    r.R2_add=repmat(R2_add,length(A),1);
                    
                    R=addstruct(R,r);
                end
            end
        end
        
        % Plot summary (prediction R2 over baseline subtruction levels)
        for h=hemis'
            figure('name','R2 comparison',...
                'position',[10 10 25*sqrt(2),25]);
            %set(gcf,'DefaultPatchLineStyle','non');
            c=1;
            for reg=regions'
                
                subset = R.hemi==h&R.region==reg;
                
                switch (arrange)
                    case 'subplot'
                        subplot(ceil(length(regions)/2),2,c);
                    case 'axes'
                        radius = 0.5*(1-axsize(1))*0.9;
                        orig = [radius/0.9,radius/0.9]+[0.01 0.01];
                        pos = [orig(1)+radius*cosd(posAngle(reg)),orig(2)+radius*sind(posAngle(reg)),...
                            axsize(1),axsize(2)];
                        axes('position',pos);
                end
                
                plotval = [R.R2_multi,R.R2_add];
                %plotval = [-R.RSS_multi+R.RSS_add];
                lineplot(R.A,plotval,'subset',subset,...
                    'style_shade','linecolor',{[1,0.1,0.1],[0.1, 0.1, 1]},...
                    'shadecolor',{[1,0.8,0.8],[0.8, 0.8, 1]},...
                    'leg',{'Multi','Add'},'leglocation','best','transp',1);
                
                drawline(0,'dir','horz');
                
                xlabel('% removal of baseline');
                ylabel('Prediction R2');
                %ylabel('RSS (predicted)')
                %set(gca,'xscale','log','yscale','linear');
                title(sprintf('%s-%s',regname{reg},hemName{h}));
                
                c=c+1;
            end
        end
        
        % check if positive prediction R2
        P = [];
        for h=hemis'
            for reg=regions'
                subset = R.hemi==h&R.region==reg&R.A==1;
                multi = R.R2_multi(subset);
                addit = R.R2_add(subset);
                diffm = multi-addit;
                
                [~,O.pm] = ttest_mc(multi,0,1,'onesample');
                [~,O.pa] = ttest_mc(addit,0,1,'onesample');
                [~,O.pd] = ttest_mc(multi,addit,2,'paired');
                O.hemis = h;
                O.region = reg;
                P = addstruct(P,O);
            end
        end
        alpha = 0.05;
        figure('position',[5 5 10 30]);c=1;
        for h=hemis'
            for reg=regions'
                idx = P.hemis==h&P.region==reg;
                if (P.pm(idx)>alpha&&P.pa(idx)>alpha)&(~ismember(reg,[3,6,7]))
                else%if P.pd(idx)<alpha
                    
                    subset = R.hemi==h&R.region==reg&R.A==1;
                    subplot(ceil(length(regions)/2),2,c);
                    
                    plotval = [R.R2_multi,R.R2_add];
                    
                    [x,m,e]=barplot([],plotval,'subset',subset,...
                        'style_bold','facecolor',{red,blue});
                    drawline(0,'dir','horz');
                    
                    if P.pd(idx)<alpha
                        drawline(max(m+e)*1.1,'dir','horz','lim',x,'linewidth',2);
                        text(0.5*sum(x),max(m+e)*1.2,num2str(P.pd(idx)),'color','r');
                    else
                        text(0.5*sum(x),max(m+e)*1.2,num2str(P.pd(idx)));
                    end
                    
                    xlabel('');
                    ylabel('Prediction R^2');
                    set(gca,'ylim',[-0.1 1],'xtick',[],'tickdir','out',...
                        'ticklength',[0.03 0.01]);
                    title(sprintf('%s-%s',regname{reg},hemName_s{h}));
                    
                    c=c+1;
                end
            end
        end
        varargout = {R};
    case 'ROI_tuning_compare_model_alpha'        % Compare additive vs multiplicative models with different level of second movement influence
        arrange='axes';
        posAngle = [135 180 225 270 90 315 45 0];
        axsize = [0.2 0.14];
        fig = 0;
        center = 1;
        
        vararginoptions(varargin,{'fig','center'});
        
        % load data
        switch center
            case 0
                fname = fullfile(analyzeDir,'bmw1_Gtuning_forplot_uncentered.mat');
            case 1
                fname = fullfile(analyzeDir,'bmw1_Gtuning_forplot_centered.mat');
        end
        T = load(fname);
        
        subjects = unique(T.subj);
        hemis = unique(T.hemi);
        regions = unique(T.region);
        A = [0:0.05:1];
        A = A(1:end-1);
        
        R = [];
        for s=subjects'
            for h=hemis'
                for reg=regions'
                    T_ = getrow(T,T.subj==s&T.hemi==h&T.region==reg);
                    Phi = T_.bi;
                    dtheta = T_.dtheta;
                    
                    % normalize Phi with so Phi(0,0) = 1
                    Phi = Phi/(Phi(dtheta==0,dtheta==0));
                    
                    % TSS
                    tmp = vec(Phi(dtheta~=0,dtheta~=0));
                    TSS = sum(tmp.*tmp);
                    
                    % loop over different value of alpha (0<=alpha<1)
                    R2_multi = [];RSS_multi=[];
                    RSS_add=[]; R2_add = [];denom = [];
                    for a=A;
                        % reconstruct phi
                        Phi_shift = (circshift(Phi,[3 3])+circshift(Phi,[-3 3])...
                            +circshift(Phi,[3 -3])+circshift(Phi,[-3 -3]))/4; % 180 degree shifted version of original Phi
                        tmp = Phi-(2*a/(1+a^2))*Phi_shift;
                        Phi_recon = tmp/((1+a^2)-(4*a^2/(1+a^2)));
                        denom(end+1) = ((1+a^2)-(4*a^2/(1+a^2)));
                        
                        % again, normalize
                        Phi_recon = Phi_recon/Phi_recon(dtheta==0,dtheta==0);
                        
                        % take center parts
                        Phi_d0 = Phi_recon(:,dtheta==0);
                        Phi_0d = Phi_recon(dtheta==0,:);
                        Phi_00 = Phi_recon(dtheta==0,dtheta==0);
                        
                        % additive prediction
                        Phi_add = bsxfun(@plus,Phi_d0,Phi_0d)-Phi_00;
                        res_add = Phi_recon-Phi_add;
                        res     = res_add(dtheta~=0,dtheta~=0);
                        RSS_add(end+1) = (sum(sum(res.*res)));
                        R2_add(end+1)  = 1-RSS_add(end)/TSS; % could calc correlation
                        
                        % multiplicative prediction
                        Phi_multi = bsxfun(@times,Phi_d0,Phi_0d)/Phi_00;
                        res_multi = Phi_recon-Phi_multi;
                        res       = res_multi(dtheta~=0,dtheta~=0);
                        RSS_multi(end+1) = (sum(sum(res.*res)));
                        R2_multi(end+1)  = 1-RSS_multi(end)/TSS;
                        
                        % plot
                        if (fig)
                            figure('name',sprintf('%s-%s-%s-%d',subj_name{s},hemName{h},regname{reg},a));
                            plot(Phi_recon)
                            %                             subplot(2,3,1);
                            %                             mycontourf(dtheta,dtheta,Phi_recon);axis square; title('true \Phi (reconstructed)');
                            %                             colorbar();
                            %
                            %                             subplot(2,3,2);
                            %                             mycontourf(dtheta,dtheta,Phi_multi);axis square; title('predicted \Phi (multi)');
                            %                             colorbar();
                            %
                            %                             subplot(2,3,3);
                            %                             mycontourf(dtheta,dtheta,Phi_add);axis square; title('predicted \Phi (add)');
                            %                             colorbar();
                            %
                            %                             subplot(2,3,4); % scatterplot
                            %                             plot(vec(Phi_recon(dtheta~=0,dtheta~=0)),vec(Phi_add(dtheta~=0,dtheta~=0)),'b+');hold on
                            %                             plot(vec(Phi_recon(dtheta~=0,dtheta~=0)),vec(Phi_multi(dtheta~=0,dtheta~=0)),'r+');hold on
                            %
                            %                             subplot(2,3,5);
                            %                             mycontourf(dtheta,dtheta,res_multi);axis square; title('residual (multi)');
                            %                             colorbar();
                            %
                            %                             subplot(2,3,6);
                            %                             mycontourf(dtheta,dtheta,res_add);axis square; title('residual (add)');
                            %                             colorbar();
                            
                            pause();close();
                        end
                    end
                    
                    r.subj=repmat(s,length(A),1);
                    r.hemi=repmat(h,length(A),1);
                    r.region=repmat(reg,length(A),1);
                    r.A = A';
                    r.RSS_multi=RSS_multi';
                    r.RSS_add=RSS_add';
                    r.R2_multi=R2_multi';
                    r.R2_add=R2_add';
                    r.denom = denom';
                    R=addstruct(R,r);
                end
            end
        end
        
        % Plot summary (prediction R2 over baseline subtruction levels)
        for h=hemis'
            figure('name','R2 comparison',...
                'position',[10 10 25*sqrt(2),25]);
            %set(gcf,'DefaultPatchLineStyle','non');
            c=1;
            for reg=regions'
                
                subset = R.hemi==h&R.region==reg;
                
                switch (arrange)
                    case 'subplot'
                        subplot(ceil(length(regions)/2),2,c);
                    case 'axes'
                        radius = 0.5*(1-axsize(1))*0.9;
                        orig = [radius/0.9,radius/0.9]+[0.01 0.01];
                        pos = [orig(1)+radius*cosd(posAngle(reg)),orig(2)+radius*sind(posAngle(reg)),...
                            axsize(1),axsize(2)];
                        axes('position',pos);
                end
                
                plotval = [R.R2_multi,R.R2_add];
                %plotval = 1./R.denom;
                %plotval = [-R.RSS_multi+R.RSS_add]; % log likelihood ratio?
                lineplot(R.A,plotval,'subset',subset,...
                    'style_shade','linecolor',{[1,0.1,0.1],[0.1, 0.1, 1]},...
                    'shadecolor',{[1,0.8,0.8],[0.8, 0.8, 1]},...
                    'leg',{'Multi','Add'},'leglocation','best','transp',1);
                
                drawline(0,'dir','horz');
                
                xlabel('alpha (magnitude of second mvement influence)');
                ylabel('Prediction R2');
                set(gca,'xtick',A(1:5:end),'xlim',[0 1])
                %ylabel('RSS (predicted)')
                %set(gca,'xscale','log','yscale','linear');
                title(sprintf('%s-%s',regname{reg},hemName{h}));
                
                c=c+1;
            end
        end
        
        varargout = {R};
        
    case 'Beh_tuning_get'                           % Calculate RDM tuning from behavioural data
        direction = [0:60:300]';
        
        glm = 10;
        center = 1;
        %vararginoptions(varargin,{'glm','center'});
        
        % load data
        T = load(fullfile(analyzeDir,'bmw1_behRDMIPM_all.mat'));
        
        subjects    = unique(T.SN);
        
        Tuning      = []; Tu = [];
        for s = subjects'
            D = getrow(T,T.SN==s);
            
            G = rsa_squareIPM(D.IPM);
            
            % symmetrize G
            G = 0.5*(G+G');
            
            % split IPM into unimanual and bimnaual parts (UL,UR,Bi)
            UL = G(1:nDirections,1:nDirections);
            UR = G(nDirections+1:2*nDirections,nDirections+1:2*nDirections);
            Bi = G(2*nDirections+1:end,2*nDirections+1:end);
            
            % center
            if center
                [UL,UL_] = doubleCenter(UL);
                [UR,UR_] = doubleCenter(UR);
                [Bi,Bi_] = doubleCenter(Bi);
            else
                UL_ = zeros(size(UL));
                UR_ = zeros(size(UR));
                Bi_ = zeros(size(Bi));
            end
            
            % sort IPM by dtheta
            [ul,diffs]  = RDM2Tuning(UL,'uni',1);
            ur          = RDM2Tuning(UR,'uni',1);
            bi          = RDM2Tuning(Bi,'bi',1);
            ul_         = RDM2Tuning(UL_,'uni',1);
            ur_         = RDM2Tuning(UR_,'uni',1);
            bi_         = RDM2Tuning(Bi_,'bi',1);
            
            % sort RDM by dtheta
            
            % get mean dist
            
            
            % tuning data for plot
            tuning.ul       = ul;
            tuning.ur       = ur;
            tuning.bi       = bi;
            tuning.ul_      = ul_;
            tuning.ur_      = ur_;
            tuning.bi_      = bi_;
            tuning.dtheta   = diffs; ndiff = length(diffs);
            tuning.subj     = repmat(s,ndiff,1);
            tuning.region   = ones(ndiff,1);
            tuning.hemi     = ones(ndiff,1);
            tuning.center   = repmat(center,ndiff,1);
            Tuning          = addstruct(Tuning,tuning);
        end
        % save result
        switch center
            case 0
                fname = fullfile(analyzeDir,'bmw1_Beh_tuning_forplot_uncentered.mat');
            case 1
                fname = fullfile(analyzeDir,'bmw1_Beh_tuning_forplot_centered.mat');
        end
        
        save(fname,'-struct','Tuning');
        
        % plot
        bmw1_imana('ROI_tuning_plot',Tuning,'Ylim_uni',[-1.5*10e-1 3*10e-1],...
            'Ylim_bi',[-1.5*10e-1 3*10e-1],'arrange','axes','mode','IPM');
        
        varargout = {Tuning,Tu,T};
        
    case 'simulate_RDM'
        directions = [0:60:300];
        
        [modelRDM,U,Bm,Ba] = getRDMTuning(directions);
        
        % plot
        figure('name','simulated RDM and tuning');
        subplot(2,3,1);
        imagesc(squareform(modelRDM.RDM{1}));axis square
        subplot(2,3,2);
        imagesc(squareform(modelRDM.RDM{2}));axis square;title('multi')
        subplot(2,3,3);
        imagesc(squareform(modelRDM.RDM{3}));axis square;title('add')
        
        subplot(2,3,4);
        plot(U.dtheta,U.tuning);axis square
        subplot(2,3,5);
        contourf(Bm.dtheta,Bm.dtheta,Bm.tuning);axis square;title('multi')
        subplot(2,3,6);
        contourf(Ba.dtheta,Ba.dtheta,Ba.tuning);axis square;title('add')
        
        varargout = {};
    otherwise
        disp('there is no such case.')
end;

% return to original
set(0,'DefaultFigureUnit',DefaultFigureUnit);

%% Local functions
    function dircheck(dir)
        % Checks existance of specified directory. Makes it if it does not exist.
        % SArbuckle 01/2016
        if ~exist(dir,'dir');
            %warning('%s didn''t exist, so this directory was created.\n',dir);
            mkdir(dir);
        end
    end

    function mysetEnv()
        disp('setting environmental variables...')
        %--- Set environmental variables for freesurfer and FSL
        path1 = getenv('PATH');
        add = [':/Applications/freesurfer/bin'...
            ':/Applications/freesurfer/fsfast/bin'...
            ':/Applications/freesurfer/mni/bin'...
            ':/Applications/caret/bin_macosx32'...
            ':/usr/local/fsl/bin'...
            ':/usr/bin'... % not exist? was ':/usr/local/bin'
            ':/usr/X11/bin'...
            ':/usr/X11R6/bin'];
        if isempty(strfind(path1,add))
            path1 = [path1, add];
            setenv('PATH', path1);
        end
        % Environment variables for freesurfer
        command_freesurfer1 = 'source $FREESURFER_HOME/SetUpFreeSurfer.sh';
        command_freesurfer2 = 'source $FREESURFER_HOME/FreeSurferEnv.sh';
        setenv('FREESURFER_HOME','/Applications/freesurfer');
        %         [~,result]=system('export FREESURFER_HOME=/Applications/freesurfer')
        [~,result]=system(command_freesurfer1);
        [~,result]=system(command_freesurfer2);
        getenv('SUBJECTS_DIR')
        getenv('MNI_DIR')
        setenv('SUBJECTS_DIR', freesurferDir);%'/Users/Atsushi/MotorControl/Atlas_templates/suit_flat');
        setenv('MNI_DIR','/Applications/freesurfer/mni');
        
        % Environment variables for FSL
        setenv('FSLDIR','/usr/local/fsl');
        setenv('FSLOUTPUTTYPE','NIFTI');
        setenv('FSLMULTIFILEQUIT','TRUE');
        setenv('FSLTCLSH','/usr/local/fsl/bin/fsltclsh');
        setenv('FSLWISH','/usr/local/fsl/bin/fslwish');
        setenv('FSLCONFDIR','/usr/local/fsl/config');
        setenv('FSLMACHTYPE','/usr/local/fsl/etc/fslconf/fslmachtype.sh');
        
        out = 1;
    end

    function setenvLinux()
        disp('setting environmental variables in Linux...')
        %--- Set environmental variables for freesurfer and FSL
        path1 = getenv('PATH');
        add = [':/usr/local/freesurfer/bin'...
            ':/usr/local/freesurfer/fsfast/bin'...
            ':/usr/local/freesurfer/mni/bin'...
            ':/etc/caret/bin_linux64/'...
            ':/usr/share/fsl/5.0/bin'];
        if isempty(strfind(path1,add))
            path1 = [path1, add];
            setenv('PATH', path1);
        end
        setenv('FREESURFER_HOME','/usr/local/freesurfer');
        % Environment variables for freesurfer
        command_freesurfer1 = 'source $FREESURFER_HOME/SetUpFreeSurfer.sh';
        command_freesurfer2 = 'source $FREESURFER_HOME/FreeSurferEnv.sh';
        %         [~,result]=system('export FREESURFER_HOME=/Applications/freesurfer')
        [~,result]=system(command_freesurfer1);
        [~,result]=system(command_freesurfer2);
        [~,result]=system('source /home/dduarte/.bashrc');
        getenv('SUBJECTS_DIR')
        getenv('MNI_DIR')
        setenv('SUBJECTS_DIR', freesurferDir);%'/Users/Atsushi/MotorControl/Atlas_templates/suit_flat');
        setenv('MNI_DIR','/usr/local/freesurfer/mni');
        
        % Environment variables for FSL
        setenv('FSLDIR','/usr/share/fsl/5.0/');
        setenv('FSLOUTPUTTYPE','NIFTI');
        %         setenv('FSLMULTIFILEQUIT','TRUE');
        %         setenv('FSLTCLSH','/usr/share/fsl/5.0/bin/fsltclsh');
        %         setenv('FSLWISH','/usr/share/fsl/5.0/bin/fslwish');
        %         setenv('FSLCONFDIR','/usr/local/fsl/config');
        %         setenv('FSLMACHTYPE','/usr/local/fsl/etc/fslconf/fslmachtype.sh');
        
        out = 1;
    end

    function [RDMvec_out] = splitUB(RDMvec_in,part)
        [nRow,nRDM] = size(RDMvec_in);
        
        % make index for unimanual left, unimanual right and
        % bimanual trials
        idxU = tril(~logical(eye(nDirections)));
        idxB = tril(~logical(eye(nDirections^2)));
        switch (part)
            case 'All'
                idxAll = tril(~logical(eye(nDirections*(nDirections+2))));
            case 'Uni'
                idxAll = blockdiag(idxU,idxU,false(nDirections^2));
            case 'UL'
                idxAll = blockdiag(idxU,false(nDirections),false(nDirections^2));
            case 'UR'
                idxAll = blockdiag(false(nDirections),idxU,false(nDirections^2));
            case 'Bi'
                idxAll = blockdiag(false(nDirections),false(nDirections),idxB);
        end
        idxAll = logical(idxAll);
        for row=1:nRow
            rdm_tmp = squareform(RDMvec_in(row,:));
            RDMvec_out(row,:) = rdm_tmp(idxAll);
        end
    end

    function [out,U,Bm,Ba] = getRDMTuning(directions)
        % get patterns
        phi = [0:359]; % preferred directions (uniform)
        tuningwidth = 30;
        amplitude = 0.99;
        for theta=1:length(directions)
            % uni
            R_uni(theta,:) = standardGaussTuning(directions(theta),phi,tuningwidth,amplitude);
            % bi
            for theta2=1:length(directions)
                R_uni2 = standardGaussTuning(directions(theta2),phi,tuningwidth,amplitude);
                % multiplicative
                R_bi_mul(theta,theta2,:) = R_uni(theta,:).*R_uni2;
                % additive
                R_bi_add(theta,theta2,:) = R_uni(theta,:)+R_uni2;
            end
        end
        R_bi_mul = reshape(R_bi_mul,nDirections^2,length(phi));
        R_bi_add = reshape(R_bi_add,nDirections^2,length(phi));
        
        % get RDM
        con_uni = [1:length(directions)];
        con_bi  = [1:nDirections^2];
        C_uni   = indicatorMatrix('allpairs',con_uni);
        C_bi    = indicatorMatrix('allpairs',con_bi);
        
        G_uni       = R_uni*R_uni';
        G_bi_mul    = R_bi_mul*R_bi_mul';
        G_bi_add    = R_bi_add*R_bi_add';
        
        dist_uni = diag(C_uni*G_uni*C_uni');
        dist_mul = diag(C_bi*G_bi_mul*C_bi');
        dist_add = diag(C_bi*G_bi_add*C_bi');
        
        dist_uni = dist_uni'/norm(dist_uni);
        dist_mul = dist_mul'/norm(dist_mul);
        dist_add = dist_add'/norm(dist_add);
        
        out.RDM{1,1} = dist_uni;
        out.RDM{2,1} = dist_mul;
        out.RDM{3,1} = dist_add;
        out.name = {'unimanual';'bimanual-multiplicative';'bimanual-additive'};
        out.subj = [0;0;0];
        out.region = [0;0;0];
        out.hemis = [1;1;1];
        
        [U.tuning,U.dtheta]   = RDM2Tuning(dist_uni,'uni');
        [Bm.tuning,Bm.dtheta]  = RDM2Tuning(dist_mul,'bi');
        [Ba.tuning,Ba.dtheta]  = RDM2Tuning(dist_add,'bi');
        
    end

    function r = standardGaussTuning(theta,phi,sigma,amplitude)
        delta = acosd(cosd(phi-theta));
        r = amplitude.*exp(-0.5*(delta).*(delta)/(sigma*sigma))+(1-amplitude);
    end

    function [out,diffs] = RDM2Tuning(in,type,is180)
        [row,col] = size(in);
        
        if row~=col
            % make square form
            in=squareform(in);
        end
        
        % get delta-direction indices for unimanual trials
        direction = [0:60:300];
        for d=1:nDirections
            Uni(:,d) = direction-direction(d);
        end
        signs = sign(sind(Uni));
        signs(signs==0)=1;
        Uni = acosd(cosd(Uni));
        Uni = Uni.*signs;
        % for bimanual trials
        BiL = kron(Uni,ones(nDirections));
        BiR = kron(ones(nDirections),Uni);
        
        diffs       = unique(Uni);
        switch (type)
            case 'uni'
                % sort distances
                for d=1:length(diffs)
                    out(d,1) = nanmean(in(Uni==diffs(d)));
                end
                if is180
                    diffs = [-diffs(end);diffs];
                    out = [out(end);out];
                end
            case 'bi'
                % sort distances
                for d=1:length(diffs)
                    for dd=1:length(diffs)
                        out(d,dd) = nanmean(in(BiL==diffs(d)&BiR==diffs(dd)));
                    end
                end
                if is180
                    diffs = [-diffs(end);diffs];
                    out = [out(end,:);out];
                    out = [out(:,end),out];
                end
        end
    end

    function [O,M] = doubleCenter(I)
        % double center Matrix to force it to have zero mean
        n = length(I);
        one = ones(length(I),1);
        Me = one*one'/n;
        Cen = eye(size(I))-Me;
        O = Cen*I*Cen';
        M = Me*I*Me';
    end

    function scatterplotMDS(varargin)
        MDS     = varargin{1};
        split   = varargin{2};
        label   = varargin{3};
        
        % Plot properties
        color               = mat2cell(hsv(6),ones(6,1),3);
        CAT.markercolor     = color;
        CAT.markerfill      = color;
        CAT.markertype      = 'o';
        CAT.markersize      = 7;
        
        % Plot
        scatterplot3(MDS(1:22,1),MDS(1:22,2),MDS(1:22,3),'split',split,...
            'label',label,'CAT',CAT); hold on;
        
        %         for i=1:2 % drawing line between unimanual conditions
        %             hold on;
        %             indx=[1:nDirections 1]'+(i-1)*nDirections;
        %             line(MDS(indx,1),MDS(indx,2),MDS(indx,3),'color','k');
        %         end;
        %
        %         % drawing line between bimanual conditions
        %         indx=[2*nDirections+1:nConditions,2*nDirections+1]';
        %         line(MDS(indx,1),MDS(indx,2),MDS(indx,3),'color','k');
        %         hold on;
        
        % plotting origin
        plot3(0,0,0,'kpentagram','markersize',3*CAT.markersize,...
            'markerfacecolor','k');
        hold off;
        axis equal;
    end



end

function makefieldmap(dataDir, subj_name, run, varargin)
prefixepi='a';
prefixfieldmap = '';
image=1;
subfolderRawdata='';
subfolderFieldmap='';
rawdataDir='';

use3D=false;
writeunwarped = 1;

vararginoptions(varargin,{'prefixepi','prefixfieldmap', 'image',...
    'subfolderRawdata', 'subfolderFieldmap', 'use3D','rawdataDir','writeunwarped'});


spm_dir= fileparts(which('spm'));
spmVer=spm('Ver');

% displaying whether using 3D files or not
% disp(['Using 3D: ' num2str(use3D)])

%_______DEFAULTS_________________________________
J.defaults.defaultsval.et = [4.92 7.38]; % gre_field_maping at UWO RRI Prisma 3T
J.defaults.defaultsval.maskbrain = 1;
J.defaults.defaultsval.blipdir = -1;
J.defaults.defaultsval.tert = 90*0.7/2;% (=[matrix size in phase-encoding direction] X [echo spacing time (in ms) / [iPAT factor]])
J.defaults.defaultsval.epifm = 0;
J.defaults.defaultsval.ajm = 0;
J.defaults.defaultsval.uflags.method = 'Mark3D';
J.defaults.defaultsval.uflags.fwhm = 10;
J.defaults.defaultsval.uflags.pad = 0;
J.defaults.defaultsval.uflags.ws = 1;
switch (spmVer)
    case 'SPM12'
        J.defaults.defaultsval.mflags.template = {fullfile(spm_dir,'canonical','avg152T1.nii')};
    case 'SPM8'
        J.defaults.defaultsval.mflags.template = {fullfile(spm_dir,'templates','T1.nii')};
end;
J.defaults.defaultsval.mflags.fwhm = 5;
J.defaults.defaultsval.mflags.nerode = 2;
J.defaults.defaultsval.mflags.ndilate = 4;
J.defaults.defaultsval.mflags.thresh = 0.5;
J.defaults.defaultsval.mflags.reg = 0.02;
J.matchvdm = 1;
J.writeunwarped = writeunwarped;
J.anat = [];
J.matchanat = 0;


if (isempty(rawdataDir))
    rawdataDir=fullfile(dataDir, 'imaging_data_raw', subj_name, subfolderRawdata);
end;

%_______for multiple run with same fieldmap_________________________________
for i=1:numel(run)
    if use3D
        J.session(i).epi ={fullfile(rawdataDir, [prefixepi subj_name,'_run_',run{i},'_',num2str(image),'.nii'])};
    else
        J.session(i).epi ={fullfile(rawdataDir, [prefixepi,subj_name,'_run_',run{i},'.nii,',num2str(image)])};
    end;
end
J.sessname = {'run_RENAMErunNumber'};%{['run_', run{i}]};
%_______change code here if we have 1 fieldmap for each run_________________________________
J.phase ={fullfile(dataDir, 'fieldmaps', subj_name, subfolderFieldmap, [prefixfieldmap,subj_name,'_phase.nii,1'])};          %,'_',num2str(run(1))
J.magnitude =  {fullfile(dataDir, 'fieldmaps', subj_name, subfolderFieldmap, [prefixfieldmap,subj_name,'_magnitude.nii,1'])};  %,'_',num2str(run(1))

matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj= J;

spm_jobman('run',matlabbatch);
end
function realign_unwarp(dataDir, subj_name, run, startTR, endTR, varargin)
% spmj_realign_unwarp(dataDir, subj_name, run, startTR, endTR)
% INPUT:
%   dataDir:    Root directory for the imaging structure (needs directories
%               imaging_data_raw and fieldmaos
%   subj_name:  Name for subdirectory and file-prefix (i.e. 's05')
%   run:        Cell array of identifiers for the run
%               {'01' '02','03','04','05','06','07','08'}
%   startTR:    First image to align
%   endTR:      Last image to align (if INF, it will use all available)
% VARARGINOPTIONS:
%   'prefix'            prefix for run name (default 'a');
%   'scanType'          sub folder in your subject directory
%   'subfolderFieldmap' subfolder in the fieldmap directory
%   'subfolderRawdata'  subfolder in the imaging_data_raw directory
%   'rawdataDir'        overwrites standard naming and forces routine to
%                       use this folder for location of raw data
% Tobias Wiestler & Joern Diedrichsen
% 06/02/2012 subfolder option replaced with two options 'subfolderFieldmap' 'subfolderRawdata'
% 23/10/2012 added rawdataDir to be able to overwrite the standard naming convention
%
prefix= 'a';
prefixfieldmap = 'r';
subfolderRawdata='';
subfolderFieldmap='';
use3D=false;
rawdataDir='';
targetImg = [];
vararginoptions(varargin,{'prefix', 'prefixfieldmap','subfolderRawdata', 'subfolderFieldmap','use3D','rawdataDir','targetImg'});

% displaying whether using 3D files or not
disp(['Using 3D: ' num2str(use3D)]);


%_______DEFAULTS_________________________________
J.eoptions.quality = 0.9;
J.eoptions.sep = 2;%4;
J.eoptions.fwhm = 5;
J.eoptions.rtm = 0; %why zero and not one
J.eoptions.einterp = 2;
J.eoptions.ewrap = [0 1 0];     %  wrap-around in the [x y z] direction during the estimation (ewrap)  wrap of the front of the head to the back of the head
J.eoptions.weight = {''};
J.uweoptions.basfcn = [12 12];
J.uweoptions.regorder = 1;
J.uweoptions.lambda = 100000;
J.uweoptions.jm = 0;
J.uweoptions.fot = [4 5];
J.uweoptions.sot = [1];
J.uweoptions.uwfwhm = 4;
J.uweoptions.rem = 1;
J.uweoptions.noi = 5;
J.uweoptions.expround = 'Average';
J.uwroptions.uwwhich = [2 1]; %[2 1] with mean image without [2 0]
J.uwroptions.rinterp = 4;
J.uwroptions.wrap = [0 1 0];  %  wrap-around in the [x y z] direction during the reslicing (wrap)
J.uwroptions.mask = 1;
J.uwroptions.prefix = 'u';


if (isempty(rawdataDir))
    rawdataDir=fullfile(dataDir, 'imaging_data_raw',subj_name,subfolderRawdata);
end;

%_______images and fieldmap definition_________________________________
for j=1:numel(run)
    if (isinf(endTR))  % All avaialble images: only works with 4d-nifits right now
        V = nifti(fullfile(rawdataDir,[prefix,subj_name,'_run_',run{j},'.nii']));
        imageNumber=startTR:V.dat.dim(4);
    else
        imageNumber= startTR:endTR;
    end;
    for i= 1:numel(imageNumber)
        if use3D
            scans{i,1}= fullfile(rawdataDir, [prefix subj_name,'_run_',run{j},'_',num2str(imageNumber(i)),'.nii']);
        else
            scans{i,1}= fullfile(rawdataDir, [prefix,subj_name,'_run_',run{j},'.nii,',num2str(imageNumber(i))]);
        end;
        
    end
    J.data(j).scans = scans;
    J.data(j).pmscan = {fullfile(dataDir, 'fieldmaps',subj_name,subfolderFieldmap,['vdm5_sc',prefixfieldmap,subj_name,'_phase_session',num2str(j),'.nii,1'])};
end

% If you want to align images to non-first volume of the first session
% (e.g., last volume of the last session)
if ~isempty(targetImg)
    J.data(1).scans = cat(1, {targetImg}, J.data(1).scans); % add copy of the specified image to top
end

matlabbatch{1}.spm.spatial.realignunwarp= J;
spm_jobman('run',matlabbatch);

% Remove the (dummy) first image
% if ~isempty(targetImg)
%     j=1;
%     if use3D % rename
%         c=0;
%         for i=2:numel(imageNumber)+1
%             c=c+1;
%             source = fullfile(rawdataDir, [J.uwroptions.prefix, prefix subj_name,'_run_',run{j},'_',num2str(i),'.nii']);
%             dest = fullfile(rawdataDir, [J.uwroptions.prefix, prefix subj_name,'_run_',run{j},'_',num2str(c),'.nii']);
%             copyfile(source,dest);
%         end;
%         removefile(source); % delete last one
%         
%     else % create new 4D file
%         c=0;
%         for i=2:numel(imageNumber)+1
%             c=c+1;
%             P{c,1} = fullfile(rawdataDir, [J.uwroptions.prefix, prefix,subj_name,'_run_',run{j},'.nii,',num2str(i)]);
%         end
%         outfilename = fullfile(rawdataDir, [J.uwroptions.prefix, prefix,subj_name,'_run_',run{j},'.nii']);
%         spm_file_merge(char(P),outfilename);
%     end    
% end

end
function realign(dataDir, subj_name, run, startTR, endTR, varargin)
% spmj_realign_unwarp(dataDir, subj_name, run, startTR, endTR)
% INPUT:
%   dataDir:    Root directory for the imaging structure (needs directories
%               imaging_data_raw and fieldmaos
%   subj_name:  Name for subdirectory and file-prefix (i.e. 's05')
%   run:        Cell array of identifiers for the run
%               {'01' '02','03','04','05','06','07','08'}
%   startTR:    First image to align
%   endTR:      Last image to align (if INF, it will use all available)
% VARARGINOPTIONS:
%   'prefix'            prefix for run name (default 'a');
%   'scanType'          sub folder in your subject directory
%   'subfolderFieldmap' subfolder in the fieldmap directory
%   'subfolderRawdata'  subfolder in the imaging_data_raw directory
%   'rawdataDir'        overwrites standard naming and forces routine to
%                       use this folder for location of raw data
% Tobias Wiestler & Joern Diedrichsen
% 06/02/2012 subfolder option replaced with two options 'subfolderFieldmap' 'subfolderRawdata'
% 23/10/2012 added rawdataDir to be able to overwrite the standard naming convention
%
prefix= 'a';
subfolderRawdata='';
use3D=false;
rawdataDir='';
targetImg = [];
vararginoptions(varargin,{'prefix','subfolderRawdata','use3D','rawdataDir','targetImg'});

% displaying whether using 3D files or not
disp(['Using 3D: ' num2str(use3D)]);


%_______DEFAULTS_________________________________
% J.eoptions.quality = 0.9;
% J.eoptions.sep = 2;%4;
% J.eoptions.fwhm = 5;
% J.eoptions.rtm = 0; %why zero and not one
% J.eoptions.einterp = 2;
% J.eoptions.ewrap = [0 1 0];     %  wrap-around in the [x y z] direction during the estimation (ewrap)  wrap of the front of the head to the back of the head
% J.eoptions.weight = {''};

J.eoptions.quality = 0.9;
J.eoptions.sep = 4;
J.eoptions.fwhm = 5;
J.eoptions.rtm = 1;
J.eoptions.interp = 2;
J.eoptions.wrap = [0 0 0];
J.eoptions.weight = {''};
J.roptions.which = [2 1];
J.roptions.interp = 4;
J.roptions.wrap = [0 0 0];
J.roptions.mask = 1;
J.roptions.prefix = 'r';

if (isempty(rawdataDir))
    rawdataDir=fullfile(dataDir, 'imaging_data_raw',subj_name,subfolderRawdata);
end;

%_______images definition_________________________________
for j=1:numel(run)
    if (isinf(endTR))  % All avaialble images: only works with 4d-nifits right now
        V = nifti(fullfile(rawdataDir,[prefix,subj_name,'_run_',run{j},'.nii']));
        imageNumber=startTR:V.dat.dim(4);
    else
        imageNumber= startTR:endTR;
    end;
    for i= 1:numel(imageNumber)
        if use3D
            scans{i,1}= fullfile(rawdataDir, [prefix subj_name,'_run_',run{j},'_',num2str(imageNumber(i)),'.nii']);
        else
            scans{i,1}= fullfile(rawdataDir, [prefix,subj_name,'_run_',run{j},'.nii,',num2str(imageNumber(i))]);
        end;
        
    end
    J.data{j,1} = scans;
    %J.data(j).pmscan = {fullfile(dataDir, 'fieldmaps',subj_name,subfolderFieldmap,['vdm5_sc',prefixfieldmap,subj_name,'_phase_session',num2str(j),'.nii,1'])};
end

% If you want to align images to non-first volume of the first session
% (e.g., last volume of the last session)
if ~isempty(targetImg)
    J.data{1} = cat(1, {targetImg}, J.data{1}); % add copy of the specified image to top
end

matlabbatch{1}.spm.spatial.realign.estwrite= J;
spm_jobman('run',matlabbatch);

% Remove the (dummy) first image
% if ~isempty(targetImg)
%     j=1;
%     if use3D % rename
%         c=0;
%         for i=2:numel(imageNumber)+1
%             c=c+1;
%             source = fullfile(rawdataDir, [J.uwroptions.prefix, prefix subj_name,'_run_',run{j},'_',num2str(i),'.nii']);
%             dest = fullfile(rawdataDir, [J.uwroptions.prefix, prefix subj_name,'_run_',run{j},'_',num2str(c),'.nii']);
%             copyfile(source,dest);
%         end;
%         removefile(source); % delete last one
%         
%     else % create new 4D file
%         c=0;
%         for i=2:numel(imageNumber)+1
%             c=c+1;
%             P{c,1} = fullfile(rawdataDir, [J.uwroptions.prefix, prefix,subj_name,'_run_',run{j},'.nii,',num2str(i)]);
%         end
%         outfilename = fullfile(rawdataDir, [J.uwroptions.prefix, prefix,subj_name,'_run_',run{j},'.nii']);
%         spm_file_merge(char(P),outfilename);
%     end    
% end

end