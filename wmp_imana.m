function varargout=wmp_imana(what,varargin)
% Script for minimal preprocessing of the wmp3/wmp4 data

%% ========================================================================================================================
% SCANNER INFO
numDummys = 3;                                                             % per run
numTRs    = 390;                                                           % per run (includes dummies)
%% ========================================================================================================================
% PATH DEFINITIONS
baseDir  = '/srv/diedrichsen/data/SequenceAndChord/WorkingMemoryPlanning/wmp3';
if ~exist(baseDir,'dir')
    baseDir  ='/Volumes/diedrichsen_data$/data/SequenceAndChord/WorkingMemoryPlanning/wmp3';
end
imagingDir      ='/imaging_data';
imagingDirRaw   ='/imaging_data_raw';
anatomicalDir   ='/anatomicals';
regDir          ='/RegionOfInterest';
fmapDir         ='/fieldmaps';
freesurferDir   ='/surfaceFreesurfer';
%% ========================================================================================================================
% PRE-PROCESSING 
subj_name = {'p01'};
loc_AC={[79;123;132]}; % Coordinates of anterior commissure in mm.  Use SPM Display.
%% ========================================================================================================================
% GLM INFO
run         = {'01','02','03','04','05','06','07','08','09','10'};
runB        = [1,2,3,4,5,6,7,8,9,10];  % Behavioural labelling of run
%% ========================================================================================================================
% Global default values for plots throughout
% colors
cbs_red = [213 94 0]/255;
cbs_blue = [0 114 178]/255;
cbs_yellow = [240 228 66]/255;
cbs_pink = [204 121 167]/255;
cbs_green = [0 158 115]/255;
blue = [49,130,189]/255;
lightblue = [158,202,225]/255;
red = [222,45,38]/255;
lightred = [252,146,114]/255;
green = [49,163,84]/255;
lightgreen = [161,217,155]/255;
orange = [253,141,60]/255;
yellow = [254,196,79]/255;
lightyellow = [255,237,160]/255;
purple = [117,107,177]/255;
lightpurple = [188,189,220]/255;
darkgray = [50,50,50]/255;
gray = [150,150,150]/255;
lightgray = [200,200,200]/255;
silver = [240,240,240]/255;
black = [0,0,0]/255;
pastelred = [255 105 97]/255;
pastelblue = [97 105 255]/255;

% plot defaults
fs = 20; % default font size for all figures
lw = 4;  % default line width for all figures
ms = 12; % default marker size for all figures

% styles
style.reset;
style.custom({blue,lightblue,red,lightred,orange,yellow,lightyellow,purple,lightpurple,darkgray,gray,lightgray,green,lightgreen,black,silver,...
    cbs_red,cbs_yellow,cbs_blue,cbs_green,cbs_pink,pastelred,pastelblue});
baselinesty_dotted = style.custom({black},'linestyle',':');
baselinesty_solid = style.custom({black},'linestyle','-');
lightgraysty = style.custom({lightgray}, 'markertype', 'none', 'linewidth', lw, 'errorwidth', lw, 'errorcap', 0, 'linestyle', '--');
%% ========================================================================================================================
switch(what)
    
    case 'ANAT:reslice_LPI'               % Reslice anatomical image within LPI coordinate systems
        % example: wmp_imana('ANAT:reslice_LPI',1)
        sn  = varargin{1}; % subjNum 
        subjs=length(sn);
        
        for s=1:subjs
            
            % (1) Reslice anatomical image to set it within LPI co-ordinate frames
            source  = fullfile(baseDir,anatomicalDir,subj_name{sn(s)},'anatomical_raw.nii');
            dest    = fullfile(baseDir,anatomicalDir,subj_name{sn(s)},'anatomical.nii');
            spmj_reslice_LPI(source,'name', dest);
            
            % (2) In the resliced image, set translation to zero
            V               = spm_vol(dest);
            dat             = spm_read_vols(V);
            V.mat(1:3,4)    = [0 0 0];
            spm_write_vol(V,dat);
            fprintf(1,'Manually retrieve the location of the anterior commissure (x,y,z) before continuing.\n');
        end
    case 'ANAT:center_AC'                 % Re-centre AC
        % Before running provide coordinates in the preprocessing section
        % 1. Use SPM Display to locate anterior commissure
        % 2. Enter coordinates to the loc_AC varible in
        %       PRE-PROCESSING section at top of this file
        % 3. AC coordinate sets should be listed in the same order
        %       as subjects in subj_name.
        % example: wmp_imana('ANAT:center_AC',1)
        sn=varargin{1}; % subjNum
        
        subjs=length(sn);
        for s=1:subjs
            img    = fullfile(baseDir,anatomicalDir,subj_name{sn(s)},'anatomical.nii');
            V               = spm_vol(img);
            dat             = spm_read_vols(V);
            oldOrig         = V.mat(1:3,4);
            V.mat(1:3,4)    = oldOrig-loc_AC{sn(s)};
            spm_write_vol(V,dat);
            fprintf('Done for %s \n',subj_name{sn(s)})
        end
    case 'ANAT:segmentation'              % Segmentation + Normalisation
        % example: wmp_imana('ANAT:segmentation',1)
        sn=varargin{1}; % subjNum
        
        subjs=length(sn);
        
        SPMhome=fileparts(which('spm.m'));
        J=[];
        for s=1:subjs
            
            J.channel.vols = {fullfile(baseDir,anatomicalDir,subj_name{sn(s)},'anatomical.nii,1')};
            J.channel.biasreg = 0.001;
            J.channel.biasfwhm = 60;
            J.channel.write = [0 1];
            
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
            J.warp.write = [1 1];
            
            matlabbatch{1}.spm.spatial.preproc=J;
            spm_jobman('run',matlabbatch);
            fprintf('Check segmentation results for %s\n', subj_name{sn(s)})
        end
    
    case 'FMAP:makefieldmap'
        % read this one
        % https://lcni.uoregon.edu/kb-articles/kb-0003
        sn  = varargin{1}; % subjNum 
        subjs=length(sn);
        
        
        et1 = 4.08; % use the magnitude of this file because has better quality
        et2 = 5.10;
        % total echo readout time
        tert = 23.0033;

        for s=1:subjs
            spmj_makefieldmap(baseDir,subj_name{sn(s)},run,'prefix','','et1',et1,'et2',et2,'tert',tert);
            fprintf('Done for %s \n',subj_name{sn(s)})
        end
        
    case 'FUNC:remDum'                    % Remove the extra dummy scans from all functional runs.
        % funtional runs have to be named as subjname_run01-r.nii, subjname_run02-r.nii ...
        % example: wmp_imana('FUNC:remDum',1)
        sn=varargin{1}; % subjNum
        subjs=length(sn);
        cwd = pwd;
        
        for s=1:subjs
            cd(fullfile(baseDir,imagingDirRaw,subj_name{sn(s)}));
            funScans = dir('*-r.nii');
            for i=1:length(funScans)  
                outfilename = sprintf('%s_run%2.2d.nii',subj_name{sn(s)},i);
                mkdir temp;
                spm_file_split(funScans(i).name,'temp');
                cd temp;
                list = dir([subj_name{sn(s)} '_run*.nii']);
                list = list(numDummys+1:end);  % Remove dummies
                V = {list(:).name};
                spm_file_merge(V,outfilename);
                movefile(outfilename,fullfile(baseDir,imagingDirRaw,subj_name{sn(s)}))
                cd ..
                rmdir('temp','s');
                fprintf('Run %d done for %s \n',i,subj_name{sn(s)});
            end
        end
        cd(cwd)
    case 'FUNC:fieldmap_RealignUnwarp'    % Realign functional images
        sn  = varargin{1}; % subjNum 
        subjs=length(sn);
        
        for s=1:subjs
            spmj_realign_unwarp(baseDir, subj_name{sn(s)}, run, 1, inf, 'prefix', '');
            fprintf('Done for %s \n',subj_name{sn(s)})
        end
    case 'FUNC:move_data'                 % Move realigned data
        % Moves image data from imaging_data_raw into imaging_data.
        % example: wmp_imana('FUNC:move_data',1,[1:10])
        sn=varargin{1}; % subjNum
        runs=varargin{2}; % runNum
        
        subjs=length(sn);
        
        for s=1:subjs
            dircheck(fullfile(baseDir,imagingDir,subj_name{sn(s)}))
            for r=1:length(runs)
                % move realigned data for each run
                source = fullfile(baseDir,imagingDirRaw,subj_name{sn(s)},sprintf('u%s_run%2.2d.nii',subj_name{sn(s)},runs(r)));
                dest = fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('%s_run%2.2d.nii',subj_name{sn(s)},runs(r)));
                copyfile(source,dest);
                
                % move realignment parameter files for each run
                source = fullfile(baseDir,imagingDirRaw,subj_name{sn(s)},sprintf('rp_%s_run%2.2d.txt',runs(r)));
                dest = fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('rp_%s_run%2.2d.txt',runs(r)));
                copyfile(source,dest);
            end
            
            % moving mean epi.
            source = fullfile(baseDir,imagingDirRaw,subj_name{sn(s)},sprintf('meanu%s_run01.nii',subj_name{sn(s)}));
            dest = fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('%umeanEPI_%s.nii',subj_name{sn(s)}));

            copyfile(source, dest);
            
            fprintf('realigned epi''s moved for %s \n',subj_name{sn(s)})
        end
    case 'FUNC:plot_motion_params'        % Quick inspection of head motion parameters
        cwd = pwd;
        sn=varargin{1}; % subjNum
        subjs=length(sn);
        for s=1:subjs
            X = cell(1,1);

            for r=1:length(run)
                X{r,1} = dlmread(fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('rp_%s_run%02d.txt',subj_name{sn(s)},r)));
            end
            X = cell2mat(X);
            clr = hsv(3); figure;
            subplot(2,1,1);
            for i = 1:3
                plot(X(:,i),'Color',clr(i,:));
                hold on;
            end
            legend('x','y','z','location','EastOutside');
            xlabel('Volume number'); ylabel('Deviation in mm'); ylim([-2 2]);
            title(subj_name{sn(s)}); set(gca,'fontsize',fs);
            hold on;
            
            % draw line marking each new run
            for r=1:length(run)
                plt.drawline(r*(numTRs-numDummys),'dir','vert','linestyle',':','linewidth',lw)
            end
            hold off;
            subplot(2,1,2);
            for j=1:3
                plot(X(:,j+3)*180/pi,'Color',clr(j,:));
                hold on;
            end
            legend('pitch','roll','yaw','location','EastOutside');
            xlabel('Volume number'); ylabel('Deviation in mm'); ylim([-2 2]);
            set(gca,'fontsize',fs); hold on;

            % draw line marking each new run
            for r=1:length(run)
                plt.drawline(r*(numTRs-numDummys),'dir','vert','linestyle',':','linewidth',lw)
            end
            hold off;
            cd(cwd);
            pause;
        end    
    case 'FUNC:meanimage_bias_correction' % Bias correct mean image prior to coregistration
        % example: wmp_imana('FUNC:meanimage_bias_correction',1)
        sn=varargin{1}; % subjNum
        subjs=length(sn);

        for s=1:subjs
            % make copy of original mean epi, and work on that
            source  = fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('%umeanEPI_%s.nii',subj_name{sn(s)}));
            dest    = fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('bumeanEPI_%s.nii',subj_name{sn(s)}));
            copyfile(source,dest);

            % bias correct mean image for grey/white signal intensities
            P{1}    = dest;
            spmj_bias_correct(P);
        end
    case 'FUNC:coreg'                     % Coregister meanEPI to anatomical image
        % example: wmp_imana('FUNC:coreg',1)
        % (1) Manually seed the functional/anatomical registration
        % - Do "coregtool" on the matlab command window
        % - Select anatomical image and meanepi image to overlay
        % - Manually adjust meanepi image and save result as rmeanepi
        sn=varargin{1}; % subjNum
        subjs=length(sn);
        
        cwd = pwd;
        
        for s=1:subjs
            cd(fullfile(baseDir,anatomicalDir,subj_name{sn(s)}));
            coregtool;
            keyboard();
            
            % (2) Automatically co-register functional and anatomical images
            J.ref    = {fullfile(baseDir,anatomicalDir,subj_name{sn(s)},'anatomical.nii')};
            J.source = {fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('rbbumeanEPI_%s.nii',subj_name{sn(s)}))};
            J.other = {''};
            J.eoptions.cost_fun = 'nmi';
            J.eoptions.sep = [4 2];
            J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            J.eoptions.fwhm = [7 7];
            matlabbatch{1}.spm.spatial.coreg.estimate = J;
            spm_jobman('run', matlabbatch);

            % (3) Manually check again
            coregtool;
            keyboard();
            cd(cwd);
        end
    case 'FUNC:make_samealign'            % Align to first image (rb*meanEPI_* of first session)
        sn=varargin{1}; % subjNum
        subjs=length(sn);
                        
        for s=1:subjs
            % Select image for reference
            P{1} = fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('rbbumeanEPI_%s.nii',subj_name{sn(s)}));
            % Select images to be realigned
            Q = {};
            for r=1:length(run)
                for i=1:numTRs-numDummys
                    Q{end+1} = fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('%s_run%02d.nii,%d',subj_name{sn(s)},r,i));
                end
            end
            % Run spmj_makesamealign_nifti to bring all functional runs into
            % same space as realigned mean epis
            spmj_makesamealign_nifti(char(P),char(Q));
            fprintf(1,'Done. Run spmj_checksamealign to check alignment.\n')
        end
    case 'FUNC:check_samealign'           % Check if all functional scans are align to the anatomical
        sn=varargin{1}; % subjNum
        subjs=length(sn);
                        
        for s=1:subjs
            % Select image for reference
            P{1} = fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('rbbumeanEPI_%s.nii',subj_name{sn(s)}));
            % Select images to be realigned
            Q = {};
            for r=1:length(run)
                for i=1:numTRs-numDummys
                    Q{end+1} = fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('%s_run%02d.nii,%d',subj_name{sn(s)},r,i));
                end
            end
            spmj_checksamealign(char(P),char(Q));
        end
    case 'FUNC:make_maskImage'            % Make mask images (noskull and gray_only)
        sn=varargin{1}; % subjNum
        subjs=length(sn);
        cwd = pwd;
        
        for s=1:subjs
            %cd(fullfile(baseDir,imagingDir,subj_name{sn(s)}));
            
            nam{1}  = fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('rbbumeanEPI_%s.nii',subj_name{sn(s)})); 
            nam{2}  = fullfile(baseDir,anatomicalDir,subj_name{sn(s)},'c1anatomical.nii');
            nam{3}  = fullfile(baseDir,anatomicalDir,subj_name{sn(s)},'c2anatomical.nii');
            nam{4}  = fullfile(baseDir,anatomicalDir,subj_name{sn(s)},'c3anatomical.nii');
            spm_imcalc(nam,'rmask_noskull.nii','i1>1 & (i2+i3+i4)>0.2')
            
            nam = {}; 
            nam{1}  = fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('rbbumeanEPI_%s.nii',subj_name{sn(s)})); 
            nam{2}  = fullfile(baseDir,anatomicalDir,subj_name{sn(s)},'c1anatomical.nii');
            spm_imcalc(nam,'rmask_gray.nii','i1>1 & i2>0.4');
            
            %cd(cwd)
        end

    case 'SURF:freesurfer'
        sn=varargin{1}; % subjNum
        subjs=length(sn);
        cwd = pwd;
        for s=1:subjs
            freesurfer_reconall(fullfile(baseDir,freesurferDir),subj_name{sn(s)},fullfile(baseDir,anatomicalDir,subj_name{sn(s)},'anatomical.nii'));
            cd(cwd)
        end
    case 'SURF:xhemireg' % cross-Register surfaces left / right hemi
        sn=varargin{1}; % subjNum
        subjs=length(sn);
        cwd = pwd;
        for s=1:subjs
            freesurfer_registerXhem(subj_name{sn(s)},fullfile(baseDir,freesurferDir),'hemisphere',[1 2]);
            cd(cwd)
        end
    case 'SURF:map_ico' % align to the new atlas surface (map icosahedron)
        sn=varargin{1}; % subjNum
        subjs=length(sn);
        cwd = pwd;
        for s=1:subjs
            freesurfer_mapicosahedron_xhem(subj_name{sn(s)},fullfile(baseDir,freesurferDir),'smoothing',1,'hemisphere',1:2);
            cd(cwd)
        end

    
end