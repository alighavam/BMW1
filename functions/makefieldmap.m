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
