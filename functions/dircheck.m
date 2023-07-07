function dircheck(dir)

% Checks existance of specified directory. Makes it if it does not exist.
% SArbuckle 01/2016
if ~exist(dir,'dir');
    %warning('%s didn''t exist, so this directory was created.\n',dir);
    mkdir(dir);
end
