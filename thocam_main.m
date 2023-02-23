% thorcam_main.m
% This is an imaging analysis script. It analyzes image taken from the PCO
% camera using our home made MATLAB code. It anaticipates loadining in .mat
% files which contain the PWA (probe with atoms) and PWOA (probe without
% atoms) images.
%
% All loaded .mat files are assumed to contain a structure with fields
%   Params
%   X
%   Y
%   Name
%   Date
%   Units

disp(repmat('-',1,60));disp([mfilename '.m']);disp(repmat('-',1,60)); 

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath))    

%% Close all non GUI figures
% Close all figures without the GUI tag.
figs=get(groot,'Children');
disp(' ');
disp('Closing all non GUI figures.');
for kk=1:length(figs)
   if ~isequal(figs(kk).Tag,'GUI')
       disp(['Closing figure ' num2str(figs(kk).Number) ' ' figs(kk).Name]);
      close(figs(kk)) 
   end
end
disp(' ');

%% Analysis Variable
% This section of code chooses the variable to plot against for aggregate
% plots.  The chosen variable MUST match a variable provided in the params
% field of the .mat file. The unit has no tangibile affect and only affects
% display properties.

% Defautl variable to plot against
thor_xVar = 'rf_freq_HF_shift';
% thor_xVar = 'ExecutionDate';

% Should the analysis attempt to automatically find the xvariable?
thor_autoXVar = 1;

% Should the analysis attempt to automatically find the unit?
thor_autoUnit = 1;

% If ixon_autoUnit=0, this will be used.
thor_overrideUnit='G'; 


%% Analysis Flags

% Standard Analysis
doODProfile = 0;
doStandard  = 1;

doAnimate   = 1;
doSave      = 1;

%%%%%%%%%%%%%%%%%%%%%%%%
% Standard Analyses
%%%%%%%%%%%%%%%%%%%%%%%%
% These analyses are "standard" in that they primarily fit the cloud to a
% typical distribution (gaussian, lattice, fermi). Additionally analyses
% on these processed data may be applied through the special flags. The
% processed data outputs of the below fits are typically <fit_type>_name

% Box Count
doBoxCount    = 1;      % Box count analysis

% Gaussian Fit
% Fit to a gaussian distribution (thermal cloud)
doGaussFit    = 0;      % Enable gauss fitting

%% GDrive Settings
GDrive_root = 'G:\My Drive\Lattice Shared\LabData';
doUpload = 0;       % Upload to google drive?

%% Select image directory
% Choose the directory where the images to analyze are stored
disp([datestr(now,13) ' Choose an image analysis folder...']);
dialog_title='Choose the root dire ctory of the images';


if getImageDir(datevec(now))
    newdir=uigetdir(getImageDir(datevec(now)),dialog_title);
    saveOpts = struct;

    if isequal(newdir,0)
        disp('Canceling.');    
        return; 
    else
        imgdir = newdir;
        saveDir = [imgdir filesep 'figures'];

        if ~exist(saveDir,'dir'); mkdir(saveDir);end    

        saveOpts.saveDir=saveDir;
        saveOpts.Quality = 'auto';

        strs=strsplit(imgdir,filesep);
        FigLabel=[strs{end-1} filesep strs{end}];
    end
else
    disp('Canceling.');
    return;
end

%% Load the data
clear atomdata
disp(['Loading data from ' imgdir]);
files=dir([imgdir filesep '*.mat']);
files={files.name};

for kk=1:length(files)
    str=fullfile(imgdir,files{kk});
    [a,b,c]=fileparts(str);      
    disp(['     (' num2str(kk) ')' files{kk}]);    
    data=load(str);     
    data=data.data;  

    % Display image properties
    try
        disp(['     Image Name     : ' data.Name]);
        disp(['     Execution Time : ' datestr(data.Date)]);
        if ~pco_autoXVar
            disp(['     ' thor_xVar ' : ' num2str(data.Params.(thor_xVar))]);
        end
        disp(' ');
    end    
    
    % Make sure executiondate is a number
    data.Params.ExecutionDate = datenum(data.Params.ExecutionDate);
    data.Params.ExecutionDateStr = datestr(data.Params.ExecutionDate);    
    data.Units.ExecutionDate =  'days';
    data.Units.ExecutionDateStr = 'str';
    
    atomdata(kk)=data;             
end
disp(' ');

atomdata = matchParamsFlags(atomdata);

%% X Variable and Units

if thor_autoXVar
    xVars = findXVars(atomdata);
    disp([' Found ' num2str(length(xVars)) ...
        ' valid variables that are changing to plot against.']);
    disp(xVars);
    
    % Select the first one
    ind = 1;    
    thor_xVar = xVars{ind};
    
    disp([' Setting ' thor_xVar ' to be the x-variable']);
    
    for kk=1:length(atomdata)
        disp([' (' num2str(kk) ') (' num2str(atomdata(kk).Params.(thor_xVar)) ') ' ...
            atomdata(kk).Name]); 
    end
    disp(' ');
end

% Grab the unit information
if thor_autoUnit && isfield(atomdata(1),'Units') 
    thor_unit=atomdata(1).Units.(thor_xVar);
else
    thor_unit=thor_overrideUnit;
end

% Sort the data by your given parameter
disp(['Sorting atomdata by the given ''' thor_xVar '''']);
x=zeros(length(atomdata),1);
for kk=1:length(atomdata)
    if isfield(atomdata(kk).Params,thor_xVar)
        x(kk)=atomdata(kk).Params.(thor_xVar) + 3;
    else
        warning(['atomdata(' num2str(kk) ') has no ''' thor_xVar '''']);
    end
end


% Sort it
[~, inds]=sort(x);
atomdata=atomdata(inds);
%% ROIS
  ROI=[450 800 250 500];

%% Aissgn the ROI

% Assign the ROI
disp(' ')
disp('Assigning ROI to data');
disp(ROI);

[atomdata.ROI]=deal(ROI);
%% Calculate Data
for kk=1:length(atomdata)
    atomdata(kk).OD = atomdata(kk).Images;
end
%% Box Count
% This section of code computes the box counts on all your data and ROIs.

boxOpts = struct;
boxOpts.doSubBG = 0;
boxOpts.bgROI = [700 790 500 600];

if doBoxCount
    disp(repmat('-',1,60));    
    disp('Performing box count analysis');
    disp(repmat('-',1,60));      
    atomdata=boxCount(atomdata,boxOpts);
    box_data = getBoxData(atomdata,thor_xVar);
    if doSave
        save([saveDir filesep 'box_data'],'box_data');
    end  
    
        
    if doSave && doUpload && exist(GDrive_root,'dir')
        gDir = [fileparts(getImageDir2(datevec(now),GDrive_root)) filesep FigLabel];
        gFile = [gDir filesep 'box_data'];        
        if ~exist(gDir,'dir')
           mkdir(gDir) 
        end
        save(gFile,'box_data');
    end
end   

%% standard

if doStandard
    thorcam_analysis_standard;
end

%% Plot the total counts
% hSize=showTotalCounts(atomdata,xVar,unit);
% % ylim([0 1E7]);
% saveFigure(atomdata, hSize, 'counts');

%% Plot the widths
% hSize=showSizes(atomdata,xVar);
% saveFigure(atomdata, hSize, 'widths');
%% Temperature Analysis
% if isequal(xVar,'tof_time') && length(atomdata)>2
%     [hTemp,fitX,fitY]=computeGaussianTemperature(atomdata);
% end

%% Plot Track Cloud Center

% if isequal(xVar,'tof_time')
%     trackCloudCenter(atomdata,xVar,'parabola')
% else
%     trackCloudCenter(atomdata,xVar);
% end

%% Animate Cloud
% animateCloud(atomdata,xVar);
