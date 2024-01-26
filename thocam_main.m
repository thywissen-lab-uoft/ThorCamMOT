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

atom_type = 1; % 0:Rb, 1: K

%% Magnification and Pixel size

% CS165MU Pixelsize is 3.45 um
pixelsize0 = 3.45E-6; 

% Magnification depends on which camera you are using
% mag = 4.06;
mag = 3;

%% Analysis Flags

% Standard Analysis
doODProfile = 0;
doStandard  = 1;

doAnimate   = 1;
doProfile = 1;
doSave      = 1;

doMask = 1;

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
doGaussFit    = 1;      % Enable gauss fitting

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

    data.PixelSize = pixelsize0*mag;
    data.Atom = atom_type;
    atomdata(kk)=data;             
end
disp(' ');

atomdata = matchParamsFlags(atomdata);

%% Apply Mask

% if doMask
%     mask = load
%     for kk=1:length(atomdata)
%         atomdata(kk).Data.
%     end
% end

%% Compute OD if absorption image
[atomdata] = computeOD(atomdata);


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
% ROI=[550 800 500 750];

% ROI = [1 1440 1 1080];

ROI = [100 1100 100 1000];
% ROI=[400 900 400 900];

%% Aissgn the ROI

% Assign the ROI
disp(' ')
disp('Assigning ROI to data');
disp(ROI);

[atomdata.ROI]=deal(ROI);

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
end   

%% Gaussian

if doGaussFit   
    disp(repmat('-',1,60));    
    disp('Performing 2D gauss fit');
    disp(repmat('-',1,60));    
    % Iterate over all images (atomdata)
    for kk=1:length(atomdata)
        disp(repmat('-',1,60));   
        disp(['(' num2str(kk) ') ' atomdata(kk).Name]);
        % Iterate over all ROIs in an image
        for nn=1:size(atomdata(kk).ROI,1)   % Iterate over all ROIs
            sROI=atomdata(kk).ROI(nn,:);     % Grab the analysis ROI
            Dx=sROI(1):sROI(2);               % X Vector
            Dy=sROI(3):sROI(4);               % Y Vector
            data=atomdata(kk).Data(Dy,Dx);    % Optical density   
            
            [fout,gof,output]=gaussFit2D(Dx,Dy,data);    % Perform the fit               
            atomdata(kk).GaussFit{nn}=fout; % Assign the fit object       
            atomdata(kk).GaussGOF{nn}=gof; % Assign the fit object  
        end
    end
    gauss_data=getGaussData(atomdata,thor_xVar);  
    if doSave
        save([saveDir filesep 'gauss_data'],'gauss_data');
    end    
end 


%% standard

if doStandard
    thorcam_analysis_standard;
end


%% OD Profiles w or w/o Fits 
profile_opts = struct;
profile_opts.Style = 'cut'; 'sum';  % Cut or sum?
% profile_opts.Style = 'sum';  % Cut or sum?

profile_opts.FigLabel = FigLabel;

clear hF_X;clear hF_Y;
hF_X=[];hF_Y=[];

if doProfile

    for rNum=1:size(atomdata(1).ROI,1)
        profile_opts.ROINum = rNum;

        hF_Xs_rNum=showProfile(atomdata,'X',thor_xVar,profile_opts);

        if doSave
            for kk=1:length(hF_Xs_rNum) 
                figure(hF_Xs_rNum(kk));
                saveFigure(hF_Xs_rNum(kk),['OD_R' num2str(rNum) '_X' num2str(kk)],saveOpts);
                pause(0.1);
            end 
        end

        hF_Ys_rNum=showProfile(atomdata,'Y',thor_xVar,profile_opts);          
    %   Save the figures (this can be slow)
        if doSave        
            for kk=1:length(hF_Ys_rNum)
                figure(hF_Ys_rNum(kk));
                saveFigure(hF_Ys_rNum(kk),['OD_R' num2str(rNum) '_Y' num2str(kk)],saveOpts);
                pause(0.1);
            end
        end
        hF_X=[hF_X; hF_Xs_rNum];
        hF_Y=[hF_Y; hF_Ys_rNum];
    end  
 
end

%% Animate Cloud
if doAnimate
    animateOpts = struct;
    animateOpts.FigLabel = FigLabel;
    animateOpts.saveDir = saveDir;
    animateOpts.StartDelay = 1;
    animateOpts.MidDelay = 0.5;
    animateOpts.EndDelay = 1;
    animateOpts.doAverage = 1;
    animateOpts.doRotate = 0;
    animateOpts.xUnit = thor_unit;
    animateOpts.CLim = 'auto';
        animateOpts.CLim = [0 4];

    animateOpts.Order= 'ascend';

    animateCloud(atomdata,thor_xVar,animateOpts);
end
