function thorcam_gui

%% Load dependencies
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath));


%% Open Existing GUI if available
guiname = 'ThorCamMOT';
a=groot;
for kk=1:length(a.Children)
    try
       if isequal(a.Children(kk).Name,guiname)
          warning('found an existing instance of gui')
          close(a.Children(kk)); 
%           return;
       end
    end
end


%% Load Libraries to run the camera

% Load TLCamera DotNet assembly.
% Directory where the dlls are located, this could depend on your comptuer
dll_dir = ['C:\Program Files\Thorlabs\Scientific Imaging\' ...
    'Scientific Camera Support\Scientific Camera Interfaces\MATLAB'];
cd(dll_dir);

dll_file=[pwd filesep 'Thorlabs.TSI.TLCamera.dll'];

% Check is the DLL's exist
if ~exist(dll_file,'file')
    error(['Could not located the DLL to run ' ...
        'the ThorCam! Exiting software.']);
end

% Load the NET assembly framework
fprintf('Loading the NET assembly ... ');
asmInfo=NET.addAssembly(dll_file);
asmInfo.Classes;
disp('Dot NET assembly loaded.');


% Open the SDK
fprintf('Opening the camera SDK...');
tlCameraSDK = Thorlabs.TSI.TLCamera.TLCameraSDK.OpenTLCameraSDK;
disp('loaded.');

%% Data Handles

tlCamera = [];

camera_settings = struct;
camera_settings.ExposureTime = 64;
camera_settings.Gain = 52;
camera_settings.PixelSize = 3.7; % size in um
camera_settings.TriggerMode = 1;

%% Graphics Options
h = 150;

%% Figure Handle
hF = figure('Name',guiname,'Color','w','units','pixels','toolbar','none',...
    'Tag','GUI','CloseRequestFcn',@closeGUI,...
    'NumberTitle','off','Position',[50 50 1200 800]);

% Callback for when the GUI is requested to be closed.
function closeGUI(fig,~)
    disp('Closing camera GUI...');
    try           
        tlCameraSDK.Dispose;                    % Unload the SDK
    catch ME
        warning(ME.message);
    end
    delete(fig);                % Delete the figure
end

function sizeChFcn(src,evt)
    hpC.Position(2) = hF.Position(4)-hpC.Position(4);
    hpAcq.Position(2) = hpC.Position(2);
end
hF.SizeChangedFcn=@sizeChFcn;
[255 204 0]/255;
%% Connection Panel
hpC=uipanel(hF,'units','pixels','backgroundcolor','w',...
    'Position',[1 hF.Position(4)-h 100 h],'title','connect');

% Refresh button
hbRefresh=uicontrol(hpC,'style','pushbutton','string','refresh','units','pixels',...
    'fontsize',8,'Position',[5 110 80 20],'backgroundcolor',[173 216 230]/255,...
    'Callback',@refreshCB,'enable','on');
hbCams = uicontrol(hpC,'Style','popupmenu','String',{'no cameras'},...
    'units','pixels','fontsize',8);
hbCams.Position = [5 80 80 20];


hbConnect=uicontrol(hpC,'style','pushbutton','string','connect','units','pixels',...
    'fontsize',8,'Position',[5 50 80 20],'backgroundcolor',[80 200 120]/255,...
    'Callback',@connectCB,'enable','off');
hbDisconnect=uicontrol(hpC,'style','pushbutton','string','disconnect','units','pixels',...
    'fontsize',8,'Position',[5 25 80 20],'backgroundcolor',[255 102 120]/255,...
    'Callback',@disconnectCB,'enable','off');

    function connectCB(src,evt)
        try
            sn = hbCams.String{hbCams.Value};        
            if isequal(sn,'no cameras')
                warning('cannot connect to no camera');
                return;
            end
            tlCamera=openCamera(sn,camera_settings);
            hbDisconnect.Enable='on';
            hbConnect.Enable='off';
        catch ME
            warning(ME.message);
        end
    end

    function disconnectCB(src,evt)
        try
            closeCamera(tlCamera);
            hbDisconnect.Enable='off';
            hbConnect.Enable='on';
        catch ME
            warning(ME.message)
        end
    end


    function refreshCB(src,evt)
        try
            sns = getCameras;
            hbCams.String = sns;
            hbConnect.Enable = 'on';
        catch ME
            warning(ME.message);
            sns = {'no cameras'};
            hbConnect.Enable = 'off';
        end
    end
%% Acquisition Panel
hpAcq=uipanel(hF,'units','pixels','backgroundcolor','w',...
    'Position',[100 hF.Position(4)-h 150 h],'title','acquisition');

ttstr='Start the camera and image acquisition';
hbstart=uicontrol(hpAcq,'style','pushbutton','string','start','units','pixels',...
    'fontsize',10,'Position',[5 110 40 20],'backgroundcolor',[80 200 120]/255,...
    'Callback',@startCamCB,'ToolTipString',ttstr,'enable','off');

% Clear the camera buffer
ttstr='Clear the camera buffer.';
hbclear=uicontrol(hpAcq,'style','pushbutton','string','clear',...
    'units','pixels','fontsize',10,'Position',[50 110 40 20],'enable','off',...
    'backgroundcolor',[255 204 0]/255,'callback',@clearBuffer,...
    'ToolTipString',ttstr);

% Stop acquisition button
ttstr='Stop the camera.';
hbstop=uicontrol(hpAcq,'style','pushbutton','string','stop',...
    'units','pixels','fontsize',10,'Position',[95 110 40 20],'enable','off',...
    'backgroundcolor',[255 102 120]/255,'callback',@stopCamCB,...
    'ToolTipString',ttstr);

% Button group for deciding what the X/Y plots show
bgAcq = uibuttongroup(hpAcq,'units','pixels','backgroundcolor','w',...
    'BorderType','None');  
bgAcq.Position(3:4)=[125 80];
bgAcq.Position(1:2)=[5 5];
    
% Radio buttons for cuts vs sum
uicontrol(bgAcq,'Style','radiobutton','String','live',...
    'Position',[0 20 100 20],'units','pixels','backgroundcolor','w','Value',1);
uicontrol(bgAcq,'Style','radiobutton','String','trigered',...
    'Position',[0 0 100 20],'units','pixels','backgroundcolor','w');


    function startCamCB(~,~)

    end

    function stopCamCB(~,~)

    end

    function clearBuffer(~,~)

    end

tbl_acq=uitable('parent',hpAcq,'units','pixels','RowName',{},'ColumnName',{},...
    'fontsize',8,'ColumnWidth',{100,40},'columneditable',[false true]);
tbl_acq.Data={...
    ['raw pixelsize (' char(956) 'm)'], raw_pixel_size*1E6;
    'magnification',mag(1)};



%% Camera Functions
  
%Get cameras
function sns = getCameras
    try
        fprintf('Refreshing cameras ...')
        sns={};
        SNs = tlCameraSDK.DiscoverAvailableCameras;
        for ii = 1:SNs.Count           
            camName=char(SNs.Item(ii-1));   
            sns{ii} = camName;
            fprintf([sns{ii} '...']);
        end           
        disp('found');
    catch ME
        disp('uh oh')
        warning(ME.message)
        sns = {};
    end
end

% Open camera with a given serial number, default to software trigger
    function tlCamera=openCamera(SN,settings)     
    fprintf(['opening ' SN '...']);
    try
        tlCamera = tlCameraSDK.OpenCamera(SN, false);        
        % Get and Set camera parameters
        tlCamera.ExposureTime_us = uint32(settings.ExposureTime);        
        % The default black level should be zero
        tlCamera.BlackLevel=uint32(0);       
        % Set the default gain level to zero
        tlCamera.Gain=tlCamera.ConvertDecibelsToGain(uint32(settings.Gain));                
        % Set operation mode to software testing
        tlCamera.OperationMode = ...
            Thorlabs.TSI.TLCameraInterfaces.OperationMode.SoftwareTriggered;    
        % Set frames per trigger to one
        tlCamera.FramesPerTrigger_zeroForUnlimited = 1;            
        % Set default ROI to max
        tlCamera.ROIAndBin.ROIHeight_pixels=tlCamera.SensorHeight_pixels;
        tlCamera.ROIAndBin.ROIWidth_pixels=tlCamera.SensorWidth_pixels;
        % Set the trigger polarity to high
        tlCamera.TriggerPolarity=...
            Thorlabs.TSI.TLCameraInterfaces.TriggerPolarity.ActiveHigh; 
        % Turn LED Off (it can add background signal)
        tlCamera.IsLEDOn=0;      
        % Arm the camera
        tlCamera.Arm;       
        disp('connected')
    catch ME
        disp('oh no')
        warning(ME.message)
        tlCamera = [];
    end
end

% Close the camera
function closeCamera(tlCamera)
    fprintf('Closing camera...');
    try
        if (tlCamera.IsArmed)
            fprintf('disarming...')
            tlCamera.Disarm;
        end
        fprintf('disposing...')
        tlCamera.Dispose;                
        delete(tlCamera);           
        disp('done')
    catch ME
        disp('oh no')
        warning(ME.message);
    end
end

% Read the hardware ROI on the camera
function [ROI,X,Y]=readROI
   % Read the X ROI
    x1=cam.ROIAndBin.ROIOriginX_pixels;
    W=cam.ROIAndBin.ROIWidth_pixels;
    
    % Read the Y ROI
    y1=cam.ROIAndBin.ROIOriginY_pixels;
    H=cam.ROIAndBin.ROIHeight_pixels;    
    
    % Read the total sensor size
    H0=cam.SensorHeight_pixels;
    W0=cam.SensorWidth_pixels;    
    
    % Redefine the ROI for a camera that is flipped
    if  isequal(sn,'10148')
        x1=W0-W+x1;
        y1=H0-H+y1;
    end   
    
    % Final processing on ROI
    ROI=double([x1+1 W+x1 y1+1 y1+H]);     
    
    % Output a pixel position vector for simplicity
    X=ROI(1):ROI(2);
    Y=ROI(3):ROI(4);
end

% Set the hardware ROI
function setROI(newROI)
    % Read the current Region of Interest (ROI)
    H0=cam.SensorHeight_pixels;
    W0=cam.SensorWidth_pixels;
    
    % Calculate the new widths
    W=newROI(2)-newROI(1)+1;
    H=newROI(4)-newROI(3)+1;
    
    % Get the new origins
    x1=newROI(1)-1;
    y1=newROI(3)-1;        
    
    % Change the speciifications if using rotated camera
    if isequal(sn,'10148')
        x1=W0-newROI(2);
        y1=H0-newROI(4);
    end       
    
   % Arm the camera
    cam.Disarm;        
            
    % Set the ROI height and width
    cam.ROIAndBin.ROIHeight_pixels=uint32(H);
    cam.ROIAndBin.ROIWidth_pixels=uint32(W);
    
    % Set the ROI origin
    cam.ROIAndBin.ROIOriginX_pixels=x1;
    cam.ROIAndBin.ROIOriginY_pixels=y1;
    
    % Read the ROI to verify
    [ROI,xVec,yVec]=readROI;      
    
    % Arm the camera
    cam.Arm;
end

end

