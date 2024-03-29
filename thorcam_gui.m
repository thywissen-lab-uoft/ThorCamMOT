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

historyDir=['C:' filesep 'ImageHistory'];

%% Load Libraries to run the camera

% Load TLCamera DotNet assembly.
% Directory where the dlls are located, this could depend on your comptuer
dll_dir = ['C:\Program Files\Thorlabs\Scientific Imaging\' ...
    'Scientific Camera Support\Scientific Camera Interfaces\MATLAB'];
dll_file=[dll_dir filesep 'Thorlabs.TSI.TLCamera.dll'];

if exist(dll_file,'file') 
    cd(dll_dir);
    addpath(dll_dir);addpath(genpath(dll_dir));

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
else
    warning('no SDK detected, energy debug mode');
    tlCameraSDK = [];
end

%% Data Handles

tlCamera = [];

camera_settings = struct;
camera_settings.ExposureTime = 708;
camera_settings.Gain_dB = 0;
camera_settings.Gain = 0;

camera_settings.PixelSize = 3.7; % size in um
camera_settings.TriggerMode = 1;
camera_settings.QuantumEfficiency = 0.25; % size in um
camera_settings.Magnification = 1;
camera_settings.SolidAngle = .2;
camera_settings.Images = [];

X=1:1392;                       % X pixel vector
Y=1:1024;                       % Y pixel vector
Z=zeros(length(Y),length(X));   % Image to show

%% Timer Objects

mytimer=timer('Name','liveupdate','executionmode','fixedspacing',...
    'period',0.001,'TimerFcn',@liveCB);
% Callback function for live update
    function liveCB(~,~)
        tlCamera.IssueSoftwareTrigger;        
        pause(0.02);
        updateImage;
    end

    function trigCB(~,~,nImages)
        if tlCamera.NumberOfQueuedFrames           
            img=grabImage;
            if isempty(camera_settings.Images)
                camera_settings.Images = img;
            else
                camera_settings.Images(:,:,end+1)=img;  
            end
        end

        if size(camera_settings.Images,3)==nImages
            newData(camera_settings.Images)
            camera_settings.Images = [];
        end   
    end


    function newData(imgs)  

        data=processImages(imgs);
        hImg.CData=data.Data;

        disp('     New Image!');
        disp(['     Image     : ' data.Name]); 


        % Update parameters table                            
        [~,inds] = sort(lower(fieldnames(data.Params)));
        params = orderfields(data.Params,inds);  
                
        fnames=fieldnames(params);
        for nn=1:length(fnames)
            tbl_params.Data{nn,1}=fnames{nn};
            val=data.Params.(fnames{nn});
            if isa(val,'double')
                tbl_params.Data{nn,2}=num2str(val);
            end

            if isa(val,'struct')
               tbl_params.Data{nn,2}='[struct]'; 
            end  
        end
                
        % Update flags table
        fnames=fieldnames(data.Flags);
        for nn=1:length(fnames)
            tbl_flags.Data{nn,1}=fnames{nn};
            val=data.Flags.(fnames{nn});
            if isa(val,'double')
                tbl_flags.Data{nn,2}=num2str(val);
            end                    
            if isa(val,'struct')
               tbl_flags.Data{nn,2}='[struct]'; 
            end                    
        end        

        saveData(data);

        % Save image to folder
        if hcSave.Value
           saveData(data,tSaveDir.UserData); 
        end                          

  end

function saveData(data,saveDir)
        if nargin==1
           saveDir=historyDir;
           filenames=dir([saveDir filesep '*.mat']);
           filenames={filenames.name};
           filenames=sort(filenames);

           % Delete old images
           if length(filenames)>200
               f=[saveDir filesep filenames{1}];
               delete(f);
           end               
        end
        fname=[data.Name '.mat']; 
        if ~exist(saveDir,'dir')
           mkdir(saveDir);
        end        
        fname=fullfile(saveDir,fname);
        fprintf('%s',[fname ' ...']);
        save(fname,'data');
        disp(' done'); 
    end

%% Graphics Options
h = 110;
h2=30;
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
    hpOptics.Position(2) = hpC.Position(2);
    hpImgProcess.Position(2) = hpC.Position(2);
    hpRaw.Position(2) = hpC.Position(2);
    hpAnl.Position(2) = hpC.Position(2);
    hpFit.Position(4) = hF.Position(4) - hpC.Position(4);
    hpSave.Position(2) = hF.Position(4) - hpC.Position(4) - hpSave.Position(4);
    hpSave.Position(3) = hF.Position(3)-hpFit.Position(1);
    hp.Position(3:4) = [hF.Position(3)-hp.Position(1) hF.Position(4)-h-h2];
%     resizePlots;
end
hF.SizeChangedFcn=@sizeChFcn;

%% Connection Panel
hpC=uipanel(hF,'units','pixels','backgroundcolor','w',...
    'Position',[1 hF.Position(4)-h 88 h],'title','connect');

% Refresh button
hbRefresh=uicontrol(hpC,'style','pushbutton','string','refresh','units','pixels',...
    'fontsize',8,'Position',[2 74 80 20],'backgroundcolor',[173 216 230]/255,...
    'Callback',@refreshCB,'enable','on');
hbCams = uicontrol(hpC,'Style','popupmenu','String',{'no cameras'},...
    'units','pixels','fontsize',8);
hbCams.Position = [2 50 80 20];

% Connect button
hbConnect=uicontrol(hpC,'style','pushbutton','string','connect','units','pixels',...
    'fontsize',8,'Position',[2 26 80 20],'backgroundcolor',[80 200 120]/255,...
    'Callback',@connectCB,'enable','off');

% Disonnect button
hbDisconnect=uicontrol(hpC,'style','pushbutton','string','disconnect','units','pixels',...
    'fontsize',8,'Position',[2 4 80 20],'backgroundcolor',[255 102 120]/255,...
    'Callback',@disconnectCB,'enable','off');

% Connect to current camera
    function connectCB(src,evt)
        try
            sn = hbCams.String{hbCams.Value};        
            if isequal(sn,'no cameras')
                warning('cannot connect to no camera');
                return;
            end
            tlCamera=openCamera(sn,camera_settings);
            
            tblUpdate;
            hbDisconnect.Enable='on';
            hbConnect.Enable='off';
            hbstart.Enable='on';
            hbstop.Enable='off';
            hbclear.Enable='on';
            tbl_acq.Enable = 'on';

            for t=1:length(bgAcq.Children)
                bgAcq.Children(t).Enable='on';
            end
            bgAcq.Children(end).Value = 1;
            mytimer.TimerFcn = @liveCB;

        catch ME
            warning(ME.message);
        end
    end

% Disconnect from current camera
    function disconnectCB(src,evt)
        try
            closeCamera(tlCamera);
            hbDisconnect.Enable='off';
            hbConnect.Enable='on';

            hbstart.Enable='off';
            hbstop.Enable='off';
            hbclear.Enable='off';
            tbl_acq.Enable = 'off';
            for t=1:length(bgAcq.Children)
                bgAcq.Children(t).Enable='off';
            end
            bgAcq.Children(end).Value = 1;
            mytimer.TimerFcn = @liveCB;

        catch ME
            warning(ME.message)
        end
    end

% Refresh camera list
    function refreshCB(src,evt)
        try
            sns = getCameras;
            hbCams.String = sns;
            hbConnect.Enable = 'on';
        catch ME
            warning(ME.message);
            sns = {'no cameras'};
            hbCams.String = sns;
            hbConnect.Enable = 'off';
        end
    end
%% Acquisition Panel
hpAcq=uipanel(hF,'units','pixels','backgroundcolor','w',...
    'Position',[hpC.Position(1)+hpC.Position(3) hF.Position(4)-h 330 h],'title','acquisition');

% Start acquisitino button
ttstr='Start the camera and image acquisition';
hbstart=uicontrol(hpAcq,'style','pushbutton','string','start','units','pixels',...
    'fontsize',10,'Position',[2 74 40 20],'backgroundcolor',[80 200 120]/255,...
    'Callback',@startCamCB,'ToolTipString',ttstr,'enable','off');

% Clear the camera buffer
ttstr='Clear the camera buffer.';
hbclear=uicontrol(hpAcq,'style','pushbutton','string','clear',...
    'units','pixels','fontsize',10,'Position',[2 52 40 20],'enable','off',...
    'backgroundcolor',[255 204 0]/255,'callback',@clearBuffer,...
    'ToolTipString',ttstr);

% Stop acquisition button
ttstr='Stop the camera.';
hbstop=uicontrol(hpAcq,'style','pushbutton','string','stop',...
    'units','pixels','fontsize',10,'Position',[2 30 40 20],'enable','off',...
    'backgroundcolor',[255 102 120]/255,'callback',@stopCamCB,...
    'ToolTipString',ttstr);

% Button group for deciding what the X/Y plots show
bgAcq = uibuttongroup(hpAcq,'units','pixels','backgroundcolor','w',...
    'BorderType','None','SelectionChangedFcn',@modeChangeCB);  
bgAcq.Position(3:4)=[140 80];
bgAcq.Position(1:2)=[50 h-(80+15)];

    function modeChangeCB(src,evt)
        if evt.NewValue.UserData == 0 % Go into live mode
            try
                tlCamera.Disarm;
                tlCamera.OperationMode=Thorlabs.TSI.TLCameraInterfaces.OperationMode.SoftwareTriggered;
                mytimer.TimerFcn = @liveCB;
            catch ME
                warning(ME.message);
            end
        else % Triggered mode
            try
                nImages = evt.NewValue.UserData;
                tlCamera.Disarm;
                tlCamera.OperationMode=Thorlabs.TSI.TLCameraInterfaces.OperationMode.HardwareTriggered;
                mytimer.TimerFcn = @(src,evt) trigCB(src,evt,nImages);
            catch ME
                warning(ME.message);
            end
        end
    end

       
        
           
    
    
% Radio buttons for cuts vs sum
uicontrol(bgAcq,'Style','radiobutton','String','live','fontsize',7,...
    'Position',[0 60 120 20],'units','pixels','backgroundcolor','w','Value',1,'enable','off',...
    'UserData',0);
uicontrol(bgAcq,'Style','radiobutton','String','trig (PWA)','fontsize',7,...
    'Position',[0 40 120 20],'units','pixels','backgroundcolor','w','enable','off',...
    'UserData',1);
uicontrol(bgAcq,'Style','radiobutton','String','trig (PWA, PWOA)','fontsize',7,...
    'Position',[0 20 120 20],'units','pixels','backgroundcolor','w','enable','off',...
    'UserData',2);
uicontrol(bgAcq,'Style','radiobutton','String','trig (PWA, PWOA, dark)','fontsize',7,...
    'Position',[0 0 120 20],'units','pixels','backgroundcolor','w','enable','off',...
    'UserData',3);
bgAcq.Children
% Start camera callback
    function startCamCB(~,~)
        tlCamera.Arm;
        start(mytimer);
        hbstart.Enable='off';
        hbstop.Enable='on';
        for t=1:length(bgAcq.Children)
            bgAcq.Children(t).Enable='off';
        end

    end

% Stop camera callback
    function stopCamCB(~,~)
        tlCamera.Disarm;

        stop(mytimer);
        hbstart.Enable='on';
        hbstop.Enable='off';
        for t=1:length(bgAcq.Children)
            bgAcq.Children(t).Enable='on';
        end
    end

% Clear buffer calblack
    function clearBuffer(~,~)

    end

tbl_acq=uitable('parent',hpAcq,'units','pixels','RowName',{},'ColumnName',{},...
    'fontsize',8,'ColumnWidth',{100,40},'columneditable',[false true],...
    'celleditcallback',@chSet,'Enable','off');
tbl_acq.Data={...
    'gain (dB)', camera_settings.Gain_dB;
    'gain', camera_settings.Gain;
    'exposure time (us)',camera_settings.ExposureTime};
tbl_acq.Position(3:4) = tbl_acq.Extent(3:4);
tbl_acq.Position(1:2)=[170 10];


 function chSet(tbl,data)
        r=data.Indices(1);
        c=data.Indices(2);
        val=data.NewData;            
        % Gain goes for 0 to 48 dB
        % Exposure goes from 64us to 51925252us

        switch r
            case 1                    
                if val>=0 && val<=48
                    val=round(val,1);
                    disp(['Changing gain to ' num2str(val) ' dB']);
                    gVal=tlCamera.ConvertDecibelsToGain(val);
                    tlCamera.Gain=gVal;                        
                    gGain=tlCamera.ConvertGainToDecibels(tlCamera.Gain);
                    tbl.Data{r,c}=gGain;                        
                else
                    warning('Valid gain is 0 dB to 48 dB')
                    tbl.Data{r,c}=data.PreviousData;                        
                end            
            case 3
                if val>=64 && val<=1E5
                    val=round(val);
                    disp(['Changing exposure to ' num2str(val) ' us']);
                    tlCamera.ExposureTime_us=uint32(val);
                    tbl.Data{r,c}=double(tlCamera.ExposureTime_us);
                else
                    tbl.Data{r,c}=data.PreviousData;    
                    warning('Valid expsure is 64u to 1e5us');
                end
        end   

        tblUpdate;
 end

    function tblUpdate
        camera_settings.Gain = tlCamera.Gain;
        camera_settings.Gain_dB=tlCamera.ConvertGainToDecibels(...
            tlCamera.Gain);
        camera_settings.ExposureTime=tlCamera.ExposureTime_us;
        tbl_acq.Data{1,2} = camera_settings.Gain_dB;
        tbl_acq.Data{2,2} = camera_settings.Gain;
        tbl_acq.Data{3,2} = camera_settings.ExposureTime;
        updateDescStr;
    end
%% Optics

hpOptics=uipanel('parent',hF,'units','pixels','backgroundcolor','w',...
    'title','optics');
hpOptics.Position=[hpAcq.Position(1)+hpAcq.Position(3) hF.Position(4)-h 160 h]; 

tbl_optics=uitable('parent',hpOptics,'units','pixels','RowName',{},'ColumnName',{},...
    'fontsize',8,'ColumnWidth',{90,40},'columneditable',[false true]);
tbl_optics.Data={...
    'Magnification', camera_settings.Magnification;
    'Pixel Size (um)',camera_settings.PixelSize;
    'QE',camera_settings.QuantumEfficiency;
    'Solid Angle (ster)',camera_settings.SolidAngle};
tbl_optics.Position(3:4) = tbl_optics.Extent(3:4);
tbl_optics.Position(1:2)=[5 10];
%% Image Process

hpImgProcess=uipanel('parent',hF,'units','pixels','backgroundcolor','w',...
    'title','processing');
hpImgProcess.Position=[hpOptics.Position(1)+hpOptics.Position(3) hF.Position(4)-h 200 h]; 


%% Image Process

hpRaw=uipanel('parent',hF,'units','pixels','backgroundcolor','w',...
    'title','raw images');
hpRaw.Position=[hpImgProcess.Position(1)+hpImgProcess.Position(3) hF.Position(4)-h 400 h]; 

%% Image Process

hpAnl=uipanel('parent',hF,'units','pixels','backgroundcolor','w',...
    'title','analysis');
hpAnl.Position=[hpRaw.Position(1)+hpRaw.Position(3) hF.Position(4)-h 200 h]; 

%% Fit Results Panel
hpFit=uitabgroup(hF,'units','pixels');
hpFit.Position=[0 0 220 hF.Position(4)-h];

tabs(1)=uitab(hpFit,'Title','params','units','pixels');
tabs(2)=uitab(hpFit,'Title','flags','units','pixels');
tabs(3)=uitab(hpFit,'Title','1','units','pixels');

% Table for run parameters
tbl_params=uitable(tabs(1),'units','normalized','RowName',{},'fontsize',8,...
    'ColumnName',{},'ColumnWidth',{135 60},'columneditable',[false false],...
    'Position',[0 0 1 1]);
% Table for run parameters
tbl_flags=uitable(tabs(2),'units','normalized','RowName',{},'fontsize',7,...
    'ColumnName',{},'ColumnWidth',{145 50},'columneditable',[false false],...
    'Position',[0 0 1 1]);

% Table for analysis outputs
tbl_analysis(1)=uitable(tabs(3),'units','normalized','RowName',{},'ColumnName',{},...
    'fontsize',8,'ColumnWidth',{60 65 65},'columneditable',false(ones(1,3)),...
    'Position',[0 0 1 1]);

%% Save Data

hpSave=uipanel('parent',hF,'units','pixels','backgroundcolor','w');
hpSave.Position(1:2)=[hpFit.Position(1)+hpFit.Position(3) hF.Position(4)-h-h2];
hpSave.Position(3:4) = [hF.Position(3)-hpFit.Position(1) h2];


% Auto Save check box
ttstr=['Enable/Disable saving to external directory. Does ' ...
    'not override saving to image history.'];
hcSave=uicontrol(hpSave,'style','checkbox','string','save?','fontsize',8,...
    'backgroundcolor','w','Position',[5 5 60 20],'callback',@saveCheck,...
    'ToolTipString',ttstr);

% Save checkbox callback
    function saveCheck(src,~)
        if src.Value
            tSaveDir.Enable='on';
            bBrowse.Enable='on';
        else
            tSaveDir.Enable='off';
            bBrowse.Enable='off';
        end
    end

% Browse button
cdata=imresize(imread('images/browse.jpg'),[20 20]);
bBrowse=uicontrol(hpSave,'style','pushbutton','CData',cdata,'callback',@browseCB,...
    'enable','off','backgroundcolor','w','position',[60 5 size(cdata,[1 2])]);

% String for current save directory
tSaveDir=uicontrol(hpSave,'style','text','string','directory','fontsize',8,...
    'backgroundcolor','w','units','pixels','horizontalalignment','left',...
    'enable','off','UserData','','Position',[85 0 hF.Position(3)-290 22]);

% Browse button callback
    function browseCB(~,~)
        str=getDayDir;
        str=uigetdir(str);
        if str
            tSaveDir.UserData=str; % Full directory to save
            str=strsplit(str,filesep);
            str=[str{end-1} filesep str{end}];
            tSaveDir.String=str; % display string
        else
            disp('no directory chosen!');
        end
    end

%% Main Image

hp=uipanel('parent',hF,'units','pixels','backgroundcolor','w');
hp.Position(1:2)=[hpFit.Position(1)+hpFit.Position(3) 0];
hp.Position(3:4) = [hF.Position(3)-hp.Position(1) hF.Position(4)-h-h2];


axImg=axes('parent',hp,'UserData','OD');cla
hImg=imagesc(X,Y,Z);
set(axImg,'box','on','linewidth',.1,'fontsize',10,'units','normalized',...
    'XAxisLocation','top');
hold on
% axImg.Position=[50 150 hp.Position(3)-200 hp.Position(4)-200];
axis equal tight
colormap(inferno);
% Box for ROI (this will become an array later)
% pROI=rectangle('position',[1 1 1392 1024],'edgecolor',co(1,:),'linewidth',2);

cBar=colorbar('fontsize',8,'units','pixels','location','eastoutside');
drawnow;


tImgDesc=text(4,4,'test','units','pixels','verticalalignment','bottom',...
    'color','r','fontweight','bold','fontsize',12);

    function updateDescStr(counts)
       str = [num2str(camera_settings.Gain_dB) ' dB ' ...
           num2str(camera_settings.ExposureTime)  ' \mus'];
       if nargin == 1
           str = [str newline sprintf('%.4e',counts)];
       end
       tImgDesc.String =str;
    end


updateDescStr

%% Helper Functions
function s3=getDayDir
    t=now;
    d=['Y:\Data'];
    s1=datestr(t,'yyyy');s2=datestr(t,'yyyy.mm');s3=datestr(t,'mm.dd');
    s1=[d filesep s1];s2=[s1 filesep s2];s3=[s2 filesep s3];

    if ~exist(s1,'dir'); mkdir(s1); end
    if ~exist(s2,'dir'); mkdir(s2); end
    if ~exist(s3,'dir'); mkdir(s3); end
end

%%


ttstr='Maximize display ROI to full image size.';
cdata=imresize(imread('images/fullLim.png'),[15 15]);
hbFullLim=uicontrol(hp,'style','pushbutton','Cdata',cdata,'Fontsize',10,...
    'Backgroundcolor','w','Position',[1 1 21 20],'Callback',@fullDispCB,...
    'ToolTipString',ttstr);

ttstr='Snap display ROI to data ROI(s).';
cdata=imresize(imread('images/snapLim.png'),[15 15]);
hbSnapLim=uicontrol(hp,'style','pushbutton','Cdata',cdata,'Fontsize',10,...
    'Backgroundcolor','w','Position',[22 1 21 20],'Callback',@snapDispCB,...
    'ToolTipString',ttstr);

% Button to enable GUI selection of display limits
ttstr='Select the display ROI.';
cdata=imresize(imread('images/target.jpg'),[15 15]);
hbSlctLim=uicontrol(hp,'style','pushbutton','Cdata',cdata,'Fontsize',10,...
    'Backgroundcolor','w','Position',[44 1 20 20],'Callback',@slctDispCB,...
    'ToolTipString',ttstr);

% Table for changing display limits
tbl_dispROI=uitable('parent',hp,'units','pixels','RowName',{},'columnname',{},...
    'ColumnEditable',[true true true true],'CellEditCallback',@tbl_dispROICB,...
    'ColumnWidth',{30 30 30 30},'FontSize',8,'Data',[1 size(Z,2) 1 size(Z,1)]);
tbl_dispROI.Position(3:4)=tbl_dispROI.Extent(3:4);
tbl_dispROI.Position(1:2)=[66 1];

    function tbl_dispROICB(src,evt)
        ROI=src.Data;        % Grab the new ROI     
        % Check that the data is numeric
        if sum(~isnumeric(ROI)) || sum(isinf(ROI)) || sum(isnan(ROI))
            warning('Incorrect data type provided for ROI.');
            src.Data(evt.Indices(2))=evt.PreviousData;
            return;
        end        
        ROI=round(ROI);      % Make sure this ROI are integers   

        % Keep the ROI within image bounds (this is hardcoded and could be
        % changed if we ever implement hardware ROI but want to keep 
        % absolute pixel positions relative to total sensor.)
        if ROI(2)<=ROI(1) || ROI(4)<=ROI(3)
           warning('Bad ROI specification given.');
           ROI(evt.Indices(2))=evt.PreviousData;
        end       
        if ROI(1)<1; ROI(1)=1; end       
        if ROI(3)<1; ROI(3)=1; end   
        if ROI(4)>size(dstruct.PWA,1); ROI(4)=size(dstruct.PWA,1);end       
        if ROI(2)>size(dstruct.PWA,2); ROI(2)=size(dstruct.PWA,2);end       
        src.Data=ROI;       
        try
            set(axImg,'XLim',ROI(1:2),'YLim',ROI(3:4));
            set(axPWA,'XLim',axImg.XLim,'YLim',axImg.YLim);
            set(axPWOA,'XLim',axImg.XLim,'YLim',axImg.YLim);
            set(axDark,'XLim',axImg.XLim,'YLim',axImg.YLim);

            
            resizePlots;
            drawnow;
%             pDisp.Position=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];           
            updateScalebar;
            drawnow;
        catch ab
            warning('Unable to change display ROI.');
            src.Data(evt.Indices)=evt.PreviousData;
        end
    end

    function fullDispCB(~,~)
       ROI=[1 size(dstruct.PWA,2) 1 size(dstruct.PWOA,1)];
       tbl_dispROI.Data=ROI;
       tbl_dispROICB(tbl_dispROI);
        set(axPWA,'XLim',axImg.XLim,'YLim',axImg.YLim);
        set(axPWOA,'XLim',axImg.XLim,'YLim',axImg.YLim);
        set(axDark,'XLim',axImg.XLim,'YLim',axImg.YLim);

       resizePlots;
       drawnow;
    end

    function snapDispCB(~,~)
       ROI=[min(tblROI.Data(:,1)) max(tblROI.Data(:,2)) ...
           min(tblROI.Data(:,3)) max(tblROI.Data(:,4))];
       tbl_dispROI.Data=ROI;
       tbl_dispROICB(tbl_dispROI);
        set(axPWA,'XLim',axImg.XLim,'YLim',axImg.YLim);
        set(axPWOA,'XLim',axImg.XLim,'YLim',axImg.YLim);
        set(axDark,'XLim',axImg.XLim,'YLim',axImg.YLim);

       resizePlots;
       drawnow;
    end

    function slctDispCB(~,~)
        disp(['Selecting display ROI .' ...
            ' Click two points that form the rectangle ROI.']);
        axes(axImg)                 % Select the OD image axis
        [x1,y1]=ginput(1);          % Get a mouse click
        x1=round(x1);y1=round(y1);  % Round to interger        
        p1=plot(x1,y1,'+','color','k','linewidth',1); % Plot it
        
        [x2,y2]=ginput(1);          % Get a mouse click
        x2=round(x2);y2=round(y2);  % Round it        
        p2=plot(x2,y2,'+','color','k','linewidth',1);  % Plot it

        % Create the ROI
        ROI=[min([x1 x2]) max([x1 x2]) min([y1 y2]) max([y1 y2])];

        % Constrain ROI to image
        if ROI(1)<1; ROI(1)=1; end       
        if ROI(3)<1; ROI(3)=1; end   
        if ROI(4)>size(dstruct.PWA,1); ROI(4)=size(dstruct.PWA,2); end       
        if ROI(2)>size(dstruct.PWA,2); ROI(2)=size(dstruct.PWA,2); end   
        
        % Try to update ROI graphics
        tbl_dispROI.Data=ROI;
        tbl_dispROICB(tbl_dispROI);
        resizePlots;       
        drawnow;        
        delete(p1);delete(p2);                   % Delete markers
    end
%%
% Text label for color limit table on OD image
climtext=uicontrol('parent',hp,'units','pixels','string','color',...
    'fontsize',8,'backgroundcolor','w','style','text');
climtext.Position(3:4)=climtext.Extent(3:4);
climtext.Position(1:2)= [1 25];

% Color limit table for OD image
climtbl=uitable('parent',hp,'units','pixels','RowName',{},'ColumnName',{},...
    'Data',[0 1000],'ColumnWidth',{40,40},'ColumnEditable',[true true],...
    'CellEditCallback',@climCB);
climtbl.Position(3:4)=climtbl.Extent(3:4);
climtbl.Position(1:2) = [30 25];

% Callback for changing the color limits table
    function climCB(src,evt)
        try
            axImg.CLim=climtbl.Data;
        catch exception
            warning('Bad OD color limits given. Using old value.');
            src.Data(evt.Indices)=evt.PreviousData;
        end
    end

axImg.CLim=climtbl.Data;
%% Image Processing
    function data=processImages(imgs) 
        tin = now;
        data = struct;
        data.Date=datevec(tin);
        data.Name=['ThorCamImage_' datestr(tin,'yyyy-mm-dd_HH-MM-SS')];
        data.X = 1:size(imgs,2);
        data.Y = 1:size(imgs,1);
        data.Images = imgs;
        data.ExposureTime = camera_settings.ExposureTime;

        if size(imgs,3)==2
            data.Data = imgs(:,:,1)-imgs(:,:,2);
        end

        data.Gain_dB = camera_settings.Gain_dB;
        data.Gain = camera_settings.Gain;
        data.PixelSize = camera_settings.PixelSize;
        [data.Params,data.Units,data.Flags]=grabSeqeunceParams; 
end

%% Param
function [vals,units,flags]=grabSeqeunceParams(src)
    if nargin~=1
    src = ['Y:\_communication\control2.mat'];
    end
    data=load(src);
    vals = data.vals;
    units=data.units;
    flags=data.flags;
end
%% Camera Functions
  

% Grab the image camera if available
function img=grabImage
    img=[];
    imageFrame = tlCamera.GetPendingFrameOrNull;
    if ~isempty(imageFrame)
        imageData = imageFrame.ImageData.ImageData_monoOrBGR;
        imageHeight = imageFrame.ImageData.Height_pixels;
        imageWidth = imageFrame.ImageData.Width_pixels;   
        img = reshape(uint16(imageData), [imageWidth, imageHeight]); 
        img = img';   
        %img=double(img);        
    end
end

function updateImage   
    % Grab the image
    img=grabImage;

    % Exit if no image to be had
    if isempty(img)
        return
    end    
    hImg.CData = img;        
    updateDescStr(sum(sum(img)));    
end

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
        tlCamera.Gain=tlCamera.ConvertDecibelsToGain(uint32(settings.Gain_dB));                
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
%         tlCamera.Arm;    
        tlCamera.Disarm;

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

