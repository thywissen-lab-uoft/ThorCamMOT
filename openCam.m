function hFMain=openCam(sn,tlCameraSDK)

% The two cameras that we connect to for fluorescence analysis are:
% 10118 - X CAMERA
% 10148 - Y CAMERA
 
% Default settings for the camera.  These are loaded upon a new connection
% to the camera.
tExp=64;                % Exposure time us
gGain=52;               % Gain in dB
ROIbg=[800 1000 1 200]; % ROI for background detection
pixelsize=3.7;          % the pixelsize in um
mag=1;                  % the magnification atomD=pixelsize/mag
cdata=imresize(imread(['browse.jpg']),[40 40]);

doDebug=0;
%% Open the camera
cam=openCamera(sn);

[ROI,xVec,yVec]=readROI;
imgBG=zeros(length(yVec),length(xVec));

cameraMode='Live';

% Stucture for live settings
live=struct;
live.Fit=false;
live.AutoBackground=false;
live.ImgBackground=imgBG;
live.BackgroundSubtract=false;

% Structure for triggered settings
trig=struct;
trig.Fit=false;
trig.AutoBackground=false;
trig.ImgBackground=imgBG;
trig.BackgroundSubtract=false;
trig.Mode=2; % 0 : one image, 1 : background then image, 2 : image then background
trig.NumImages=0;
trig.Images={};

%% Initialize analysis figures
% Initialize graphics for analysis.  This makes it easier on graphical
% resources and prevents the new figures from taking over.
    function [hF,img1,img2,img3,tbl]=initQuick
        hF=figure(str2num(sn)+1);
        clf
        set(hF,'color','w','toolbar','none','windowstyle','docked','Tag','GUI');

        subplot(221,'parent',hF);
        img1=imagesc(xVec,yVec,zeros(length(yVec),length(xVec)));
        title('image 1');
        axis equal tight; hold on;

        subplot(223,'parent',hF);
        img2=imagesc(xVec,yVec,zeros(length(yVec),length(xVec)));
        title('image 2');
        axis equal tight; hold on;        

        subplot(222,'parent',hF);
        img3=imagesc(xVec,yVec,zeros(length(yVec),length(xVec)));
        title('raw data');
        axis equal tight; hold on   
        
        a=subplot(224);
        pos=a.Position;
        delete(a);
        
       % Table for adjusting gain and exposure
        tbl=uitable(hF,'units','pixels','ColumnName',{},...
            'ColumnEditable',[false],'ColumnWidth',{120},...
            'ColumnFormat',{'char'},...
            'Data',{datestr(now,'yyyy-mm-dd'); datestr(now,'HH:MM:SS'); gGain; tExp; 1},...
            'RowName',{'date','time','gain (dB)','exposure (us)','pixelsize (um)',},...
            'fontsize',12);
        tbl.Position(3:4)=tbl.Extent(3:4);  
        tbl.Units='normalized';
        tbl.Position(1:2)=pos(1:2); 
    end

%% Gaussian Analysis Figure
% Gaussian analysis automatically performed on images. This figure holds
% the graphical objects for this information.

hF3=figure(str2num(sn)+2); % Create the figure.
co=get(gca,'colororder');
clf
set(hF3,'color','w','toolbar','none','windowstyle','docked','Tag','GUI');

tt=uicontrol('style','text','units','pixels','fontsize',12,...
    'backgroundcolor','w','String',datestr(now,'yyyy-mm-dd_HH-MM-SS'));
tt.Position(1:2)=[2 2];
tt.Position(3:4)=tt.Extent(3:4);

% Subplot for the fluorescence image data
subplot(211);
imgGauss=imagesc(xVec,yVec,zeros(length(yVec),length(xVec)));
axis equal tight; hold on;
pRetGauss=plot(0,0,'r-','linewidth',2);

% X axis fit
axX=subplot(425,'box','on','linewidth',1,'fontsize',10);hold on;
pX=plot(0,0,'-','linewidth',1,'color',co(1,:));
pXF=plot(0,0,'-','linewidth',2,'color','r');
text(0,.98,'~$n(x,y_c)$','units','normalized','fontsize',12,...
    'verticalalignment','top','horizontalalignment','left',...
    'interpreter','latex');

% Y axis fit
axY=subplot(426,'box','on','linewidth',1,'fontsize',10);hold on;
pY=plot(0,0,'-','linewidth',1,'color',co(1,:));
pYF=plot(0,0,'-','linewidth',2,'color','r');
text(0,.98,'~$n(y,x_c)$','units','normalized','fontsize',12,...
    'verticalalignment','top','horizontalalignment','left',...
    'interpreter','latex');

% Radial fit
axR=subplot(427,'box','on','linewidth',1,'fontsize',10);hold on;
pR=plot(0,0,'-','linewidth',2,'color',co(1,:));
pRF=plot(0,0,'-','linewidth',2,'color','r');
text(0,.01,'~$n(r)$','units','normalized','fontsize',12,...
    'verticalalignment','bottom','horizontalalignment','left',...
    'interpreter','latex');

% Summary table
d={ ['gain    ' num2str(gGain) ' dB'];...
    ['exp     ' num2str(tExp) ' us'];...
    ['px size ' num2str(pixelsize) ' um'];...
    ['(xC,yC) ' '(0,0)'];...
    ['(xS,yS) ' '(0,0)'];...
    ['nbg     ' '0'];...
    ['Amp     ' '0'];...
    ['Ntot    ' '0']};
gTbl=uitable(hF3,'units','pixels','ColumnName',{},'RowName',{},...
    'ColumnEditable',[false],'ColumnWidth',{200},'FontSize',9,...
    'ColumnFormat',{'char'},'Data',d,'fontname','monospaced');
drawnow;
gTbl.Position(3:4)=gTbl.Extent(3:4);  
drawnow;
gTbl.Units='normalized';
gTbl.Position(1)=(axY.Position(1)+axY.Position(3)/2)-gTbl.Position(3)/2;
gTbl.Position(2)=(axR.Position(2)-gTbl.Position(4)+axR.Position(4));

% Function for updating the gaussian analysis upon receiving new image data
    function dstruct=updateGauss(dstruct)
        % Grab data from the image structure
        x=dstruct.X;    % X vector
        y=dstruct.Y;    % Y vector
        z=dstruct.Data; % N counts
        
        % Update the graphical data
        set(imgGauss,'XData',x,'YData',y,'CData',z);
        set(pX,'XData',x,'YData',sum(z,1));axX.XLim=[min(x) max(x)];
        set(pY,'XData',y,'YData',sum(z,2));axY.XLim=[min(y) max(y)];

        
        % Perform the fit
        fout=gaussfit2D(x,y,z);
        c=coeffvalues(fout);        
        dstruct.GaussFit=fout;
        try
            % Do the fit again but better        
            fout=gaussfit2DFine(x,y,z,fout);    
            c=coeffvalues(fout);        
            dstruct.GaussFit=fout;        
        catch exception
            warning('The more exhaustive 2D gauss fit failed. Using the quick one');
        end
        % Plot gaussian fit reticle as 1/e^2
        t=linspace(0,2*pi,400);                   
        xR=c(2)+2*c(3)*cos(t);
        yR=c(4)+2*c(5)*sin(t);    
        set(pRetGauss,'XData',xR,'YData',yR,'linewidth',2);
  
        % Evaluvate the fit for doing numerical projection
        [xx,yy]=meshgrid(x,y);
        zzF=feval(fout,xx,yy);     

        % Calculate and plot radial residual.
        try  % Use a try because converting to polar requires good fit
            % Region to examine is 5sigma, should compare this number to
            % the imaging sensor edge
            s=round(5*max([fout.Xs fout.Ys]));
            
            % Get distances to edges of images
            x1=fout.Xc-min(dstruct.X);
            x2=max(dstruct.X)-fout.Xc;            
            y1=fout.Yc-min(dstruct.Y);
            y2=max(dstruct.Y)-fout.Yc;        
            
            % Restrict polar fit to lie within the image
            s=min([s x1 x2 y1 y2]);

            % Evenly spaced polarmesh, so not even in cartersian space
            theta=linspace(0,2*pi,360);
            rho=linspace(0,s,floor(s));

            % Convert 1D vectors into a 2D mesh over the disk
            [t,r]=meshgrid(theta,rho);

            % Convert mesh points to cartersian
            [xx,yy]=pol2cart(t,r);
            xx=xx+round(fout.Xc);   % offset x-center by fit
            yy=yy+round(fout.Yc);   % offset y-center by fit

            % Interpolate the gridded fit and raw data.
            % Grid polar mesh is incomensurate with grid cartersian mesh
            vq=interp2(dstruct.X,dstruct.Y,double(dstruct.Data),xx,yy);
            rhoData=sum(vq,2);
            vq=interp2(dstruct.X,dstruct.Y,double(zzF),xx,yy);
            rhoFit=sum(vq,2);        

            % Update radial graphical objercts
            set(pR,'XData',rho,'YData',rhoData);
            set(pRF,'XData',rho,'YData',rhoFit);
            axR.XLim=[0 max(rho)];
            set(axR,'color','w');
            drawnow;
            

            
        catch exception
            % In event of error, do not updte graphics
            disp(['Radial analysis failed! Probably from bad fit']);
            set(pR,'XData',0,'YData',0);
            set(pRF,'XData',0,'YData',0);
            set(axR,'color','r');
        end
        
        % Plot the cartersian data cuts with fits. This can through an
        % error if the fit is not well behaved
        try 
            % Update cartersian plots using a cut
            zXCut=z(floor(fout.Yc),:); % Z(Xc,x)
            zYCut=z(:,floor(fout.Xc)); % Z(x,Yc)            
            set(pX,'XData',x,'YData',zXCut);
            set(pY,'XData',y,'YData',zYCut);

            % Update cartesion plots
            zXCutF=zzF(floor(fout.Yc),:); % Zfit(Xc,x)
            zYCutF=zzF(:,floor(fout.Xc)); % Zfit(x,Yc)      
            set(pXF,'XData',x,'YData',zXCutF);
            set(pYF,'XData',y,'YData',zYCutF);
                        
            x1=max([min(x) fout.Xc-6*fout.Xs]);
            x2=min([max(x) fout.Xc+6*fout.Xs]);
            
            y1=max([min(y) fout.Yc-6*fout.Ys]);
            y2=min([max(y) fout.Yc+6*fout.Ys]);
            
            axX.XLim=[x1 x2];           
            axY.XLim=[y1 y2];
            
        
        catch
            warning('Issue with cartersian graphical plot. Showing the sum');
            % Update cartersian plots using a cut            
            set(pX,'XData',x,'YData',sum(z,1));
            set(pY,'XData',y,'YData',sum(z,2));

            % Update cartesion plots
            set(pXF,'XData',x,'YData',sum(zzF,1));
            set(pYF,'XData',y,'YData',sum(zzF,2));

            axY.XLim=[min(y) max(y)];
            axX.XLim=[min(x) max(x)];
        end
        

        
        % Update summary table and other text objects
        tt.String=datestr(dstruct.Date,'yyyy-mm-dd_HH-MM-SS');
        dd={ ['gain       ' num2str(dstruct.Gain) ' dB'];...
            ['exposure   ' num2str(dstruct.ExposureTime) ' us'];...
            ['pixelsize  ' num2str(dstruct.PixelSize) ' um'];...
            ['xC,yC      ' num2str(round(fout.Xc,1)) ',' num2str(round(fout.Yc,1)) ''];...
            ['xS,yS      ' num2str(round(fout.Xs,1)) ',' num2str(round(fout.Ys,1)) ''];...
            ['nbg        ' num2str(round(fout.nbg,1))];...
            ['Amp        ' num2str(round(fout.A,1))];...
            ['Ntot       ' num2str(2*pi*fout.Xs*fout.Ys*fout.A,3)]};
        gTbl.Data=dd;
        
    end
    %%

[hFRaw,himg1,himg2,himg3,tbl]=initQuick;
%% Create Graphical Interface

% Initialize the figure GUI
hFMain=figure(str2num(sn));
set(hFMain,'Color','w','Toolbar','None','CloseRequestFcn',@closeCB,...
    'WindowStyle','docked','Tag','GUI');
clf
if isequal(sn,'10118')
    hFMain.Name=['MOT CAMERA SN - ' sn ' X'];
else
    hFMain.Name=['MOT CAMERA SN - ' sn ' Y'];
end

% Callback function for closing the camera GUI
    function closeCB(fig,~)        
        stop(timerLive); stop(timerTrig);       % stop timers
        pause(1);                               % Wait
        delete(timerLive); delete(timerTrig);   % delete timer
        closeCamera(cam);                       % disonnect camera 
        disp('done');    
        delete(fig)                             % close figure
    end


%%%%%%%%%%%%%%%% Initialize image and axis %%%%%%%%%%%%%%%%%%%%%%%
ax=axes;
cla
hImg=imagesc(xVec,yVec,imgBG);
set(ax,'XAxisLocation','top','fontsize',14,'fontname','arial',...
    'CLim',[0 1024]);
xlabel('x pixels');ylabel('y pixels');
axis equal tight
colormap parula
cbar=colorbar;
cbar.Label.String='counts';
hold on

% Initialize the fit reticle
pRet=plot(0,0,'r-','Visible','off');
pp=[];

% Text objects at bottom of axis for summary of settings
textCounts=text(2,-2,'test','units','pixels','verticalalignment','top',...
    'color','k','fontweight','bold','fontsize',12);
textExp=text(100,-2,[num2str(tExp) ' us'],'units','pixels','color','k',...
    'verticalalignment','top','fontweight','bold','fontsize',12);
textGain=text(170,-2,[num2str(gGain) ' dB'],'units','pixels','color','k',...
    'verticalalignment','top','fontweight','bold','fontsize',12);
textFit=text(2,2,'boop','units','pixels','verticalalignment','bottom',...
    'color','r','fontweight','bold','fontsize',12,'visible','off');

%%%%%%%%%%%%%%%% Graphics for settings %%%%%%%%%%%%%%%%%%%%%%%

% Radio button group for operation mode
bg = uibuttongroup('units','pixels','backgroundcolor','w',...
    'position',[0 0 80 40],'SelectionChangedFcn',@chCameraMode);        
% Create three radio buttons in the button group.
b1=uicontrol(bg,'Style','radiobutton','String','Live',...
    'Position',[0 0 80 20],'units','pixels','backgroundcolor','w');
b2=uicontrol(bg,'Style','radiobutton','String','Triggered',...
    'Position',[0 20 80 20],'units','pixels','backgroundcolor','w');

% Button for start/stop
bEnable=uicontrol('style','pushbutton','String','start',...
    'enable','on','units','pixels','fontsize',8,'Callback',@buttCB,...
    'backgroundcolor',[80 200 120]/255);
bEnable.Position(1)=bg.Position(3);
bEnable.Position(2)=0;
bEnable.Position(3:4)=[60 20];

% Button for start/stop
bClear=uicontrol('style','pushbutton','String','clear',...
    'enable','on','units','pixels','fontsize',8,'Callback',@clearCB);
bClear.Position(1)=bg.Position(3);
bClear.Position(2)=20;
bClear.Position(3:4)=[60 20];

% Button for settings GUI
bSettings=uicontrol('style','pushbutton','String','settings',...
    'enable','on','units','pixels','fontsize',8,'Callback',@settingsGUI);
bSettings.Position(1)=bEnable.Position(1)+bEnable.Position(3);
bSettings.Position(2)=0;
bSettings.Position(3:4)=[60 20];

% Checkmark for saving
cSave=uicontrol('style','checkbox','String','autosave?','units','pixels',...
    'value',0,'fontsize',8,'callback',@saveCheck,'backgroundcolor','w');
cSave.Position=bSettings.Position;
cSave.Position(3)=70;
cSave.Position(1)=cSave.Position(1)+bSettings.Position(3);

% Button for file selection of the sequenece file
strr='Choose a directory to save the recorded data';
bBrowse=uicontrol('style','pushbutton','CData',cdata,...
    'backgroundcolor','w','tooltipstring',strr,'enable','off');
bBrowse.Position(3:4)=size(cdata,[1 2]);
bBrowse.Position(1)=cSave.Position(1);
bBrowse.Position(2)=cSave.Position(4)+1;
bBrowse.Callback=@browseCB;

% String for current save directory
tSaveDir=uicontrol('style','text','string','directory','backgroundcolor','w',...
    'fontsize',8,'units','pixels','horizontalalignment','left',...
    'enable','off');
tSaveDir.Position(3)=500;
tSaveDir.Position(4)=bBrowse.Position(4)-3;
tSaveDir.Position(1)=bBrowse.Position(1)+bBrowse.Position(3)+2;
tSaveDir.Position(2)=bBrowse.Position(2);

    function saveCheck(event,~)
        if event.Value
            tSaveDir.Enable='on';
            bBrowse.Enable='on';
        else
            tSaveDir.Enable='off';
            bBrowse.Enable='off';
        end
    end

    function browseCB(~,~)
        str=getDayDir;
        str=uigetdir(str);
        if str
            if isequal(cameraMode,'Triggered') && ~isequal(tSaveDir.String,str)
                clearBuffer;
            end
            tSaveDir.String=str;
        else
            disp('no directory chosen!');
        end
    end

    function buttCB(~,~)
        if isequal(cameraMode,'Live')
            if isequal(timerLive.Running,'on')
                stop(timerLive);
                bEnable.String='start';
                bEnable.BackgroundColor=[80 200 120]/255;
                b1.Enable='on';
                b2.Enable='on';
                disp('stopping live acquisition');
           else
                start(timerLive);
                bEnable.String='stop';
                bEnable.BackgroundColor=[255 102 120]/255;
                b1.Enable='off';
                b2.Enable='off';
                disp('starting live acquisition');
           end
        end
        if isequal(cameraMode,'Triggered')
            if isequal(timerTrig.Running,'on')
                stop(timerTrig)
                bEnable.String='start';
                bEnable.BackgroundColor=[80 200 120]/255;
                b1.Enable='on';
                b2.Enable='on';
                disp('stopping triggered acquisition');
                clearBuffer;
            else
                clearBuffer;
                disp('starting triggered acquisition');
                start(timerTrig)
                bEnable.String='stop';
                bEnable.BackgroundColor=[255 102 120]/255;
                b1.Enable='off';
                b2.Enable='off';
            end           
        end        
    end

    function clearBuffer
        disp('Clearing trigger buffer.')
        imageFrame = cam.GetPendingFrameOrNull; 
        trig.NumImages=0;
        trig.Images={};
        if ~isempty(imageFrame)
           clearBuffer;
        end       
    end

    function clearCB(~,~)
        if isequal(cameraMode,'Triggered')
            clearBuffer;
        end

        % Clear stored live data
        T=[];
        Y=[];
        t0=now;
        cBGs=[];
    end
        
% Callback function for when the mode of the camera is changed
    function chCameraMode(~,event)                
        str=event.NewValue.String;
        switch str
            % Switch the camera to live mode
            case 'Live'     
                disp('Switching camera to live mode');
                cameraMode='Live';
                pause(0.5);                 % Wait
                hImg.CData=hImg.CData*0;    % Clear the display image
                textCounts.String='0';      % Clear the display counts
                cam.Disarm;                 % Stop the camera
                cam.OperationMode=...       % Change to software trigger
                    Thorlabs.TSI.TLCameraInterfaces.OperationMode.SoftwareTriggered;
                cam.Arm;                    % Start the camera
                trig.NumImages=0;
                trig.Images={};
                bEnable.String='Start';
            % Switch the camera to triggered mode
            case 'Triggered'
                cameraMode='Triggered';
                disp('Switching camera to triggered mode');
                pause(0.5);                 % Wait
                hImg.CData=hImg.CData*0;    % Clear the display image
                textCounts.String='0';      % Clear the display counts
                cam.Disarm;                 % Stop the camera
                cam.OperationMode=...       % Change to hardware trigger
                    Thorlabs.TSI.TLCameraInterfaces.OperationMode.HardwareTriggered;
                cam.Arm;                    % Start the camera
                trig.NumImages=0;
                trig.Images={};
        end      
    end

%% Timer Objects

% Timers for updates (this should change to nly one timer)
timerTrig=timer('Name','TrigChecker','executionmode','fixedspacing',...
    'period',0.1,'TimerFcn',@trigCB);
timerLive=timer('Name','liveupdate','executionmode','fixedspacing',...
    'period',0.001,'TimerFcn',@liveCB);

t0=now;
T=[];
Y=[];
cBGs=[];

%% Settings callback
    function settingsGUI(~,~)
        % Read the camera settings
        gGain=double(cam.ConvertGainToDecibels(cam.Gain));
        tExp=double(cam.ExposureTime_us);
        
        % Store the camera settings
        [ROI,xVec,yVec]=readROI;
        textExp.String=[num2str(tExp) ' \mus'];   
        textGain.String=[num2str(gGain) ' dB'];        
        
        % Make new figure for settings
        str=[hFMain.Name ' Settings'];
        hFSet=figure('Name',str,'Toolbar','none','menubar','none',...
            'resize','off','color','w');
        hFSet.Position(3:4)=[600 300];
        hFSet.Position(1:2)=hFMain.Position(1:2)+hFMain.Position(3:4)/2-...
            hFSet.Position(3:4)/2;        
        
        %%%%%%%%%%%%%%%%%%%% Global settings panel %%%%%%%%%%%%%%%%%%%%%%%%
        % Panel for adjusting main settings
        hpMain=uipanel('parent',hFSet,'units','pixels','title',...
            'global settings','position',[10 10 250 280],...
            'backgroundcolor','w');  
        
        % Table for adjusting gain and exposure
        tblAcq=uitable(hpMain,'units','pixels','ColumnName',{},...
            'ColumnEditable',[true true],'CellEditCallback',@chSet,...
            'Data',[gGain; tExp],'RowName',{'gain (dB)','exposure (us)'},...
            'ColumnWidth',{50 50});
        tblAcq.Position(3:4)=tblAcq.Extent(3:4);        
        tblAcq.Position(1:2)=[10 10];        
        
        % Table for adjusting hardware ROI
        tblROI=uitable(hpMain,'units','pixels','RowName',{},'Data',ROI,...
            'ColumnEditable',[true true true true],'CellEditCallback',@chROI,...
            'ColumnName',{'x1','x2','y1','y2'},'ColumnWidth',{50 50 50 50});
        tblROI.Position(3:4)=tblROI.Extent(3:4);        
        tblROI.Position(1:2)=[10 70];
           
        % Table for adjusting color limits on plot
        tblCLIM=uitable(hpMain,'units','pixels','RowName',{},...
            'Data',ax.CLim,'ColumnWidth',{50 50},'ColumnName',{'c1','c2'},...
            'ColumnEditable',[true true],'CellEditCallback',@chCLIM); 
        tblCLIM.Position(3:4)=tblCLIM.Extent(3:4);        
        tblCLIM.Position(1:2)=[10 120];   
        
        % Checkbox for debug mode
        uicontrol('parent',hpMain,'style','checkbox','units','pixels',...
            'string','debug mode','callback',@cDebugCB,...
            'Position',[10 200 150 20],'backgroundcolor','w',...
            'Value',doDebug);
        
        % Callback for editing debug mode
        function cDebugCB(cb,~)
           if cb.Value
               doDebug=1;
           else
               doDebug=0;
           end
        end
        
        %%%%%%%%%%%%%%%%%%%% Triggered settings panel %%%%%%%%%%%%%%%%%%%%%%
        hpTrig=uipanel(hFSet,'units','pixels','title','trigger settings',...
            'position',[270 10 320 140],'backgroundcolor','w');        

        % Checkbox for automated fitting
%         cFitTrig=uicontrol('parent',hpTrig,'style','checkbox','string','fit?',...
%             'units','pixels','position',[220 10 50 20],...
%             'backgroundcolor','w');
        
        % Radio button group for operation mode
        bgTrig = uibuttongroup('parent',hpTrig,'units','pixels','backgroundcolor','w',...
            'position',[10 10 160 80],'SelectionChangedFcn',@chTrigMode,...
            'Title','Trigger Mode');
        % Create three radio buttons in the button group.
        a=uicontrol(bgTrig,'Style','radiobutton','String','image only',...
            'Position',[0 0 160 20],'units','pixels',...
            'backgroundcolor','w','UserData',0);
        b=uicontrol(bgTrig,'Style','radiobutton','String','background then image',...
            'Position',[0 20 160 20],'units','pixels',...
            'backgroundcolor','w','UserData',1);
        c=uicontrol(bgTrig,'Style','radiobutton','String','image then background',...
            'Position',[0 40 160 20],'units','pixels',...
            'backgroundcolor','w','UserData',2);  
        
        switch trig.Mode
            case 0
                a.Value=1;
            case 1
                b.Value=1;
            case 2
                c.Value=1;
        end
        
        function chTrigMode(~,b)
            disp(['Changing triggered mode to ' b.NewValue.String]);
            trig.Mode=b.NewValue.UserData;           
        end
        
        uicontrol('parent',hpTrig,'style','pushbutton','string','Clear Trigger Buffer',...
            'units','pixels','position',[10 100 120 20],...
            'callback',@trigResetCB);   
        
        function trigResetCB(~,~)
            disp('Clearing stored image buffer.');
            trig.NumImages=0;
            trig.Images={};
        end
        
        %%%%%%%%%%%%%%%%%%%% Live mode settings panel %%%%%%%%%%%%%%%%%%%%%%
        hpLive=uipanel(hFSet,'units','pixels','title','Live Mode Settings',...
            'position',[270 150 320 140],'backgroundcolor','w');  
        
        % Pushbutton for viewing the live number of counts
        uicontrol('parent',hpLive,'style','pushbutton','string',...
            'open live counts','units',...
            'pixels','position',[10 10 100 20],'Callback',@bPDCB);
        
        % Checkbox for auto background search
        uicontrol('parent',hpLive,'style','checkbox','units','pixels',...
            'string','auto-background search','callback',@cBkgdCB,...
            'Position',[10 30 150 20],'backgroundcolor','w',...
            'value',live.AutoBackground);
        
        % Checkbox for subtracting the background
        cSubLive=uicontrol('parent',hpLive,'style','checkbox','units','pixels',...
            'string','background subtract','position',[10 50 150 20],...
            'backgroundcolor','w','callback',@cSubCB,...
            'Value',live.BackgroundSubtract);
        
        % Checkbox for automated fitting
        cFitLive=uicontrol('parent',hpLive,'style','checkbox','string','fit?',...
            'units','pixels','position',[10 70 50 20],'Callback',@cFitCB,...
            'backgroundcolor','w','Value',live.Fit);
        
        % Callback for check box on background subtract
        function cSubCB(cb,~)
            if cb.Value
                live.BackgroundSubtract=true;
            else
                live.BackgroundSubtract=false;
            end
        end
        
        % Callback for check box on auto background search
        function cBkgdCB(cb,~)
            if cb.Value
                live.AutoBackground=true;
                imgBG=1024*ones(5000,5000);
                live.BackgroundSubtract=false;
                set(cSubLive,'Value',0,'Enable','off');
                set(cFitLive,'value',0,'enable','off');
            else
                live.AutoBackground=false;
                set(cSubLive,'Enable','on');
                set(cFitLive,'enable','on');
            end
        end

        % Callback for check box on engaging the fit
        function cFitCB(cb,~)
            if cb.Value
                live.Fit=true;
                pRet.Visible='on';
                textFit.Visible='on';
            else
                live.Fit=false;
                pRet.Visible='off';
                textFit.Visible='off';
            end
        end   

        function bPDCB(~,~)
           hFPD=figure('Name','photodiode','toolbar','none','menubar','none');
           hFPD.Color='w';
           axes;
           pp=plot(0,0);
           xlabel('time (s)');
           ylabel('counts');
           set(gca,'FontSize',14);       
           uicontrol('style','pushbutton','string','clear data',...
               'units','pixels','position',[0 0 70 20],'callback',@bResetCB);       
            function bResetCB(~,~)
               T=[];
               Y=[];
               t0=now;
               cBGs=[];
            end       
        end

        function chROI(tbl,~)
            disp(['Changing the hardware ROI. This needs to stop the ' ...
                'camera and the software timer.']);
            stop(timerLive);
            pause(0.5);
            newROI=tbl.Data;     
            setROI(newROI);          
            set(hImg,'XData',xVec,'YData',yVec,'CData',...
                zeros(length(yVec),length(xVec)));
            start(timerLive);
        end
                
        function chSet(tbl,data)
            r=data.Indices(1);
            val=data.NewData;            
            % Gain goes for 0 to 48 dB
            % Exposure goes from 64us to 51925252us
            
            switch r
                case 1                    
                    if val>=0 && val<=48
                        val=round(val,1);
                        disp(['Changing gain to ' num2str(val) ' dB']);
                        gVal=cam.ConvertDecibelsToGain(val);
                        cam.Gain=gVal;                        
                        gGain=cam.ConvertGainToDecibels(cam.Gain);
                        tbl.Data(data.Indices)=gGain;                        
                        textGain.String=[num2str(gGain) ' dB'];
                    else
                        tbl.Data(data.Indices)=data.PreviousData;                        
                    end            
                case 2
                    if val>=64 && val<=1E5
                        val=round(val);
                        disp(['Changing exposure to ' num2str(val) ' us']);
                        cam.ExposureTime_us=uint32(val);
                        tbl.Data(r)=double(cam.ExposureTime_us);
                        tExp=tbl.Data(r);
                        textExp.String=[num2str(tExp) ' \mus'];                        
                    else
                        tbl.Data(r)=data.PreviousData;                        
                    end
            end          
        end
        
        function chCLIM(tbl,data)
            cThis=data.Indices(2);            
            cOther=mod(cThis,2)+1;            
            val=data.NewData;            
            if cThis==1
                if val>=0 && val<=tbl.Data(cOther)
                   ax.CLim=[val tbl.Data(cOther)]; 
                else
                    tbl.Data(cThis)=data.PreviousData;
                end
            else
                if val<=1024 && val>=tbl.Data(cOther)
                   ax.CLim=[tbl.Data(cOther) val]; 
                else
                    tbl.Data(cThis)=data.PreviousData;
                end                
            end
        end
    end


%% Functions

% Callback function for checking if the camera was triggered
tDebug=now;
    function trigCB(~,~)
        if doDebug       
            dT=(now-tDebug)*24*60*60;
            
            if dT>5          
                % Grab image
                img=grabImage;
                disp([datestr(now,13) ' Trigger (' num2str(trig.NumImages+1) ')']); 
                trig.Images{1}=img;
                trig.NumImages=1;
                
                % Grab image
                img=grabImage;
                disp([datestr(now,13) ' Trigger (' num2str(trig.NumImages+1) ')']); 
                trig.Images{2}=img;
                
                % Process
                processTriggeredImages;
                
                % reset
                trig.NumImages=0;
                trig.Images={}; 
                
                tDebug=now;
            end
        end
                
        if cam.NumberOfQueuedFrames && ~doDebug
            img=grabImage;
            disp([datestr(now,13) ' Trigger (' num2str(trig.NumImages+1) ')']); 
            switch trig.NumImages
                case 0
                    trig.Images{1}=img;
                    trig.NumImages=1;
                case 1
                    trig.Images{2}=img;
                    processTriggeredImages;
                    trig.NumImages=0;
                    trig.Images={};
                otherwise
                    disp('hinonono');
                    % Clear triggered data
                    trig.NumImages=0;
                    trig.Images={};
            end
        end
                
    end

% Process the triggered images
    function processTriggeredImages     
        % Create the subtracted image data
        data=(double(trig.Images{1})-double(trig.Images{2}))*(-1)^(trig.Mode); 
        data=int16(data);
        % Assign the data and metadata to the structure
        dstruct=struct;
        dstruct.Date=datevec(now);
        dstruct.Data=data;dstruct.X=xVec;dstruct.Y=yVec;
        dstruct.ExposureTime=tExp;
        dstruct.Gain=gGain;    
        dstruct.ROI=ROI;
        dstruct.PixelSize=pixelsize;
        
        % Raw Data Figure
        set(himg1,'XData',xVec,'YData',yVec,'CData',trig.Images{1});
        set(himg2,'XData',xVec,'YData',yVec,'CData',trig.Images{2});
        set(himg3,'XData',xVec,'YData',yVec,'CData',data);
        
        tbl.Data{1}=datestr(dstruct.Date,'yyyy-mm-dd');
        tbl.Data{2}=datestr(dstruct.Date,'HH:MM:SS');
        tbl.Data{3}=gGain;
        tbl.Data{4}=tExp;
        tbl.Data{5}=pixelsize;
        
        hFRaw.Color=[220 237 200]/255;        
        drawnow;
        disp([datestr(now,13) ' New image!']);
        pause(0.15);
        hFRaw.Color='w';
        
        % Perform the default gaussian fit
        dstruct=updateGauss(dstruct);
        
        % Grab and append the sequence parameters
        out=grabSequenceParams;
        dstruct.Params=out;
        disp(out);
        
%         % Write to results (this will change in the future
%         str='Rb_molasses_det';
%         cvals=coeffvalues(dstruct.GaussFit);
%         cvals = [dstruct.Params.(str), cvals];   
%         try
%             a=['results(end+1,:)=[' num2str(cvals) '];'];
%             evalin('base',a)
%         catch exception
%             a=['results=[' num2str(cvals) '];'];
%             evalin('base',a)
%         end
%         
%         % Save the raw image (this will change in the future)
%         data=double(dstruct.Data);
%          save_data = 0;
%         if (save_data)
%         dir = 'Y:\Data\2020\2020.09\29 September 2020\G_Rb_molasses_TOF_40ms_vary_molasses_time_Cam10148';
%         save(fullfile(dir,[str,'_',num2str(dstruct.Params.(str)),'_',datestr(now,'HHMMss'),'.mat']),'data');
%         end       
%         
        
        % Save image and data to file
        if cSave.Value
            fname=[datestr(dstruct.Date,'yyyy-mm-dd_HH-MM-SS_') sn];  
            str=fullfile(tSaveDir.String,fname);
            disp(['Saving to ' str]);         
            save(str,'dstruct');
        end

        %% Clear triggered data
        trig.NumImages=0;
        trig.Images={};
        
    end

% Callback function for live update
    function liveCB(~,~)
        cam.IssueSoftwareTrigger;        
        pause(0.02);
        updateImage;
    end

% Grab the image camera if available
    function img=grabImage
        img=[];
        imageFrame = cam.GetPendingFrameOrNull;
        if ~isempty(imageFrame)
            imageData = imageFrame.ImageData.ImageData_monoOrBGR;
            imageHeight = imageFrame.ImageData.Height_pixels;
            imageWidth = imageFrame.ImageData.Width_pixels;   
            img = reshape(uint16(imageData), [imageWidth, imageHeight]);  
            img = img';
            if isequal(sn,'10148')
                img = imrotate(img,180);       
            end
            img=double(img);
        end
        
        if doDebug && isequal(cameraMode,'Live')
           a=datevec(now);
           a=a(6);           
           t=mod(a,8);
           N0=800*(1+rand*.05)*(1-exp(-t/2));  
           xC=mean(xVec)+rand*10;
           yC=mean(yVec)+rand*10;     
           yS=100*(1+rand*.2);
           xS=200*(1+rand*.2);
           [xx,yy]=meshgrid(xVec,yVec);           
           foo=@(x,y) N0*exp(-(x-xC).^2/(2*xS^2)).*exp(-(y-yC).^2/(2*yS^2));
           data=foo(xx,yy);          
           noise=50*rand(length(yVec),length(xVec));           
           img=data+noise;             
        end
          
        if doDebug && isequal(cameraMode,'Triggered')
            N0=400*(1+rand*.05);  
            xC=mean(xVec)+rand*10;
            yC=mean(yVec)+rand*10;     
            yS=100*(1+rand*.05);
            xS=200*(1+rand*.05);
            [xx,yy]=meshgrid(xVec,yVec);           
            foo=@(x,y) N0*exp(-(x-xC).^2/(2*xS^2)).*exp(-(y-yC).^2/(2*yS^2));
            data=foo(xx,yy);          
            noise=300*rand(length(yVec),length(xVec));  
            switch trig.Mode
                case 0
                    img=data;
                case 1
                    if trig.NumImages==0
                       img=noise;
                    else
                        img=data+noise;
                    end
                case 2
                    if trig.NumImages==0
                       img=data+noise;
                    else
                        img=noise;
                    end                    
            end      
        end
    end

    function updateImage   
        % Grab the image
        img=grabImage;
        
        % Exit if no image to be had
        if isempty(img)
            return
        end
           
        % Subtract the background image
        if live.BackgroundSubtract
            hImg.CData=img-imgBG;
        else  
            hImg.CData=img;  
        end        
        c=sum(sum(img));  
        textCounts.String=sprintf('%.4e',c);             

        if live.Fit
            if live.BackgroundSubtract
                data=img-imgBG;
            else
                data=img;
            end
            
           fout=gaussfit2D(xVec,yVec,data);
           cvals=coeffvalues(fout);
           cvals(2:5)=cvals(2:5);
           textFit.String=['cen : (' num2str(round(cvals(2))) ',' ...
               num2str(round(cvals(4))) '), \sigma : (' ...
               num2str(round(cvals(3))) ',' num2str(round(cvals(5))) ')'];
           tVec=linspace(0,2*pi,200);                   
           xvec=cvals(2)+2*cvals(3)*cos(tVec);
           yvec=cvals(4)+2*cvals(5)*sin(tVec);                   
           set(pRet,'XData',xvec,'YData',yvec);      
        end

        t=(now-t0)*24*60*60;
        c=sum(sum(img));                
        T(end+1)=t;Y(end+1)=c;   
        cBG=sum(sum(img(ROIbg)));
        cBGs(end+1)=cBG;

        if length(T)>3E4
            T=[];
            Y=[];
            cBGs=[];
        end

        if live.AutoBackground                            
            if cBG>.5*max(cBGs) && c<sum(sum(imgBG))
                figure(12)
                clf
                imgBG=img;                        
                imagesc(imgBG);
            end
        end               

        try
            set(pp,'XData',T,'YData',Y);
        end
     
    end

%% Camera Functions

% Open camera with a given serial number, default to software trigger
    function tlCamera=openCamera(SN)     
        disp(['Opening camera ' SN]);
        tlCamera = tlCameraSDK.OpenCamera(SN, false);        
        % Get and Set camera parameters
        tlCamera.ExposureTime_us = uint32(tExp);        
        % The default black level should be zero
        tlCamera.BlackLevel=uint32(0);       
        % Set the default gain level to zero
        tlCamera.Gain=tlCamera.ConvertDecibelsToGain(uint32(gGain));                
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
    end

% Close the camera
    function closeCamera(tlCamera)
        disp('Closing camera');
        if (tlCamera.IsArmed)
            tlCamera.Disarm;
        end
        tlCamera.Dispose;        
        delete(tlCamera);           
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

  
function fout=gaussfit2D(Dx,Dy,data)
% Make and X and Y pixel vectors

data=double(data);              % data is double and not uint32 
data=imresize(data,0.25);
Dx=imresize(Dx,.25);
Dy=imresize(Dy,.25);


dSmooth=imgaussfilt(data,15);   % Smooth data
N0=max(max(dSmooth));           % Extract peak

% Get rid of noise
Z=dSmooth;
Z(dSmooth<N0*.5)=0;



% Get the profiles
X=sum(Z,1);
Y=sum(Z,2)';

% Get the total number of counts
Nx=sum(X);
Ny=sum(Y);

% Find the Center
Xc=mean(Dx(X>.9*max(X)));
Yc=mean(Dy(Y>.9*max(Y)));

% Calculate sigma in X and Y
Xs=1.5*sqrt(sum((Dx-Xc).^2.*X)/Nx);
Ys=1.5*sqrt(sum((Dy-Yc).^2.*Y)/Ny);

% Make a mesh grid for fitting
[xx,yy]=meshgrid(Dx,Dy);

% Make an initial guess
Zguess=N0*exp(-(xx-Xc).^2./(2*Xs)^2).*exp(-(yy-Yc).^2./(2*Ys)^2);

% Copy the data
data2=data;
xx2=xx;
yy2=yy;

% Elminate data points below a threshold to reduce fitting space
xx2(Zguess<.15*N0)=[];
yy2(Zguess<.15*N0)=[];
data2(Zguess<.15*N0)=[];

% Calculate the appropriate background
bg=sum(sum(data-Zguess))/(length(X)*length(Y));

% Create fit object
myfit=fittype('A*exp(-(xx-Xc).^2./(2*Xs^2)).*exp(-(yy-Yc).^2./(2*Ys^2))+nbg',...
    'independent',{'xx','yy'},'coefficients',{'A','Xc','Xs','Yc','Ys','nbg'});
opt=fitoptions(myfit);
opt.StartPoint=[N0 Xc Xs Yc Ys bg];
opt.Lower=[N0/10 10 1 10 1 0];
opt.Upper=[5*N0 max(X) range(X) max(Y) range(Y) N0];
opt.Weights=[];

% Check that the upper and lower bounds make sense
badInds=opt.Upper<opt.Lower;
if sum(badInds)
    warning(['Generated lower bounds for gaussian fit exceed the upper ' ...
        'bounds for ' num2str(sum(badInds)) ' parameters. This ' ...
        'may be caused by no atoms.']);
    opt.Lower=[0 0 0 0 0 0];
    opt.Upper=[];
    opt.StartPoint=[100 mean(Dx) range(Dx)/10 mean(Dy) range(Dy)/10 ...
        10];
end



% Perform the fit
[fout,gof,output]=fit([xx2(:) yy2(:)],data2(:),myfit,opt);
% keyboard


end


  
function fout=gaussfit2DFine(Dx,Dy,data,fin)
% Given a pretty good initial gaussian fit, this does a finer fit
% Make and X and Y pixel vectors
data=double(data);              % data is double and not uint32 
dSmooth=imgaussfilt(data,5);   % Smooth data
data2=dSmooth;

data2=imresize(data,0.5);
Dx=imresize(Dx,.5);
Dy=imresize(Dy,.5);

% Only analyze data with 4*sigma of the original fit

% Make a mesh grid for fitting
[xx,yy]=meshgrid(Dx,Dy);
xx=xx(:);
yy=yy(:);
data2=data2(:);


% Remove data points beyond a 3 sigma radius from X center
inds=abs(xx-fin.Xc)>(2.5*fin.Xs);
xx(inds)=[];
yy(inds)=[];
data2(inds)=[];

% Remove data points beyond a 3 sigma radius from Y center
inds=abs(yy-fin.Yc)>(2.5*fin.Ys);
xx(inds)=[];
yy(inds)=[];
data2(inds)=[];

% Create fit object
myfit=fittype('A*exp(-(xx-Xc).^2./(2*Xs^2)).*exp(-(yy-Yc).^2./(2*Ys^2))+nbg',...
    'independent',{'xx','yy'},'coefficients',{'A','Xc','Xs','Yc','Ys','nbg'});
opt=fitoptions(myfit);
opt.StartPoint=coeffvalues(fin);
opt.Lower=[fin.A*.1 max([fin.Xc-50 1]) fin.Xs/5 max([fin.Yc-50 1]) fin.Ys/5 0];
opt.Upper=[fin.A*10 min([fin.Xc+50 max(Dx)]) fin.Xs*5 min([fin.Yc+50 max(Dy)]) fin.Ys*5 fin.A/10];

opt.Weights=[];

% Perform the fit
[fout,gof,output]=fit([xx yy],data2,myfit,opt);


end

function out=grabSequenceParams(src)
if nargin~=1
    src='Y:\_communication\control.txt';
end
disp(['Opening information from from ' src]);
  
out=struct;
% Open the control file
[fid,errmsg] = fopen(src,'rt');
if ~isempty(errmsg)
   warning('Unable to read control.txt. Aborting association'); 
   return
end

% Read the first six lines (and throw them away)
for i = 1:6
    fgetl(fid);
end

% Read the parameters (each line looks like "k_cMOT_detuning: 5")
params = textscan(fid,'%[^:] %*s %s');

% Close the file
fclose(fid);    

% Convert the string into a structure
out=cell2struct(num2cell(str2double(params{2})),params{1});
end



function s3=getDayDir
t=now;

d=['Y:\Data'];
s1=datestr(t,'yyyy');
s2=datestr(t,'yyyy.mm');
s3=datestr(t,'mm.dd');

s1=[d filesep s1];
s2=[s1 filesep s2];
s3=[s2 filesep s3];

if ~exist(s1,'dir')
    mkdir(s1);
end

if ~exist(s2,'dir')
    mkdir(s2);
end

if ~exist(s3,'dir')
    mkdir(s3);
end
end
