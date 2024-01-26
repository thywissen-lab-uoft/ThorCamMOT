function [hF,hF2]=showCenter(data,xVar,opts)

if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end

%% Default options
if nargin==2
    opts=struct;
    opts.CenterSineFit = 0;
    opts.CenterDecaySineFit = 0;
    opts.CenterParabolaFit = 0;
    opts.CenterLinearFit = 0;
    opts.angleTrack = 0;
end

%% Get Data
Xc = data.Xc;
Yc = data.Yc;

params = [data.Params];
xvals = [params.(xVar)];

PixelSize = data.PixelSize;


%% Make Figure

hF=figure('Name',[pad([data.FitType ' centre'],20) FigLabel],...
    'units','pixels','color','w','numbertitle','off');
hF.Position(1)=510;
hF.Position(2)=380;
hF.Position(3)=500;
hF.Position(4)=600;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

resizeFig(hF,t)

%% Track X
hax1=subplot(221);
set(hax1,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');
for nn=1:size(Xc,2)
        
    if median(Yc(:,nn))>1092
        m = 's';     
    else
        m = 'o';   
    end
    
    
   plot(xvals,Xc(:,nn),m,'color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

str='X centre (px)';
text(0.02,.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');


%% Table X

a=subplot(222);
pos=get(a,'position');
delete(a)

sTblX=uitable('FontSize',8,'RowName',{},'ColumnName',{},...
    'ColumnEditable',[false false],'units','normalized');
sTblX.ColumnWidth={100 60};
sTblX.Position=pos;

r = max(Xc(:,nn)) - min(Xc(:,nn));
sTblX.Data={[char(0x0394) 'X (px)'],num2str(round(r,1))};
drawnow;


%% Track Y

hax2=subplot(223);
set(hax2,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');

yyaxis left
set(gca,'YColor','k');

yyaxis right
set(gca,'YColor','k');
yL = 1024;
yH = 1;
for nn=1:size(Yc,2)
    
    if median(Yc(:,nn))>1092
        yyaxis right
        m = 's';
        
        yL = min([yL min(Yc(:,nn))-1024]);
        yH = max([yH max(Yc(:,nn))-1024]);

    else
        yyaxis left        
        m = 'o';
        yL = min([yL min(Yc(:,nn))]);
        yH = max([yH max(Yc(:,nn))]);

    end
    
   plot(xvals,Yc(:,nn),m,'color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end

yLCam = [yL yH]+(yH-yL)*.1*[-1 1];

yyaxis left
ylim(yLCam);

yyaxis right
ylim(yLCam+1024);

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

str='Y centre (px)';
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');

%% Table Y

a=subplot(224);
pos=get(a,'position');
delete(a)

sTblY=uitable('FontSize',8,'RowName',{},'ColumnName',{},...
    'ColumnEditable',[false false],'units','normalized');
sTblY.ColumnWidth={100 60};
sTblY.Position=pos;
 r = max(Yc(:,nn)) - min(Yc(:,nn));
sTblY.Data={[char(0x0394) 'Y (px)'],num2str(round(r,1))};
drawnow;

%% Fits
    tbl_data={};


if isequal(opts.CenterDecaySineFit,1) && length(xvals)>4
    tVec=linspace(min(xvals),max(xvals),100);    
    
    % X Fit
    axes(hax1);
    fit1=makeSineDecayFit(xvals',Xc(:,nn));
    plot(tVec,feval(fit1,tVec),'r-');  

    sTblX.ColumnWidth={60 60 60};
    sTblY.ColumnWidth={60 60 60};

    cX=coeffvalues(fit1);
    cIntX=confint(fit1);
    
    
    %tbl_data{1,3}=range(cInt(:,1))/2;
    %tbl_data{2,3}=range(cInt(:,2))/2;

    %tbl_data{3,3}=1./(range(cInt(:,2))/2);
    %tbl_data{4,3}=range(cInt(:,3))/2;
    %tbl_data{5,3}=range(cInt(:,4))/2;
    %tbl_data{6,3}=range(cInt(:,5))/2;

    
    sTblX.tbl_data={};
    tbl_data{1,1}='amp (px)';
    tbl_data{2,1}='period';
    tbl_data{3,1}='freq';

    tbl_data{4,1}='phase (rad)';
    tbl_data{5,1}='offset (px)';
    tbl_data{6,1}='tau ';

    tbl_data{1,2}=cX(1);
    tbl_data{2,2}=cX(2);
    tbl_data{3,2}=1/cX(2);

    tbl_data{4,2}=cX(3);
    tbl_data{5,2}=cX(4);
    tbl_data{6,2}=cX(5);
    
    tbl_data{7,1}='<HTML> &Delta;X (px)</HTML>';
    tbl_data{7,2}=range(Xc(:,nn));
    tbl_data{8,1}='<HTML> Mean(x) </HTML>';
    tbl_data{8,2}=mean(Xc(:,nn));
    
    sTblX.tbl_data=tbl_data;
    sTblX.Position(3)=sTblX.Extent(3);
    sTblX.Position(4)=sTblX.Extent(4); 
    
    % Y Fit
    tbl_data={};
    
    
    axes(hax2);
    fit2=makeSineDecayFit(xvals',Yc(:,nn));
    cIntY=confint(fit2);
    
    %tbl_data{1,3}=range(cInt(:,1))/2;
    %tbl_data{2,3}=range(cInt(:,2))/2;

    %tbl_data{3,3}=1./(range(cInt(:,2))/2);
    %tbl_data{4,3}=range(cInt(:,3))/2;
    %tbl_data{5,3}=range(cInt(:,4))/2;
    %tbl_data{6,3}=range(cInt(:,5))/2;
    
    
    plot(tVec,feval(fit2,tVec),'r-');  

    cY=coeffvalues(fit2);
    
    sTblY.tbl_data={};
    tbl_data{1,1}='amp (px)';
    tbl_data{2,1}='period';
    tbl_data{3,1}='freq';
    tbl_data{4,1}='phase (rad)';
    tbl_data{5,1}='offset (px)';
    tbl_data{6,1}='tau ';

    tbl_data{1,2}=cY(1);
    tbl_data{2,2}=cY(2);
    tbl_data{3,2}=1/cY(2);

    tbl_data{4,2}=cY(3);
    tbl_data{5,2}=cY(4);
    tbl_data{6,2}=cY(5);
    
    tbl_data{7,1}='<HTML> &Delta;Y (px)</HTML>';
    tbl_data{7,2}=range(Yc(:,nn));
    tbl_data{8,1}='<HTML> Mean(y) </HTML>';
    tbl_data{8,2}=mean(Yc(:,nn));
    
    sTblY.tbl_data=tbl_data;
    sTblY.Position(3)=sTblY.Extent(3);
    sTblY.Position(4)=sTblY.Extent(4); 
    drawnow;
end


if opts.CenterParabolaFit && length(xvals)>1
    tVec=linspace(min(xvals),max(xvals),100);   
    
    tbl_dataX={};
    tbl_dataY={};

    for nn=1:size(Xc,2)
        D1=Xc(:,nn);    
        D2=Yc(:,nn);    
        % X Fit
        axes(hax1);
        fit1=polyfit(xvals',D1,2);
        plot(tVec,polyval(fit1,tVec),'r-','linewidth',1);  

        sTblX.Data={};

        tbl_dataX{1+6*(nn-1),1}='curvature (px/var^2)';
        tbl_dataX{2+6*(nn-1),1}='curvature (um/var^2)';
        tbl_dataX{3+6*(nn-1),1}='slope (px/var)';
         tbl_dataX{4+6*(nn-1),1}='slope (um/var)';
        tbl_dataX{5+6*(nn-1),1}='intercept (px)';
         tbl_dataX{6+6*(nn-1),1}='intercept (um) ';

        tbl_dataX{1+6*(nn-1),2}=fit1(1);
        tbl_dataX{2+6*(nn-1),2}=fit1(1)*PixelSize*1e6;
        tbl_dataX{3+6*(nn-1),2}=fit1(2);
        tbl_dataX{4+6*(nn-1),2}=fit1(2)*PixelSize*1e6;
        tbl_dataX{5+6*(nn-1),2}=fit1(3);
        tbl_dataX{6+6*(nn-1),2}=fit1(3)*PixelSize*1e6;    

    %     tbl_data{5,1}='<HTML> &Delta;X (px)</HTML>';
    %     tbl_data{5,2}=range(Xc(:,nn));
    %     tbl_data{6,1}='<HTML> Mean(x) </HTML>';
    %     tbl_data{6,2}=mean(Xc(:,nn));

        sTblX.Data=tbl_dataX;
        sTblX.Position(3)=sTblX.Extent(3);
        sTblX.Position(4)=sTblX.Extent(4); 

        % X Fit
        axes(hax2);
        fit2=polyfit(xvals',D2,2);
        
        if median(D2)>1024
            yyaxis right
        else
            yyaxis left
        end
        plot(tVec,polyval(fit2,tVec),'r-','linewidth',1);  

        sTblY.Data={};
        tbl_dataY{1+6*(nn-1),1}='curve (px/var^2)';
        tbl_dataY{2+6*(nn-1),1}='curve (um/var^2)';
        tbl_dataY{3+6*(nn-1),1}='slope (px/var)';
        tbl_dataY{4+6*(nn-1),1}='slope (um/var)';
        tbl_dataY{5+6*(nn-1),1}='intercept (px)';
        tbl_dataY{6+6*(nn-1),1}='intercept (um) ';

        tbl_dataY{1+6*(nn-1),2}=fit2(1);
        tbl_dataY{2+6*(nn-1),2}=fit2(1)*PixelSize*1e6;
        tbl_dataY{3+6*(nn-1),2}=fit2(2);
        tbl_dataY{4+6*(nn-1),2}=fit2(2)*PixelSize*1e6;
        tbl_dataY{5+6*(nn-1),2}=fit2(3);
        tbl_dataY{6+6*(nn-1),2}=fit2(3)*PixelSize*1e6;   

        tbl_dataY{5+6*(nn-1),1}='<HTML> &Delta;Y (px)</HTML>';
        r = max(Yc(:,nn)) - min(Yc(:,nn));
        tbl_dataY{5+6*(nn-1),2}=r;
        tbl_dataY{6+6*(nn-1),1}='<HTML> Mean(y) </HTML>';
        tbl_dataY{6+6*(nn-1),2}=mean(Yc(:,nn));

        sTblY.Data=tbl_dataY;
        sTblY.Position(3)=sTblY.Extent(3);
        sTblY.Position(4)=sTblY.Extent(4); 
    end
end



if opts.CenterLinearFit && length(xvals)>1
    tVec=linspace(min(xvals),max(xvals),100);   
    
    D1=Xc(:,nn);    
    D2=Yc(:,nn);

    
    % X Fit
    axes(hax1);
    fit1=polyfit(xvals',D1,1);
    plot(tVec,polyval(fit1,tVec),'r-','linewidth',1);  
    
    sTblX.Data={};
    tbl_data{1,1}='slope (px/var)';
    tbl_data{2,1}='slope (um/var)';
    tbl_data{3,1}='intercept (px)';
    tbl_data{4,1}='intercept (um) ';

    tbl_data{1,2}=fit1(1);
    tbl_data{2,2}=fit1(1)*PixelSize*1e6;
    tbl_data{3,2}=fit1(2);
    tbl_data{4,2}=fit1(2)*PixelSize*1e6;
    
    tbl_data{5,1}='<HTML> &Delta;X (px)</HTML>';
    r = max(Xc(:,nn)) - min(Xc(:,nn));
    tbl_data{5,2}=r;
    tbl_data{6,1}='<HTML> Mean(x) </HTML>';
    tbl_data{6,2}=mean(Xc(:,nn));
    
    sTblX.Data=tbl_data;
    sTblX.Position(3)=sTblX.Extent(3);
    sTblX.Position(4)=sTblX.Extent(4); 
    
     % Y Fit
    axes(hax2);
    if median(D2)>1024
        yyaxis right
    else
        yyaxis left
    end
    fit2=polyfit(xvals',D2,1);
    plot(tVec,polyval(fit2,tVec),'r-','linewidth',1);  
    
    sTblY.Data={};
    tbl_data{1,1}='slope (px/var)';
    tbl_data{2,1}='slope (um/var)';
    tbl_data{3,1}='intercept (px)';
    tbl_data{4,1}='intercept (um) ';

    tbl_data{1,2}=fit2(1);
    tbl_data{2,2}=fit2(1)*PixelSize*1e6;
    tbl_data{3,2}=fit2(2);
    tbl_data{4,2}=fit2(2)*PixelSize*1e6;
    
    tbl_data{5,1}='<HTML> &Delta;Y (px)</HTML>';
    r = max(Yc(:,nn)) - min(Yc(:,nn));
    tbl_data{5,2}=r;
    tbl_data{6,1}='<HTML> Mean(y) </HTML>';
    tbl_data{6,2}=mean(Yc(:,nn));
    
    sTblY.Data=tbl_data;
    sTblY.Position(3)=sTblY.Extent(3);
    sTblY.Position(4)=sTblY.Extent(4); 
    
    

end

% hax1.Position(4)=hax1.Position(4)-15;
% hax2.Position(4)=hax1.Position(4);

%% Make Figure

hF2=figure('Name',[pad([data.FitType ' centre 2'],20) FigLabel],...
    'units','pixels','color','w','numbertitle','off');
hF2.Position(1)=1020;
hF2.Position(2)=380;
hF2.Position(3)=350;
hF2.Position(4)=350;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF2.Position(3);
t.Position(1:2)=[5 hF2.Position(4)-t.Position(4)];

uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

hax1=axes;
set(hax1,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on
xlabel(['x centre'],'interpreter','none');
ylabel(['y centre (circle)'],'interpreter','none');

% yyaxis right

co=get(gca,'colororder');

yyaxis left
set(gca,'YColor','k');
yyaxis right
set(gca,'YColor','k');
ylabel(['y centre (square)'],'interpreter','none');

for nn=1:size(Xc,2)    
    if median(Yc(:,nn))>1092
        yyaxis right
        m = 's';        
        yL = min([yL min(Yc(:,nn))-1024]);
        yH = max([yH max(Yc(:,nn))-1024]);
    else
        yyaxis left        
        m = 'o';
        yL = min([yL min(Yc(:,nn))]);
        yH = max([yH max(Yc(:,nn))]);
    end
    
    
   plot(Xc(:,nn),Yc(:,nn),m,'color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);   
end

yLCam = [yL yH]+(yH-yL)*.1*[-1 1];

yyaxis left
ylim(yLCam);

yyaxis right
ylim(yLCam+1024);

% if opts.angleTrack
%     fit_xy=polyfit(Xc(:,1),Yc(:,1),1);
%     xx=linspace(min(Xc(:,1)),max(Xc(:,1)),100);
% 
%     theta = atan(1/fit_xy(1))*180/pi;
% 
%     pf=plot(xx,polyval(fit_xy,xx),'r-','linewidth',1);  
%     str = ['$\theta = ' num2str(round(theta,2)) '^\circ$'];
%     legend(pf,{str},'location','best','interpreter','latex');
%     
% %     text(5,5,str,'units','pixels','interpreter','latex','verticalalignment','bottom',...
% %         'backgroundcolor',);
% end

resizeFig(hF2,t)
end

function fitResult=makeSineDecayFit(X,Y,W)

% Guess the amplitude and offset\
r = max(Y) - min(Y);

gA=0.5*r;
gD=(max(Y)+min(Y))*.5;

% Guess the period
iHigh=find((Y-gD)/gA>.8,1);
iLow=find((Y-gD)/gA<-.8,1);
gB=abs(X(iHigh)-X(iLow))*2.2;



minValues=X(Y==min(Y));
maxValues=X(Y==max(Y));
% gB=1*abs(maxValues(1)-minValues(1));
% gB=range(X)/2;




gC=maxValues(1);
gC=pi;
gD=0.5*(max(Y)+min(Y));

gC=pi;
gE = range(X);

cosFit=fittype('A*cos(2*pi*t/B+C)*exp(-t/E)+D','independent',{'t'},...
    'coefficients',{'A','B','C','D','E'});
options=fitoptions(cosFit);          
        set(options, 'TolFun', 1E-14);
        set(options,'Lower', [0.25*gA,...
            .1*gB,...
            0, ...
            0.75*gD, ...
            0]);
        set(options, 'Upper', [5*gA, ...
            20*gB,...
            2*pi, ...
            1.5*gD, ...
            inf]);            
        set(options, 'StartPoint', [gA, gB,...
            gC,gD, gE]);     
        set(options, 'MaxIter',3000);
        set(options, 'MaxFunEvals',3000);
        set(options,'TolFun',10^-9);
        
        if nargin==3
           set(options,'Weights',W); 
        end        
        
        
        fitResult=fit(X,Y,cosFit,options);      
        
disp(fitResult)
end

