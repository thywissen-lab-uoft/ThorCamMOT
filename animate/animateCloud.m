function animateCloud(atomdata,xVar,opts)


if isfield(opts,'saveDir')
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
end

%% Animate Settings
startDelay=opts.StartDelay;
midDelay=opts.MidDelay;
endDelay=opts.EndDelay;

%% Make Filename

filename='animate'; 


if ~exist(opts.saveDir,'dir')
   mkdir(opts.saveDir); 
end

% Make the figure name with the location
filename=fullfile(opts.saveDir,[filename '.gif']);

%% Sort the Data

params = [atomdata.Params];
xvals = [params.(xVar)];
direction = opts.Order;

if isequal(direction,'ascend')
    [~,inds]=sort(xvals,'ascend');    
else
    [~,inds]=sort(xvals,'descend');    

end

atomdata=atomdata(inds);
params =[atomdata.Params];
xvals =[params.(xVar)];

ROI = atomdata(1).ROI;

%% Grab all data


% Initiate X,Y, and Z data over the display ROI
X = ROI(1,1):ROI(1,2);
Y = ROI(1,3):ROI(1,4);
Z = zeros(length(Y),length(X),length(atomdata));


% Grab all optical densities in the display ROI
for kk=1:length(atomdata)
    Z(:,:,kk) = atomdata(kk).Data(Y,X);   
end

if opts.doAverage
    xvals=unique(xvals);    
    Zu = zeros(length(Y),length(X),length(xvals));    
    for kk=1:length(xvals)        
        inds = xvals(kk) == [params.(xVar)];        
        Zu(:,:,kk)=mean(Z(:,:,inds),3);     
    end
else
    Zu = Z;
end

% %% Grab initialize data
% X=atomdata(1).X;
% Y=atomdata(1).Y;
% Z=atomdata(1).Data;
% [xx,yy]=meshgrid(X,Y);
% 
% % Find max
% N=[0 1];
% 
% for kk=1:length(atomdata)
%     N(2)=max([max(max(atomdata(kk).Data)) N(2)]);
% %     N(1)=min([min(min(atomdata(kk).Data)) N(1)]);
% end
% 
% % For the animation, change this prefactor to change the color limits to
% % see signals visually.  Here "N" is the maxmium count on a single pixel,
% % but hot pixels can sometimes mnake this too high, so that's lwhy I leave
% % the prefactor in front.  Feel free to change this at will in order to
% % make your gifs look sensible.
% N=[0 500];


%% Make Figure

H=ROI(4)-ROI(3);
W=ROI(2)-ROI(1);

L = 700;

if W>H
    H = (H/W)*L;
    W = L;
else
    W = (W/H)*L;
    H = L;
end


lgap = 50;
rgap = 50;

bgap = 20;
tgap = 40;

w = W + lgap + rgap;
h = H + bgap + tgap;

hF=figure('Name',[FigLabel ' : Animate Cloud'],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'WindowStyle','modal');
hF.Position(1)=100;
hF.Position(2)=100;
hF.Position(3)=w;
hF.Position(4)=h;
drawnow;

t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

hAxImg=axes('parent',hF,'units','pixels','Box','on','XGrid','on',...
    'YGrid','on','YDir','reverse','XAxisLocation','bottom');
hAxImg.Position(2)=bgap;
hAxImg.Position(1)=lgap;
hAxImg.Position(3:4)=[W H];
drawnow;


t=text(5,5,'hi','units','pixels','fontsize',16,'color','w',...
    'interpreter','none','verticalalignment','bottom');

colormap inferno;
hold on
hImg=imagesc(X,Y,Z(:,:,1));
axis equal tight
caxis(opts.CLim);
hold on
colorbar


%% Iterate   

for kk=1:length(xvals)
    t.String=[xVar ' = ' num2str(xvals(kk)) ' (' opts.xUnit ')'];          % Variable string   
    set(hImg,'XData',X,'YData',Y,'CData',Zu(:,:,kk));          

    drawnow
    frame = getframe(hF);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);           

    if kk == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',startDelay);
    else
        if kk==length(xvals)
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',endDelay);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',midDelay);
        end
    end
        
end
close;
    
end

