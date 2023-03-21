
function [fout,gof,output]=gaussFit2D(Dx,Dy,data)

% mask = isnan(data);

% data(mask)=0;

 % Ensure data type is double
data=double(data);Dx=double(Dx);Dy=double(Dy);

% Rescale images for fitting speed (Do this adaptively? or on option?)
sc=0.4; % Scale factor
data=imresize(data,sc);Dx=imresize(Dx,sc);Dy=imresize(Dy,sc);

% mask = imresize(mask,sc);

dSmooth=imgaussfilt(data,2);    % Smooth data
N0=max(max(dSmooth));           % Extract peak, amplitude guess

% Remove low data points
Z=dSmooth;Z(dSmooth<N0*.3)=0;

% Calculate guesses for center and size
X=sum(Z,1);Y=sum(Z,2)';             % Get X and Y sum profiles
Nx=sum(X);Ny=sum(Y);                % Get the total number of counts
Xc=mean(Dx(X>.9*max(X)));           % X center (use >90% SNR)
Yc=mean(Dy(Y>.9*max(Y)));           % Y center (use >90% SNR)
Xs=1.5*sqrt(sum((Dx-Xc).^2.*X)/Nx); % X standard deviation * 1.5
Ys=1.5*sqrt(sum((Dy-Yc).^2.*Y)/Ny); % Y standard deviation * 1.5

% Make a mesh grid for fitting
[xx,yy]=meshgrid(Dx,Dy);

% Make an initial guess
Zguess=N0*exp(-(xx-Xc).^2./(2*Xs)^2).*exp(-(yy-Yc).^2./(2*Ys)^2);

% Copy the data
data2=data;xx2=xx;yy2=yy;

% Elminate data points below a threshold to reduce # points to fit
th=0.1;
% th=-.1;
xx2(Zguess<th*N0)=[];yy2(Zguess<th*N0)=[];data2(Zguess<th*N0)=[];

 %CHEATING
inds = [data2 == 0];
xx2(inds)=[];
yy2(inds)=[];
data2(inds)=[];

% xx2(mask)=[];yy2(mask)=[];data2(mask)=[];

% Calculate the appropriate background
bg=sum(sum(data-Zguess))/(length(X)*length(Y));

bg = min(min(data));

% Create fit object
myfit=fittype('A*exp(-(xx-Xc).^2./(2*Xs^2)).*exp(-(yy-Yc).^2./(2*Ys^2))+nbg',...
    'independent',{'xx','yy'},'coefficients',{'A','Xc','Xs','Yc','Ys','nbg'});
opt=fitoptions(myfit);
opt.StartPoint=[N0 Xc Xs Yc Ys bg];
opt.Lower=[N0/10 10 1 10 1 -.1];
% opt.Upper=[1.5*N0 1.5*max(Dx) range(Dx) 1.5*max(Dy) range(Dy) 0.1];
opt.Upper=[1.5*N0 1.5*max(Dx) range(Dx) 1.5*max(Dy) range(Dy) inf];

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

% Display initial guess
str1=['(Xc0,Yc0)=(' num2str(round(Xc)) ',' num2str(round(Yc)) ');'];
str2=['(Xs0,Ys0)=(' num2str(round(Xs)) ',' num2str(round(Ys)) ')'];
fprintf([str1 str2 ';']);

% Perform the fit
fprintf(' gauss fitting...');
t1=now;
[fout,gof,output]=fit([xx2(:) yy2(:)],data2(:),myfit,opt);
t2=now;
disp([' done (' num2str(round((t2-t1)*24*60*60,1)) ' sec.).']);

% Perform the fit again if it's bad, assume zero atoms for starting point
if gof.rsquare<0.5
    opt.StartPoint=[0 mean(Dx) 10 mean(Dy) 10 mean(data2(:))];  
    opt.Upper=[std(data2(:)) max(Dx) .3*range(Dy) max(Dy) .3*range(Dy) mean(data2(:))+std(data2(:))];  
    opt.Lower=[0 min(Dx) 0 min(Dy) 0 mean(data2(:))-std(data2(:))];  

    fprintf(' fitting...');
    t1=now;
    [fout,gof,output]=fit([xx2(:) yy2(:)],data2(:),myfit,opt);
    t2=now;
    disp([' done (' num2str(round((t2-t1)*24*60*60,1)) ' sec.).']);
end

end