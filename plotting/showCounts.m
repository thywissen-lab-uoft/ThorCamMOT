function hF = showCounts(data,xVar,opts)

if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end

%% Sort the data by the parameter given
params=[data.Params];
X=[params.(xVar)]';
Ncounts = data.Ncounts;

%% Exponential Decay Fit

if isfield(opts,'NumberExpFit') && opts.NumberExpFit && size(Ncounts,1)>2
    myfit=fittype('A*exp(-t/tau)','coefficients',{'A','tau'},...
    'independent','t');
    opt=fitoptions(myfit);
    
    % Get some initial guesses
    tau0=max(X)/2;
%     tau0=0.1;
    
    fout_exp={};
    for nn=1:size(Ncounts,2)  
        A0=max(Ncounts(:,nn));
        
        % Assign start point
        opt.StartPoint=[A0 tau0];
        fout_exp{nn}=fit(X,Ncounts(:,nn),myfit,opt);
    end
end

if isfield(opts,'NumberExpOffsetFit') && opts.NumberExpOffsetFit && size(Ncounts,1)>3
    myfit=fittype('A*exp(-t/tau)+B','coefficients',{'A','tau','B'},...
    'independent','t');
    opt=fitoptions(myfit);
    
    fout_exp={};
    
    for nn=1:size(Ncounts,2)
       A0=range(Ncounts(:,nn));
       B0=min(Ncounts(:,nn));
       tau0=range(X)/4;
       
       opt.StartPoint=[A0 tau0 B0];
       opt.Lower=[0 0 0];
       fout_exp{nn}=fit(X,Ncounts(:,nn),myfit,opt);       
    end       
end

%% Lorentzian Fits
% Addding lorentzian fits, find a better way to put this in the code in
% the future, future me. Thanks <3

doLorentzianFit=opts.NumberLorentzianFit;
fouts_lorentz={};
if doLorentzianFit
    for rr=1:size(Ncounts,2)
        myfit=fittype('A*((x-x0).^2+(G/2).^2).^(-1).*(G/2)^2','coefficients',{'A','G','x0'},...
            'independent','x');
        opt=fitoptions(myfit);
        A0=max(Ncounts(:,rr));
        G0=range(X)/2;

        inds=[Ncounts(:,rr)>.8*max(Ncounts(:,rr))];
        x0=mean(X(inds));       
        opt.StartPoint=[A0 G0 x0];    
        fout_lorentz=fit(X,Ncounts(:,rr),myfit,opt);
        fouts_lorentz{rr}=fout_lorentz;
    end
end


%% Make Figure

hF=figure('Name',[pad([data.FitType ' number'],20) FigLabel],...
    'units','pixels','color','w','Menubar','figure','Resize','on',...
    'numbertitle','off');
hF.Position=[5 380 500 300];clf;

% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
drawnow;
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

% Draw PCO label
uicontrol('style','text','string','MOT THORCAM','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

% Make axis
hax=axes;co=get(gca,'colororder');
set(hax,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel([data.FitType ' counts']);
        
% Plot the data
for nn=1:size(Ncounts,2)
   plot(X,Ncounts(:,nn),'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end

if (opts.NumberExpFit || opts.NumberExpOffsetFit) && size(Ncounts,1)>3
    strs={};
    xx=linspace(0,max(X),1000);
    
    for nn=1:size(Ncounts,2)
        pExp(nn)=plot(xx,feval(fout_exp{nn},xx),'-','linewidth',1,...
            'color',0.8*co(nn,:)); 
        
        if opts.NumberExpFit        
            fstr=['$N_0 = ' num2str(round(fout_exp{nn}.A),'%.2e') '$' newline ...
                '$\tau = ' num2str(round(fout_exp{nn}.tau,1),'%.2e') ' $'];
        else
            fstr=['$\Delta N = ' num2str(round(fout_exp{nn}.A),'%.2e') '$' newline ...
                '$\tau = ' num2str(round(fout_exp{nn}.tau,1),'%.2e') ' $' newline ...
                '$N_0 = ' num2str(round(fout_exp{nn}.B),'%.2e') '$'];
        end
        strs{nn}=fstr;
    end       
    legend(pExp,strs,'interpreter','latex','location','best','fontsize',8);
    hax.YLim(1)=0;
end

if doLorentzianFit
    xx=linspace(min(X),max(X),100);
    legStr={};
    
    for rr=1:length(fouts_lorentz)
        fout_lorentz=fouts_lorentz{rr};
        pFs(rr)=plot(xx,feval(fout_lorentz,xx),'-','linewidth',2,'color',...
        co(rr,:)*.8);

        str=['$N_0 = ' num2str(round(fout_lorentz.A),'%.2e') '$' newline ...
            '$\mathrm{FWHM} = ' num2str(round(fout_lorentz.G,3)) ' $' newline ...
            '$x_0 = ' num2str(round(fout_lorentz.x0,3)) '$'];
        legStr{rr}=str;
    end
    legend(pFs,legStr,'interpreter','latex','location','best','fontsize',8);

    hax.YLim(1)=0;
end

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

resizeFig(hF,t,[hax]);
end

