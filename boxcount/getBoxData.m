function output = getBoxData(atomdata,xVar)
%GETERFDATA Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
   xVar = 'ExecutionDate'; 
end

%% Sort the data by the parameter given
params=[atomdata.Params];
X=[params.(xVar)];

[X,inds]=sort(X,'ascend');
atomdata=atomdata(inds);

% Make sure its Nx1
X = reshape(X,[length(X) 1]);

%% Grab the Erf Fit outputs


for kk=1:length(atomdata)
   for nn=1:length(atomdata(kk).BoxCount)
        bc=atomdata(kk).BoxCount(nn);               % Grab the fit
        

        N(kk,nn)=bc.Ncounts;    % Number of counts (w/ bkgd removed)
        Nraw(kk,nn)=bc.Nraw;          % Raw of number of counts
        Nbg(kk,nn)=bc.Nbkgd;          % Bakcground number of counts
        nbg(kk,nn)=bc.nbkgd;          % Background counts/px
        Xc(kk,nn)=bc.Xc;              % X center of mass
        Yc(kk,nn)=bc.Yc;              % Y center of mass
        Xs(kk,nn)=bc.Xs;              % X standard deviation
        Ys(kk,nn)=bc.Ys;              % Y standard deviation                      
        Zs(kk,nn)=bc.Xs;              % Y standard deviation                      

        Natoms(kk,nn)=N(kk,nn);  % Atom number      
   end        
end

output = struct;

output.FileNames    = {atomdata.Name}';

% if isfield(atomdata(1),'Flags')
%     output.Atom         = atomdata(1).Flags.image_atomtype;
% else
%     output.Atom         = NaN;
% end

output.xVar         = xVar;
output.X            = X;
output.Params       = params;
output.Units        = [atomdata.Units];
output.Flags        = [atomdata.Flags];
output.FitType      = 'box';

% Assign fit outputs
output.Ncounts      = Natoms;
output.Xc           = Xc;
output.Yc           = Yc;
output.Xs           = Xs;
output.Ys           = Ys;
output.Zs           = Zs;
output.nbg          = nbg;

end

