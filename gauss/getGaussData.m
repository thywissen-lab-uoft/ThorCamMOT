function output = getGaussData(atomdata,xVar)
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

PixelSize = atomdata(1).PixelSize;

for kk=1:length(atomdata)
   for nn=1:length(atomdata(kk).GaussFit)
        fout=atomdata(kk).GaussFit{nn};               % Grab the fit
        fits{kk,nn}=fout;
        GOFs{kk,nn}=atomdata(kk).GaussGOF{nn};
        R2s(kk,nn) = GOFs{kk,nn}.rsquare;
        SSEs(kk,nn) = GOFs{kk,nn}.sse;
        
        Xc(kk,nn)=fout.Xc;Yc(kk,nn)=fout.Yc;        % X and Y center
        Xs(kk,nn)=fout.Xs;Ys(kk,nn)=fout.Ys;        % X and Y sigma
        Zs(kk,nn)=fout.Xs;                          % ASSUME wZ=wX;    
        A(kk,nn)=fout.A;                            % Amplitude
        nbg(kk,nn)=fout.nbg;                        % Background

        N(kk,nn)=2*pi*Xs(kk,nn)*Ys(kk,nn)*A(kk,nn);     % Number of counts
   end        
end

output = struct;

output.FileNames    = {atomdata.Name}';

output.PixelSize    = PixelSize;
output.xVar         = xVar;
output.X            = X;
output.Params       = params;
output.Units        = [atomdata.Units];
output.Flags        = [atomdata.Flags];
output.FitType      = 'gauss';
output.Fits         = fits;
output.FitGOFs      = GOFs;
output.FitR2        = R2s;
output.FitSSE       = SSEs;

% Assign fit outputs
output.Ncounts       = N;
output.Xc           = Xc;
output.Yc           = Yc;
output.Xs           = Xs;
output.Ys           = Ys;
output.Zs           = Zs;
output.A            = A;
output.nbg          = nbg;

end

