function saveFigure(hF,filename,opts)

if nargin==2
    opts=struct;
    opts.Quality = '-r120';
    opts.Format = '-dpng';
    opts.saveDir = pwd;
    ext='.png';
end

if ~isfield(opts,'Quality') || isequal(opts.Quality,'auto')
    save_qual='-r120';
else
    save_qual=opts.Quality;
end

if ~isfield(opts,'Format') || isequal(opts.Format,'auto')
    imgformat='-dpng';
else
    imgformat=opts.Format;
end

switch imgformat
    case '-dpng'
        ext = '.png';
    case '-dpcx256'
        ext = '.pcx';
    case '-djpeg'
        ext = '.jpg';
    case '-dbmpmono'
        ext = '.bmp';
end

figDir = opts.saveDir;



for kk=1:length(hF)
    if length(hF)>1
       filename = [filename num2str(kk)]; 
    end
    
    % Make the figure name with the location
    saveLocation=fullfile(figDir,[filename ext]);
    % saveLocation='C:

    % Save the figure and the png
    fprintf([datestr(now,13) ' Saving figure handle to ']);
    fprintf([filename ext ' ... ']);
    
    figure(hF(kk))
    set(0,'CurrentFigure', hF(kk));
    set(hF(kk),'PaperPositionMode','auto');
    print(imgformat,save_qual,saveLocation);
    disp('Saved!');
end

% savefig(the_figure,saveLocation);
end

