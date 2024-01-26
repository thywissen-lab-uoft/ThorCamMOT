function [dirDay] = getImageDir(mydatevec)

dirSource = 'C:\Data';
dirYear  = [dirSource filesep num2str(mydatevec(1))];
dirMonth = [dirYear filesep num2str(mydatevec(1)) '.' sprintf('%2.2d',mydatevec(2))];
dirDay   = [dirMonth filesep sprintf('%2.2d',mydatevec(2)) '.' sprintf('%2.2d',mydatevec(3))];

if ~exist(dirSource)
   warning('No data server found.');
   dirDay = 0;
   return;  
end

if exist(dirYear)~=7
   mkdir(dirYear); 
end

if exist(dirMonth)~=7
    mkdir(dirMonth);
end

if exist(dirDay)~=7
   mkdir(dirDay); 
end

end



