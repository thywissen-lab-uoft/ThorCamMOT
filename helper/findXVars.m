function [xVars] = findXVars(atomdata)
% Find all changing variables in the atomdata that are single valued
% numbers.

params = [atomdata.Params];
fnames = fieldnames(params);

xVars={};

for kk=1:length(fnames)
    fname = fnames{kk};
    
    val_0 = params(1).(fname);
    
    % Ignore strings and arrays
    if ~isnumeric(val_0) || length(val_0)~=1
        continue
    end
    
    % Ignore timestamp
    if isequal(fname,'timestamp') 
        continue;
    end
    
    vals = [params.(fnames{kk})];    
    u_vals = unique(vals);

    if length(u_vals)>1
        xVars{end+1}=fnames{kk};
    end
end

% Find the index where it is execution date
ind = find(strcmp(xVars,'ExecutionDate')==1,1);

% Remove Execution date
xVars(ind)=[];

% Move it to the end
xVars{end+1}='ExecutionDate';

end

