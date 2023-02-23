function [data_out] = chDataXVar(data_in,xVar)
% Given a data structure such as erf_data, gauss_data, or box_data, this
% function changes the X variable to match the x variable provided.

goodVar = 1;
for kk=1:length(data_in.Params)
    if ~isfield(data_in.Params(kk),xVar)
       goodVar = 0; 
    end
end

if ~goodVar
   error('not all data has the variable requested. Aborting'); 
end

X = zeros(length(data_in.X),1);
for kk = 1 :length(data_in.Params)
   X(kk) = data_in.Params(kk).(xVar); 
end

data_out = data_in;
data_out.X = X;

end

