function atomdata = matchParamsFlags(atomdata)

pAll = sort(fieldnames(atomdata(1).Params));
fAll = sort(fieldnames(atomdata(1).Flags));

badParams = 0;
badFlags  = 0;

% Find Master list of all params and flags
for kk=1:length(atomdata)
    pThis = sort(fieldnames(atomdata(kk).Params));      
    fThis = sort(fieldnames(atomdata(kk).Flags));
    
    if ~isequal(pThis,pAll)
        badParams = 1;
        for nn=1:length(pThis)      
            ind=findStr(pAll,pThis{nn});             
            if isempty(ind)
                pAll{end+1}=pThis{nn}; 
            end  
       end
    end 
    
    if ~isequal(fThis,fAll)
        badFlags = 1;
        for nn=1:length(fThis)      
            ind=findStr(fAll,fThis{nn});             
            if isempty(ind)
                fAll{end+1}=fThis{nn}; 
            end  
       end
    end 

end

if badParams
   warning('Unequal number of parameters detected. Adding NaN'); 
end

if badFlags
   warning('Unequal number of flags detected. Adding NaN'); 
end

% Add blank values from the master list
for kk=1:length(atomdata)
    pThis = sort(fieldnames(atomdata(kk).Params));      
    fThis = sort(fieldnames(atomdata(kk).Flags));
    
    if ~isequal(pAll,pThis)        
       for nn=1:length(pAll)      
            ind=findStr(pThis,pAll{nn});   
            if isempty(ind)   
                atomdata(kk).Params.(pAll{nn}) = NaN ;
                atomdata(kk).Units.(pAll{nn}) = NaN ;
            end       
       end
    end
    
    if ~isequal(fAll,fThis)        
       for nn=1:length(fAll)      
            ind=findStr(fThis,fAll{nn});   
            if isempty(ind)
               atomdata(kk).Flags.(fAll{nn}) = NaN ;
            end       
       end 
    end 
    
end

if badParams
    for kk=1:length(atomdata)
        atomdata(kk).Params = orderfields(atomdata(kk).Params);
        atomdata(kk).Units = orderfields(atomdata(kk).Units);

    end
end

if badFlags
    for kk=1:length(atomdata)
        atomdata(kk).Flags = orderfields(atomdata(kk).Flags);
    end
end

end

function ind=findStr(cellstr,str)

ind=[];
for kk=1:length(cellstr)
   if isequal(str,cellstr{kk})
       ind = kk;
   end
end

end

