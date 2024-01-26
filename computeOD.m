function [atomdata] = computeOD(atomdata)

for kk=1:length(atomdata)
    pwa = double(atomdata(kk).Images(:,:,1));
    pwoa = double(atomdata(kk).Images(:,:,2));
    od = log(pwoa./pwa);
    atomdata(kk).Data = od;

end

end

