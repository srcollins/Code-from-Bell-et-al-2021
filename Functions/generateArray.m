function arrayOut=generateArray(CtrPt,NrPts,Res)
% 2021-03-27: This function will build an evenly spaced array around a
% given centerpoint. I use this function for cropping images for movie or
% montage generation. I borrowed this from the mathworks documentation page. 

%Inputs: CtrPt=center point of the array. NrPts=total lenght of the array.
%Res=step size between points. 

%%
arrayOut=linspace(CtrPt-Res*fix(NrPts/2), CtrPt+Res*fix(NrPts/2), NrPts);
end