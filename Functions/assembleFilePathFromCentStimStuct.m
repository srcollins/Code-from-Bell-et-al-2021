function fileNamePath=assembleFilePathFromCentStimStuct(root,res,resInd)
%% This function is designed to read the data stored in the Center 
%stimulation experiment structure array (res), and assemble a filepath for the
%the image files stored in the imData structure. The imData files for one
%stimulation type are stored in the analyzed data folder for the
%corresponding experimental folder. All of the files are not in the same
%folder.

parentFold=res(resInd).parentFolder;
fName=cellfun(@(x) sprintf('%s-imgData.mat',x),res(resInd).fileName,'UniformOutput',false);

fileNamePath=cellfun(@(x,y) [root 'Analyzed Data' filesep x filesep y],parentFold,fName,'UniformOutput',false);



end