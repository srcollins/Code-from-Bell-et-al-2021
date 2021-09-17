function printFigurePDF(outname,printOptions)
% Saves the current figure as an EPS file
if ~isempty(printOptions)
    
    if ~boolRegExp({outname},'\.pdf$')
        outname=[outname '.pdf'];
    end
    
    print('-dpdf', '-r300',printOptions, outname);
    
    
else
    if ~boolRegExp({outname},'\.pdf$')
        outname=[outname '.pdf'];
    end
    
    print('-dpdf', '-r300', outname);
end