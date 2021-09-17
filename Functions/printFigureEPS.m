function printFigureEPS(outname)
% Saves the current figure as an EPS file

if ~boolRegExp({outname},'\.eps$')
    outname=[outname '.eps'];
end

print('-depsc2', '-r300', outname);
