function matOut = shiftMatrix(mat,dx,dy)
%SHIFTMATRIX returns the same matrix, but shifted by dx and dy pixels in the
%two coordinate directions. It fills in zeros when needed at the edges.


[nr, nc]=size(mat);

matOut=mat;
if dx>=0
    matOut=[zeros(nr,dx) matOut];
    matOut=matOut(:,1:nc);
else
    matOut=matOut(:,(1-dx):end);
    matOut=[matOut zeros(nr,-1*dx)];
end

if dy>=0
    matOut=[zeros(dy,nc); matOut];
    matOut=matOut(1:nr,:);
else
    matOut=matOut((1-dy):end,:);
    matOut=[matOut; zeros(-1*dy,nc)];
end


end

