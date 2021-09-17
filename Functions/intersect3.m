function [Com,ia,ib,ic] = intersect3(A,B,C)
[C1,ia,ib] = intersect(A,B);
[Com,ic1,ic] = intersect(C1,C);
%~ ic is okay
ia = ia(ic1);
ib = ib(ic1);
end