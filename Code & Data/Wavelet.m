function [A, D]=DWT(x,nLevel,wname)

[C, L]=wavedec(x,nLevel,wname);

A=appcoef(C,L,wname,nLevel);
A=cell(nLevel,1);
D=cell(nLevel,1);
for i=1:nLevel
    A{i}=detcoef(C,L,i);
    D{i}=detcoef(C,L,i);
end
end