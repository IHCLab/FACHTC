function [C,d,Yd]= PCA_ASF(Y,p)
[l,n]= size(Y);
d= mean(Y,2);
U= Y-d*ones(1,n);
R= U*U';
[eV,~]= eig(R);
C= eV(:,l-p+2:end);
Yd= C'*(Y-d*ones(1,n));
