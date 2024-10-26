% given a full dimensional simplex conv(A) in R^n, where A is an
% n-by-(n+1) matrix, this function returns the John ellipsoid E(F,c)=FB+c of the
% simplex; here B is the unit ball of R^n
function [F,c]= JE_simplex(A)
[n,p]= size(A); %p=n+1
[~,~,B]= PCA_ASF(eye(p),p); %JE of conv(B) is the ball B(r), where r is computed below 
u= ones(p,1)/p; v= ones(p,1)/n; v(end)=0; r= norm(u-v,2);
B=B/r; %JE of conv(B) is now the unit ball B(1)=E(G,d) with G=I_n and d=0

c= mean(A,2);
H= B(:,1:n);
J= A(:,1:n)-c;
F= J*inv(H); 
F= sqrtm(F*F');
return;