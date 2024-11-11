function [C_est d_est index] = RASF(X,N,per)
%---- function of Robust Affine Set Fitting------
%---- inputs : X is the data input, N is the number of endmember, per
%             is the percentage(%) of outlier estimate.

%---- outputs: (C_est,d_est) is the estimated affine set parameter (C,d)
%             index is the estimated outlier index.
% programmers: 
%   Hao-En Huang (email: b23004705@hotmail.com) 
%   Tsung-Han Chan (email: thchan@ieee.org)
% Date: May 5, 2012 
%-------------------------------------------------------


%---Given the convergence tolerance, hyperspectral data, number of
%nedmember, and number of outlier-----%

[M L]     = size( X );   
o_percent = L*per/100;
tol       = 1 ;



%---- step 1 : initialize all the outlier estimate vector to "0" ----%
z         = zeros(M,L);
zz       = ones(M,L);
%--------------------------------------------------------------------%



%--- enter the while loop for testing

 while tol >1.0e-005
    
    
%---- step 2 : fix z_n, update d, C, and X_n

X_n = X - z;    
d = mean(X_n,2);                      %--- d is the mean data cloud 
U = X_n-d*ones(1,L);
OPTS.disp = 0;
[C D] = eigs(U*U',N-1,'LM',OPTS);    %--- C is the set which we desired
alpha = C'*(X_n-repmat(d,[1,L]));    %--- alpha is the optimal solution


%---  step 3 : fix Xn and update "z"  
       
    z = X - ( C*alpha + repmat(d,[1,L]) );
    e_norm = sqrt(sum(z.^2));
    
%--- find the L*per%  largest error in "e_norm" and its corresponding index  
   
   [eorror_norm error_index]=sort(e_norm,'descend');
   error_index=error_index(1:L*per/100);  


  
%------  fine the estimated C and d  -------- 

X_n = X - z;  
d_est = mean(X_n,2);                         %--- d is the mean data cloud 
U = X_n-d_est*ones(1,L);
OPTS.disp = 0;
[C_est D_est] = eigs(U*U',N-1,'LM',OPTS);    %--- C is the set which we desired

%---------------- step 4 : calculate the fitting error  ------------------- 

a_est = C_est'*( X - z - repmat(d_est,[1,L]) );
fit_error =  sum( sqrt(sum((X - (C_est*a_est + repmat(d_est,[1,L]) ) ).^2)));

%---------------- step 5 :  find the comparative tolerance ----------------
  
  tol=sum(sqrt(sum((zz-fit_error).^2)))/(sum(sqrt(sum(zz.^2))));
  zz = fit_error; 
  z(:,setdiff(1:L,error_index))=0;
  
  
   
  
 
 end
 


index=error_index;





