%=====================================================================
% Programmer: Chia-Hsiang Lin (Steven)
% E-mail: chiahsiang.steven.lin@gmail.com
% Date: January 6, 2021
% Website: https://sites.google.com/view/chiahsianglin/
% -------------------------------------------------------
% Reference:
% ``Nonnegative blind source separation for ill-conditioned mixtures via John ellipsoid,"
% accepted by IEEE Transactions on Neural Networks and Learning Systems, 2020.
%======================================================================
% function [M_est,X_est,time]= HiSun(Y,p)
%======================================================================
%  Input
%  Y is M-by-L data matrix, where M is the number of spectral bands and L is the number of pixels.
%  p is the number of endmembers.
%----------------------------------------------------------------------
%  Output
%  M_est is M-by-N mixing matrix whose columns are estimated endmember signatures.
%  X_est is N-by-L source matrix whose rows are estimated abundance maps.
%  time is the computation time (in secs).
%========================================================================
function [M_est,X_est,time]= HiSun(Y,p)
t0 = clock;
%% dimension reduction (DR)
[l,n]= size(Y);
d= mean(Y,2);
U= Y-d*ones(1,n);
R= U*U';
[eV,~]= eig(R);
C= eV(:,l-p+2:end);
Yd= C'*(Y-d*ones(1,n));
%% coompute H-polytope
disp(['start computing (B,h)'])
%
% %------------compute (B,h) for new data--------------
 [b,h,index] = vert2con_with_index(Yd'); B= b'; % method 2

disp(['end computing (B,h)'])
%% pre-conditioning
disp(['start computing (F,c)'])
[F,c]= algo1(B,h,Yd);

pre_Y= inv(F)*(Yd-c*ones(1,n));
disp(['end computing (F,c)'])
%% optimization
disp(['start blind source separation'])
[M_tilde,X_est]= algo3(pre_Y,p);
M_est= C*(F*M_tilde+c*ones(1,p))+d*ones(1,p);
disp(['end blind source separation'])
time = etime(clock,t0);

%% subprogram 1
%=====================================================================
% Programmer: Chia-Hsiang Lin (Steven)
% E-mail: chiahsiang.steven.lin@gmail.com
% Date: January 6, 2021
% Website: https://sites.google.com/view/chiahsianglin/
% -------------------------------------------------------
% Reference:
% ``Nonnegative blind source separation for ill-conditioned mixtures via John ellipsoid,"
% accepted by IEEE Transactions on Neural Networks and Learning Systems, 2020.
%======================================================================
% function [F,c]= algo1(B,h,Yd)
%======================================================================
%  Input
%  Yd collects data points as its columns (V-polytope).
%  (B,h) are the H-polytope representation of Yd (cf. Equation (3) of the reference paper).
%----------------------------------------------------------------------
%  Output
%  (F,c) describe the John ellipsoid of the convex hull of the input data points in Yd (cf. Equation (2) of the reference paper).
%========================================================================
function [F,c]= algo1(B,h,Yd)
%% parameters
p= size(B,1)+1; % number of endmembers
m= size(B,2); % number of facets
mu= 0.01; % ADMM parameter
iters= 80; % ADMM iterations
%% ADMM
% auxiliary matrices or vectors
bi_2norm_square_vec= (sum(B.*B))'; % column vector
I_plus_sum_bibiT_inverse= inv(eye(p-1)+B*B');
% smart initialization
[M0]= SPA(Yd,p);
[F,c]= JE_simplex(M0);
V0= zeros(p-1,p-1);
D0= V0;
Vi= zeros(2*(p-1),m);
Di= Vi;
for i=1:iters
    % update V0
    [Q,Lambda]= eig(0.5*((F-D0)+(F-D0)'));
    lambda_vec= diag(Lambda);
    lambda_vec= 0.5*(lambda_vec+sqrt((lambda_vec.^2)+(4/mu)));
    V0= Q*diag(lambda_vec)*Q';
    % update Vi
    Oi= [F*B; c*ones(1,m)]- Di;
    Vi= Oi; % case 1
    oi1_2norm_vec= (sqrt(sum((Oi(1:p-1,:)).^2)))'; % column vector
    hi_minus_biToi2_vec= h- (sum(B.*Oi(p:end,:)))'; % column vector
    index_case2= find( oi1_2norm_vec > hi_minus_biToi2_vec );
    t_star_vec= max(((bi_2norm_square_vec(index_case2).*oi1_2norm_vec(index_case2))+ hi_minus_biToi2_vec(index_case2))./(1+bi_2norm_square_vec(index_case2)),0);
    Vi(1:p-1,index_case2)= Vi(1:p-1,index_case2)*diag( t_star_vec./oi1_2norm_vec(index_case2) );
    Vi(p:end,index_case2)= Vi(p:end,index_case2)+ B(:,index_case2)*diag( (hi_minus_biToi2_vec(index_case2)-t_star_vec)./bi_2norm_square_vec(index_case2) );
    % update F
    Vi_plus_Di= Vi+Di;
    F= (V0+D0+Vi_plus_Di(1:p-1,:)*B')*I_plus_sum_bibiT_inverse;
    % update c
    c= mean(Vi_plus_Di(p:end,:),2);
    % update D0 and Di
    D0= D0+V0-F; Di= Vi_plus_Di- [F*B;c*ones(1,m)];
end
return;

%% subprogram 2
function [M,X]= algo3(pre_Y,p)
SPLX= eye(p)-(1/p)*ones(p,p);
[M0,~]= eig(SPLX*SPLX'); M0(:,1)=[]; M0= (((factorial(p-1))/(sqrt(p)))^(1/(p-1)))*M0';
alpha= 1*((1/(factorial(p-1)))*(p^(0.5*p))*((p-1)^(0.5*(p-1))))^(1/(p-1));
[M]= SPA(pre_Y,p); X= (1/p)*ones(p,size(pre_Y,2));
X= sunsal([M;ones(1,p)],[pre_Y;ones(1,size(pre_Y,2))],'POSITIVITY','yes','ADDONE','yes','lambda',0,'AL_ITERS',100,'TOL',1e-6,'X0',X);
for i=1:10,

    [V2,~,V1]= svd(pre_Y*X'*alpha*M0'); U= V1*V2';
    M= U'*alpha*M0;
    X = sunsal([M;ones(1,p)],[pre_Y;ones(1,size(pre_Y,2))],'POSITIVITY','yes','ADDONE','yes','lambda',0,'AL_ITERS',100,'TOL',1e-6,'X0',X);
end
return;

%% subprogram 3
function [A_est]= SPA(X,N)
con_tol= 1e-8; % the convergence tolence in SPA
num_SPA_itr= N; % number of iterations in post-processing of SPA
N_max= N; % max number of iterations
[M,L]= size(X);
d= mean(X,2);
U= X-d*ones(1,L);
[eV,~]= eig(U*U');
C= eV(:,M-N+2:end);
Xd= C'*(X-d*ones(1,L));
A_set=[]; Xd_t= [Xd;ones(1,L)]; index= [];
[val,ind]= max(sum(Xd_t.^2));
A_set= [A_set,Xd_t(:,ind)];
index= [index,ind];
for i=2:N,
    XX= (eye(N_max)-A_set*pinv(A_set))*Xd_t;
    [val,ind]= max(sum(XX.^2));
    A_set= [A_set,Xd_t(:,ind)];
    index= [index,ind];
end
alpha_tilde= Xd(:,index);
current_vol= det(alpha_tilde(:,1:N-1)-alpha_tilde(:,N)*ones(1,N-1));
for j=1:num_SPA_itr,
    for i=1:N,
        b(:,i)= compute_bi(alpha_tilde,i,N);
        b(:,i)= -b(:,i);
        [const,idx]= max(b(:,i)'*Xd);
        alpha_tilde(:,i)= Xd(:,idx);
    end
    new_vol= det(alpha_tilde(:,1:N-1)-alpha_tilde(:,N)*ones(1,N-1) );
    if (new_vol-current_vol)/current_vol<con_tol,
        break;
    end
end
A_est= C*alpha_tilde+d*ones(1,N);
return;

%% subprogram 4
function [bi]= compute_bi(a0,i,N)
Hindx= setdiff([1:N],[i]);
A_Hindx= a0(:,Hindx);
A_tilde_i= A_Hindx(:,1:N-2)-A_Hindx(:,N-1)*ones(1,N-2);
bi= A_Hindx(:,N-1)-a0(:,i);
bi= (eye(N-1)-A_tilde_i*(pinv(A_tilde_i'*A_tilde_i))*A_tilde_i')*bi;
bi= bi/norm(bi);
return;
