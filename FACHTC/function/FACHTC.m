%=====================================================================
% Programmer: Yangrui Liu
% E-mail: chiahsiang.steven.lin@gmail.com; q36083030@gs.ncku.edu.tw; tonyquek@sutd.edu.sg;
% Date: March 23, 2023
% Website: https://sites.google.com/view/chiahsianglin/
% -------------------------------------------------------
% Reference:
% Chia-Hsiang Lin, Yangrui Liu, Chong-Yung Chi, Chih-Chung Hsu, Hsuan Ren and Tony Q.S.Quek, 
% "Hyperspectral Tensor Completion Using Low-Rank Modeling and Convex Functional Analysis,"
% accepted by IEEE Transactions on Neural Networks and Learning Systems, 2023
%======================================================================
% function image_out=FACHTC(Y,p,Mask)
%======================================================================
%  Input
%  Y is M-by-L data matrix, where M is the number of spectral bands and L is the number of pixels.
%  p is the number of endmembers.
%  Mask is 3-D binary tensor with L1*L2*M size, where L1 and L2 represent spatial size and M is the number of spectral bands.
%----------------------------------------------------------------------
%  Output
%  image_out is the recovered hyperspectral data cube of dimension L1*L2*M.
%========================================================================
function image_out=FACHTC(Y,p,Mask)
[row,col,band]=size(Y);
N=col*row;
M_all_Pre=ones(row,col);
for i=1:band
    M_all_Pre=M_all_Pre.*Mask(:,:,i);% 2D mask
end

thre = 191; % How far pixels can be used for nearest-inpainting
%% edited by Yangrui on 2021/12/12
if sum(M_all_Pre(:))==0
    [Y,Mask,~]= Nearest_entry_inpainting(Y,Mask,thre);
end
%%
M_all=ones(row,col);
ib = [];
for i=1:band
    if sum(sum(Mask(:,:,i))) <N
        ib = [ib,i];
    end
    M_all=M_all.*Mask(:,:,i);% 2D mask
end

[m,n]=find(M_all==0);
pixel_cut=(n-1)*row+m; % Corrupted pixels
band_cut=ib;% Corrupted bands
Y=reshape(Y,N,band)';% 3-D to 2-D

X1=Y;
X2=Y;

X1(:,pixel_cut')=[];% Remove corrupted pixels
[C d O_index] = RASF(X1,p, 10/225);
X1(: , O_index) = [] ;% Remove outliers
[M_est,X_est,time]= HiSun(X1,p); % Use conplete pixels to compute signatures
A_est =  M_est;

A_i =  A_est;
A_i(band_cut,:)=[];
X2(band_cut,:)=[];
S=(1/p)*ones(p,size(X2,2)); 

% Compute abundance using complete bands
S_est=sunsal(A_i,X2,'POSITIVITY','yes','ADDONE','no','lambda',0,'AL_ITERS',100,'TOL',1e-6,'X0',(1/p)*ones(p,size(Y,2)));


X_out=A_est*S_est;% Reconstruction

image_out=reshape(X_out',row,col,band);% 2D to 3D data
end