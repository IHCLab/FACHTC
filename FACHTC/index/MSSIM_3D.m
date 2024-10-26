function k = MSSIM_3D(Y3D_ref, Y3D_rec) 

%Input : (1) Y3D_ref: 3D clean HSI data
%        (2) Y3D_rec: 3D reconstructed HSI data computed by algorithm
%Output: (1) k: mean ssim of each bands


[row, col, bands] = size(Y3D_ref);%row, column, band
K = [0.01 0.03];

%window = fspecial('gaussian', 11, 1.5);
window = ones(64);

Y2D_rec=reshape(Y3D_rec,[],bands)';
Y2D_ref=reshape(Y3D_ref,[],bands)';

% input size : Bands*observation
L=max(max(Y2D_rec(:),max(Y2D_ref(:))));
for i=1:bands 
    k_tmp(i)  = ssim_index(reshape(Y2D_ref(i,:),[row, col]), reshape(Y2D_rec(i,:),[row, col]), K, window, L);
end
k=mean(k_tmp);