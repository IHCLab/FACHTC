function [mssim,psnr]=assessment(ref,tar)

%% SSIM
mssim = MSSIM_3D(ref,tar);
%% PSNR
[~,~,bands]=size(ref);
ref=reshape(ref,[],bands);
tar=reshape(tar,[],bands);
msr=mean((ref-tar).^2,1);
max2=max(tar,[],1).^2;
psnr=mean(10*log10(max2./msr));
end