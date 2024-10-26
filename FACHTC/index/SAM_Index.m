function sam = SAM_Index(ground_truth, estimated,Mask)
[row,col,band]=size(ground_truth);
N=col*row;


M_all_Pre=ones(row,col);
for i=1:band
    M_all_Pre=M_all_Pre.*Mask(:,:,i);%2D mask
end
[m,n]=find(M_all_Pre==0);
pixel_cut=(n-1)*row+m;


% SAM
y = reshape(ground_truth,N,band)';x = reshape(estimated,N,band)';
y = y(:,pixel_cut);x = x(:,pixel_cut);
num = sum(x .* y, 1);
den = sqrt(sum(x.^2, 1) .* sum(y.^2, 1));
index=find(den==0);
den(index)=[];
num(index)=[];

sam = sum(sum(acosd(num ./ den)))/(N-length(index));
end
