function [Y,Mask,Label]= Nearest_entry_inpainting(Y,Mask,thre)
[col,row,band] = size(Y);
N =col*row;
Y=reshape(Y,N,band)';%3-D to 2-D
Mask = reshape(Mask,N,band)';%3-D to 2-D
Label = ones(N,band)';
m_val = sum(Y,2)./sum(Mask,2); 
for kk = 1:N
    L_tmp = 1:band;L_tmp = L_tmp';
    P_tmp = Y(:,kk);
    M_tmp = Mask(:,kk);
    Cp = find(M_tmp~=0);
    Ip = find(M_tmp==0);
    IpM = Ip*ones(1,length(Cp));
    CpM = ones(length(Ip),1)*Cp';
    M = abs(IpM-CpM);
    Dis_p = min(M,[],2);
    Over_p = (Dis_p>thre);
    Dis_p(Over_p) = [];
    Ip(Over_p) = [];
    PP_up_tmp = Ip-Dis_p;
    PP_down_tmp = Ip+Dis_p;
    PP_up_tmp(PP_up_tmp<1)=PP_down_tmp(PP_up_tmp<1);
    PP_up = M_tmp((PP_up_tmp)');
    PP = PP_up.*PP_up_tmp+(1-PP_up).*PP_down_tmp;
    P_tmp(Ip) = P_tmp(PP).*(m_val(Ip)./m_val(PP));
    M_tmp(Ip) = M_tmp(PP);
    L_tmp(Ip) = L_tmp(PP);
    Y(:,kk) = P_tmp;
    Mask(:,kk) = M_tmp;
    Label(:,kk) = L_tmp;
end
Y=reshape(Y',col,row,band);%3-D to 2-D
Mask=reshape(Mask',col,row,band);%3-D to 2-D

end