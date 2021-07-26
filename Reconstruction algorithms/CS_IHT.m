function hat_x=CS_IHT(y,T_Mat,m)  %%%% 这个有点问题！！！！！！
% y=T_Mat*x, T_Mat is n-by-m
% y - measurements
% T_Mat - combination of random matrix and sparse representation basis
% m - size of the original signal
% the sparsity is length(y)/4

hat_x_tp=zeros(m,1);         % initialization with the size of original 
s=floor(length(y)/4);        % sparsity
u=0.5;                       % impact factor

% T_Mat=T_Mat/sqrt(sum(sum(T_Mat.^2))); % normalizae the whole matrix

for times=1:s
    
    x_increase=T_Mat'*(y-T_Mat*hat_x_tp);
    
    hat_x=hat_x_tp+u*x_increase;
    
    [val,pos]=sort((hat_x),'descend');  % why? worse performance with abs()
    
    hat_x(pos(s+1:end))=0;   % thresholding, keeping the larges s elements

    hat_x_tp=hat_x;          % update

end
