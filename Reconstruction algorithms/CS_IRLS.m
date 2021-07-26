function hat_x=CS_IRLS(y,T_Mat,m)
% y=T_Mat*x, T_Mat is n-by-m
% y - measurements
% T_Mat - combination of random matrix and sparse representation basis
% m - size of the original signal
% the sparsity is length(y)/4

%hat_x_tp=(T_Mat'*T_Mat)^(-1)*T_Mat'*y; % initialization, leading to warning
hat_x_tp=T_Mat'*y; 

epsilong=1;

p=1;                                   % solution for l-norm p
times=1;
while (epsilong>10e-9) && (times<length(y)/4)
    
    weight=(hat_x_tp.^2+epsilong).^(p/2-1); 
    Q_Mat=diag(1./weight,0);
    
    hat_x=Q_Mat*T_Mat'*inv(T_Mat*Q_Mat*T_Mat')*y;
    
    if(norm(hat_x-hat_x_tp,2) < sqrt(epsilong)/100)
        epsilong=epsilong/10;
    end
    
    hat_x_tp=hat_x;
    times=times+1;
end