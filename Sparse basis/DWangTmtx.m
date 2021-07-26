function [ WN ] = DWangTmtx( N,dwt_type )  %离散W（Wang）变换基 不能用傅里叶矩阵
%DWangTmtx Summary of this function goes here
%   Detailed explanation goes here
%   N is the dimension of WN
%   dwt_type decides DWangT types
%   WN is the Discrete W Transform matrix
[k,n] = meshgrid(0:N-1);
if dwt_type==2
    alpha = 1/2;
    beta = 0;
elseif dwt_type==3
    alpha = 0;
    beta = 1/2;
elseif dwt_type==4
    alpha = 1/2;
    beta = 1/2;
else
    alpha = 0;
    beta = 0;
end
WN = sqrt(2/N) * sin(pi/4 + 2*pi/N * (n+alpha) .* (k+beta));
end