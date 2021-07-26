
% 压缩感知重构算法之迭代硬阈值(Iterative Hard Thresholding,IHT)！！！！！！
clear all;close all;clc;      
M = 64;%观测值个数      
N = 256;%信号x的长度      
K = 30;%信号x的稀疏度      
Index_K = randperm(N);      
x = zeros(N,1);      
% x(Index_K(1:K)) = 5*randn(K,1);%x为K稀疏的，且位置是随机的    
x(Index_K(1:K)) = 1;
Psi = eye(N);%x本身是稀疏的，定义稀疏矩阵为单位阵x=Psi*theta      
Phi = randn(M,N);%测量矩阵为高斯矩阵  
Phi = orth(Phi')';    
A = Phi * Psi;%传感矩阵    
% sigma = 0.005;    
% e = sigma*randn(M,1);  
% y = Phi * x + e;%得到观测向量y      
y = Phi * x;%得到观测向量y    
%% 恢复重构信号x      
tic  
theta = CS_IHT_2(y,A,K); 
% theta = IHT_Basic(y,A,K); 
% theta = cs_iht(y,A,size(A,2));
% theta = hard_l0_Mterm(y,A,size(A,2),round(1.5*K),'verbose',true);
x_r = Psi * theta;% x=Psi * theta      
toc      
%% 绘图      
figure;      
plot(x_r,'k.-');%绘出x的恢复信号      
hold on;      
plot(x,'r');%绘出原信号x      
hold off;      
legend('Recovery','Original')      
fprintf('\n恢复残差：');      
norm(x_r-x)%恢复残差  
[n1,r1] = biterr(x,round(x_r));