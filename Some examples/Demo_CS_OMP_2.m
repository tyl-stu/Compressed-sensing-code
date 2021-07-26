%压缩感知重构算法测试
% 压缩感知重构算法之正交匹配追踪(OMP)！！！！！！！！！！！！！！！！！！！！！！
%%%% 恢复残差：ans = 3.0390 太大！！！！！！！！！！！！！！！！！
clear all;close all;clc;
M = 64*2;%观测值个数
N = 256*2;%信号x的长度
K = 10*2;%信号x的稀疏度
Index_K = randperm(N);
x = zeros(N,1);
% x(Index_K(1:K)) = 5*randn(K,1);%x为K稀疏的，且位置是随机的
x(Index_K(1:K)) = 1;
Psi = eye(N);%x本身是稀疏的，定义稀疏矩阵为单位阵x=Psi*theta
Phi = randn(M,N);%测量矩阵为高斯矩阵
A = Phi * Psi;%传感矩阵
y = Phi * x;%得到观测向量y
%% 恢复重构信号x
tic
theta = CS_OMP(y,A,K);
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
% [n1,r1] = biterr(x,round(x_r));