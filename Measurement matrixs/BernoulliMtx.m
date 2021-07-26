function [ Phi ] = BernoulliMtx( M,N ) %生成伯努利测量矩阵
%BernoulliMtx Summary of this function goes here
%   Generate Bernoulli matrix 
%   M -- RowNumber
%   N -- ColumnNumber
%   Phi -- The Bernoulli matrix

%% (1)Generate Bernoulli matrix(The first kind)
% 1--P=0.5   -1--P=0.5
Phi = randi([0,1],M,N);%If your MATLAB version is too low,please use randint instead
%生成数值在[0,1]区间内整数的 M*N矩阵
    Phi(Phi==0) = -1;
    %Phi = Phi/sqrt(M);
% %% (2)Generate Bernoulli matrix(The second kind)
% % 1--P=1/6   -1--P=1/6  0--2/3
%     Phi = randi([-1,4],M,N);%If your MATLAB version is too low,please use randint instead
%     Phi(Phi==2) = 0;%P=1/6
%     Phi(Phi==3) = 0;%P=1/6
%     Phi(Phi==4) = 0;%P=1/6
%     %Phi = Phi*sqrt(3/M);
end
