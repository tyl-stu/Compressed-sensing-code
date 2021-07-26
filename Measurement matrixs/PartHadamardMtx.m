function [ Phi ] = PartHadamardMtx( M,N )   %生成部分哈达玛矩阵
%PartHadamardMtx Summary of this function goes here
%   Generate part Hadamard matrix 
%   M -- RowNumber
%   N -- ColumnNumber
%   Phi -- The part Hadamard matrix

%% parameter initialization
%Because the MATLAB function hadamard handles only the cases where n, n/12,
%or n/20 is a power of 2
    L_t = max(M,N);%Maybe L_t does not meet requirement of function hadamard
    L_t1 = (12 - mod(L_t,12)) + L_t; % mod() 相当于 rem()
    L_t2 = (20 - mod(L_t,20)) + L_t; 
    L_t3 = 2^ceil(log2(L_t)); % ceil() 朝正无穷大四舍五入

    L = min([L_t1,L_t2,L_t3]);%Get the minimum L
%% Generate part Hadamard matrix   
    Phi = [];
Phi_t = hadamard(L);% hadamard() 每行列皆正交
RowIndex = randperm(L);%行索引
Phi_t_r = Phi_t(RowIndex(1:M),:);% 取 RowIndex 的前 M 个数值作为行索引，取出原Phi_t中对应的行
ColIndex = randperm(L); % 列索引
Phi = Phi_t_r(:,ColIndex(1:N)); % 取 ColIndex 的前 N 个数值作为列索引，取出对应的列
end