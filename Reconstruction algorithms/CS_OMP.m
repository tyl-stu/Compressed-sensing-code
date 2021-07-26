function [ theta ] = CS_OMP( y,A,t )
%CS_OMP Summary of this function goes here
%Version: 1.0 written by jbb0523 @2015-04-18
%   Detailed explanation goes here
%   y = Phi * x
%   x = Psi * theta
%	y = Phi*Psi * theta
%   令 A = Phi*Psi, 则y=A*theta
%   现在已知y和A，求theta
    [y_rows,y_columns] = size(y);
    if y_rows<y_columns
        y = y';%y should be a column vector
    end
    [M,N] = size(A);%传感矩阵A为M*N矩阵
    theta = zeros(N,1);%用来存储恢复的theta(列向量)
    
    At = zeros(M,t);%用来迭代过程中存储A被选择的列
    Pos_theta = zeros(1,t);%用来迭代过程中存储A被选择的列序号
    r_n = y;%初始化残差(residual)为y
    for ii=1:t%迭代t次，t为输入参数
        product = A'*r_n;%传感矩阵A各列与残差的内积
        [val,pos] = max(abs(product));%找到最大内积绝对值，即与残差最相关的列
        At(:,ii) = A(:,pos);%存储这一列
        Pos_theta(ii) = pos;%存储这一列的序号
        A(:,pos) = zeros(M,1);%清零A的这一列，其实此行可以不要，因为它与残差正交
        %y=At(:,1:ii)*theta，以下求theta的最小二乘解(Least Square)
        theta_ls = (At(:,1:ii)'*At(:,1:ii))^(-1)*At(:,1:ii)'*y;%最小二乘解
        %At(:,1:ii)*theta_ls是y在At(:,1:ii)列空间上的正交投影
        r_n = y - At(:,1:ii)*theta_ls;%更新残差        
    end
    theta(Pos_theta)=theta_ls;%恢复出的theta
end
