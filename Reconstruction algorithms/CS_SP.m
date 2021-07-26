function [ theta ] = CS_SP( y,A,K )
%CS_SP Summary of this function goes here
%Version: 1.0 written by jbb0523 @2015-05-01
%   Detailed explanation goes here
%   y = Phi * x
%   x = Psi * theta
%	y = Phi*Psi * theta
%   令 A = Phi*Psi, 则y=A*theta
%   K is the sparsity level
%   现在已知y和A，求theta
%   Reference:Dai W，Milenkovic O．Subspace pursuit for compressive sensing
%   signal reconstruction[J]．IEEE Transactions on Information Theory，
%   2009，55(5)：2230-2249.
    [y_rows,y_columns] = size(y);
    if y_rows<y_columns
        y = y';%y should be a column vector
    end
    [M,N] = size(A);%传感矩阵A为M*N矩阵
    theta = zeros(N,1);%用来存储恢复的theta(列向量)
    Pos_theta = [];%用来迭代过程中存储A被选择的列序号
    r_n = y;%初始化残差(residual)为y
    for kk=1:K%最多迭代K次
        %(1) Identification
        product = A'*r_n;%传感矩阵A各列与残差的内积
        [val,pos]=sort(abs(product),'descend');
        Js = pos(1:K);%选出内积值最大的K列
        %(2) Support Merger
        Is = union(Pos_theta,Js);%Pos_theta与Js并集
        %(3) Estimation
        %At的行数要大于列数，此为最小二乘的基础(列线性无关)
        if length(Is)<=M
            At = A(:,Is);%将A的这几列组成矩阵At
        else%At的列数大于行数，列必为线性相关的,At'*At将不可逆
            break;%跳出for循环
        end
        %y=At*theta，以下求theta的最小二乘解(Least Square)
        theta_ls = (At'*At)^(-1)*At'*y;%最小二乘解
        %(4) Pruning
        [val,pos]=sort(abs(theta_ls),'descend');
        %(5) Sample Update
        Pos_theta = Is(pos(1:K));
        theta_ls = theta_ls(pos(1:K));
        %At(:,pos(1:K))*theta_ls是y在At(:,pos(1:K))列空间上的正交投影
        r_n = y - At(:,pos(1:K))*theta_ls;%更新残差 
        if norm(r_n)<1e-6%Repeat the steps until r=0
            break;%跳出for循环
        end
    end
    theta(Pos_theta)=theta_ls;%恢复出的theta
end
