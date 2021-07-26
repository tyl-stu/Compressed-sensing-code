function [ theta ] = CS_ROMP( y,A,K )
%CS_ROMP Summary of this function goes here
%Version: 1.0 written by jbb0523 @2015-04-24
%   Detailed explanation goes here
%   y = Phi * x
%   x = Psi * theta
%	y = Phi*Psi * theta
%   令 A = Phi*Psi, 则y=A*theta
%   现在已知y和A，求theta
%   Reference:Needell D，Vershynin R．Signal recovery from incomplete and
%   inaccurate measurements via regularized orthogonal matching pursuit[J]．
%   IEEE Journal on Selected Topics in Signal Processing，2010，4(2)：310—316.
    [y_rows,y_columns] = size(y);
    if y_rows<y_columns
        y = y';%y should be a column vector
    end
    [M,N] = size(A);%传感矩阵A为M*N矩阵
    theta = zeros(N,1);%用来存储恢复的theta(列向量),即重构稀疏信号
    At = zeros(M,3*K);%用来迭代过程中存储A被选择的列
    Pos_theta = zeros(1,2*K);%用来迭代过程中存储A被选择的列序号
    Index = 0;
    r_n = y;%初始化残差(residual)为y
    %Repeat the following steps K times(or until |I|>=2K)
    for ii=1:K%迭代K次
        product = A'*r_n;%传感矩阵A各列与残差的内积
        %[val,pos] = max(abs(product));%找到最大内积绝对值，即与残差最相关的列
        [val,pos] = Regularize(product,K);%按正则化规则选择原子
        At(:,Index+1:Index+length(pos)) = A(:,pos);%存储这几列
        Pos_theta(Index+1:Index+length(pos)) = pos;%存储这几列的序号
        if Index+length(pos)<=M%At的行数大于列数，此为最小二乘的基础(列线性无关)
            Index = Index+length(pos);%更新Index，为下次循环做准备
        else%At的列数大于行数，列必为线性相关的,At(:,1:Index)'*At(:,1:Index)将不可逆
            break;%跳出for循环
        end
        A(:,pos) = zeros(M,length(pos));%清零A的这几列(其实此行可以不要因为它们与残差正交)
        %y=At(:,1:Index)*theta，以下求theta的最小二乘解(Least Square)
        theta_ls = (At(:,1:Index)'*At(:,1:Index))^(-1)*At(:,1:Index)'*y;%最小二乘解
        %At(:,1:Index)*theta_ls是y在At(:,1:Index)列空间上的正交投影
        r_n = y - At(:,1:Index)*theta_ls;%更新残差
        if norm(r_n)<1e-6%Repeat the steps until r=0
            break;%跳出for循环
        end
        if Index>=2*K%or until |I|>=2K
            break;%跳出for循环
        end
    end
    theta(Pos_theta(1:Index))=theta_ls;%恢复出的theta
end

%% 选择u中适合的列集合
function [val,pos] = Regularize(product,Kin)
%Regularize Summary of this function goes here
%   Detailed explanation goes here
%   product = A'*r_n;%传感矩阵A各列与残差的内积 ，维护为 N*1
%   K为稀疏度
%   pos为选出的各列序号
%   val为选出的各列与残差的内积值
%   Reference:Needell D，Vershynin R. Uniform uncertainty principle and
%   signal recovery via regularized orthogonal matching pursuit. 
%   Foundations of Computational Mathematics, 2009,9(3): 317-334.  
    productabs = abs(product);%取绝对值
    [productdes,indexproductdes] = sort(productabs,'descend');%降序排列
    for ii = length(productdes):-1:1 %从小到大找非零值
        if productdes(ii)>1e-6  %判断productdes中非零值个数
            break;
        end
    end
    %Identify:Choose a set J of the K biggest coordinates
    if ii>=Kin  %从内积向量u中选择k个最大值
        J = indexproductdes(1:Kin);%集合J，保存对应内基向量的序号
        Jval = productdes(1:Kin);%集合J对应的序列值，保存对应内基向量的值
        K = Kin;
    else%or all of its nonzero coordinates,whichever is smaller，选择所有非零值
        J = indexproductdes(1:ii);%集合J
        Jval = productdes(1:ii);%集合J对应的序列值
        K = ii;
    end
    %Regularize:Among all subsets J0∈J with comparable coordinates
    MaxE = -1;%循环过程中存储最大能量值
    for kk = 1:K
        J0_tmp = zeros(1,K);iJ0 = 1;
        J0_tmp(iJ0) = J(kk);%以J(kk)为本次寻找J0的基准(最大值)
        Energy = Jval(kk)^2;%本次寻找J0的能量
        for mm = kk+1:K
            if Jval(kk)<2*Jval(mm)%找到符合|u(i)|<=2|u(j)|的
                iJ0 = iJ0 + 1;%J0自变量增1
                J0_tmp(iJ0) = J(mm);%更新J0
                Energy = Energy + Jval(mm)^2;%更新能量
            else%不符合|u(i)|<=2|u(j)|的
                break;%跳出本轮寻找，因为后面更小的值也不会符合要求
            end
        end
        if Energy>MaxE%本次所得J0的能量大于前一组
            J0 = J0_tmp(1:iJ0);%更新J0
            MaxE = Energy;%更新MaxE，为下次循环做准备
        end
    end
    pos = J0;
    val = productabs(J0);
end