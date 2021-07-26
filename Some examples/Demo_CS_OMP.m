function Demo_CS_OMP()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the DCT basis is selected as the sparse representation dictionary
% instead of seting the whole image as a vector, I process the image in the
% fashion of column-by-column, so as to reduce the complexity.
% 选择DCT基作为稀疏表示字典，而不是将整个图像设为向量，采用逐列的方式对图像进行处理，以降低复杂度

% Author: Chengfu Huo, roy@mail.ustc.edu.cn, http://home.ustc.edu.cn/~roy
% Reference: J. Tropp and A. Gilbert, “Signal Recovery from Random
% Measurements via Orthogonal Matching Pursuit,” 2007.
% 文献：基于正交匹配追踪的随机测量信号恢复 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------ read in the image --------------
img=imread('lena.bmp');     % testing image
img=double(img);
[height,width]=size(img);


%------------ form the measurement matrix and base matrix ---------------  形成测量矩阵和基矩阵
Phi=randn(floor(height/3),width);  % only keep one third of the original data
Phi = Phi./repmat(sqrt(sum(Phi.^2,1)),[floor(height/3),1]); % normalize each column 归一化每一列
% B = repmat(A,n) 返回一个数组，该数组在其行维度和列维度包含 A 的 n 个副本。A 为矩阵时，B 大小为 size(A)*n。
mat_dct_1d=zeros(256,256);  % building the DCT basis (corresponding to each column) 建立DCT基础(对应每一列)
for k=0:1:255
    dct_1d=cos([0:1:255]'*k*pi/256);
    if k>0
        dct_1d=dct_1d-mean(dct_1d); % mean() -- 求取数组的平均数
    end
    mat_dct_1d(:,k+1)=dct_1d/norm(dct_1d);
end


%--------- projection ---------
img_cs_1d=Phi*img;          % treat each column as a independent signal 每一列都作为一个独立的信号，Phi为测量矩阵
                            % img_cs_1d=Phi*img,即为 y=Φx，Φ为测量矩阵，x为原始信号，y为观测值

%-------- recover using omp ------------
%%%%%%%由 y=Θs，感知矩阵Θ，可获得稀疏信号的逼近值，再利用x’= Ψs’，Ψ为稀疏基，可获得原始信号的逼近值，即为目标结果

sparse_rec_1d=zeros(height,width);  % sparse - 稀疏  ，height,width 为读入照片的行、列 ,256*256
Theta_1d=Phi*mat_dct_1d;  % 相当于 Theta_1d = 测量矩阵 * DCT稀疏基 ， 即 Theta_1d 为感知矩阵或回复矩阵
for i=1:width  % 1:256
    column_rec=cs_omp(img_cs_1d(:,i),Theta_1d,height);
    sparse_rec_1d(:,i)=column_rec';         % sparse representation 稀疏表示，只有信号是K稀疏的（且K小于M远远小于N），
    %才有可能在观测M个观测值时，从K个较大的系数重建原始长度为N的信号，
    %被恢复的数据矩阵大小或感知矩阵的大小为 M * N
end %此处程序是利用观测值 y 与 感知矩阵求出了稀疏信号的逼近值
img_rec_1d=mat_dct_1d*sparse_rec_1d;          % inverse transform 逆变换，恢复后的数据 img_rec_1d = DCT稀疏基 * 稀疏恢复数据


%------------ show the results --------------------
figure(1)
subplot(2,2,1),imagesc(img),title('original image')
subplot(2,2,2),imagesc(Phi),title('measurement mat')   % 测量矩阵
subplot(2,2,3),imagesc(mat_dct_1d),title('1d dct mat') % 稀疏基
psnr = 20*log10(255/sqrt(mean((img(:)-img_rec_1d(:)).^2)))
subplot(2,2,4),imagesc(img_rec_1d),title(strcat('1d rec img ',num2str(psnr),'dB'))

disp('over');


%************************************************************************%
function hat_x=cs_omp(y,T_Mat,m)
    % y=T_Mat*s, T_Mat is n-by-m ,T_Mat是为感知矩阵
    % y - measurements
    % T_Mat - combination of random matrix and sparse representation basis
    % m - size of the original signal，或者说即为原始图像的行数
    % the sparsity is length(y)/4 稀疏度是观测值长度的四分之一

    n=length(y); %此也为输入信号的长度
    s=floor(n/4);                                     %  测量值维数或稀疏度
    hat_x=zeros(1,m);                                 %  待重构的谱域(变换域)向量
    Aug_t=[];                                         %  增量矩阵(初始值为空矩阵)
    r_n=y;                                            %  残差值 残差在数理统计中是指实际观察值与估计值(拟合值)之间的差

    for times=1:s                                  %  迭代次数(稀疏度是测量的1/4)

        product=abs(T_Mat'*r_n);   % y = T_Mat*s,T_Mat 是 M*N ,y 是 M*1  ==> product 或 s 为 N*1

        [val,pos]=max(product);                       %  最大投影系数对应的位置
        Aug_t=[Aug_t,T_Mat(:,pos)];                   %  矩阵扩充
        T_Mat(:,pos)=zeros(n,1);                      %  选中的列置零（实质上应该去掉，为了简单将其置零）
        aug_x=(Aug_t'*Aug_t)^(-1)*Aug_t'*y;           %  最小二乘,使残差最小,稀疏信号的最小二乘估计量
        r_n=y-Aug_t*aug_x;                            %  残差
        pos_array(times)=pos;                         %  纪录最大投影系数的位置

    end
    hat_x(pos_array)=aug_x;                           %  重构的向量



