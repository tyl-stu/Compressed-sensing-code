function ww=DWT(N)  %生成正交小波变换矩阵，正交小波变换基

[h,g]= wfilters('sym8','d');       %  分解低通和高通滤波器  

% N=256;                           %  矩阵维数(大小为2的整数幂次)
L=length(h);                       %  滤波器长度
rank_max=log2(N);                  %  最大层数
rank_min=double(int8(log2(L)))+1;  %  最小层数
ww=1;   %  预处理矩阵

%  矩阵构造
for jj=rank_min:rank_max
    
    nn=2^jj;
    
    %  构造向量
    p1_0=sparse([h,zeros(1,nn-L)]);
    p2_0=sparse([g,zeros(1,nn-L)]);
    
    %  向量圆周移位
    for ii=1:nn/2
        p1(ii,:)=circshift(p1_0',2*(ii-1))';
        p2(ii,:)=circshift(p2_0',2*(ii-1))';
    end
    
    %  构造正交矩阵
    w1=[p1;p2];
    mm=2^rank_max-length(w1);
    w=sparse([w1,zeros(length(w1),mm);zeros(mm,length(w1)),eye(mm,mm)]);
    ww=ww*w;
    
    clear p1;clear p2;
end

