function Demo_CS_SP()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the DCT basis is selected as the sparse representation dictionary
% instead of seting the whole image as a vector, I process the image in the
% fashion of column-by-column, so as to reduce the complexity.

% Author: Chengfu Huo, roy@mail.ustc.edu.cn, http://home.ustc.edu.cn/~roy
% Reference: W. Dai and O. Milenkovic, ¡°Subspace Pursuit for Compressive
% Sensing Signal Reconstruction,¡± 2009.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------ read in the image --------------
img=imread('lena.bmp');     % testing image
img=double(img);
[height,width]=size(img);


%------------ form the measurement matrix and base matrix ---------------
Phi=randn(floor(height/3),width);  % only keep one third of the original data  
Phi = Phi./repmat(sqrt(sum(Phi.^2,1)),[floor(height/3),1]); % normalize each column


mat_dct_1d=zeros(256,256);  % building the DCT basis (corresponding to each column)
for k=0:1:255 
    dct_1d=cos([0:1:255]'*k*pi/256);
    if k>0
        dct_1d=dct_1d-mean(dct_1d); 
    end;
    mat_dct_1d(:,k+1)=dct_1d/norm(dct_1d);
end


%--------- projection ---------
img_cs_1d=Phi*img;          % treat each column as a independent signal


%-------- recover using omp ------------
sparse_rec_1d=zeros(height,width);            
Theta_1d=Phi*mat_dct_1d;
for i=1:width
    column_rec=cs_sp(img_cs_1d(:,i),Theta_1d,height);
    sparse_rec_1d(:,i)=column_rec';           % sparse representation
end
img_rec_1d=mat_dct_1d*sparse_rec_1d;          % inverse transform


%------------ show the results --------------------
figure(1)
subplot(2,2,1),imagesc(img),title('original image')
subplot(2,2,2),imagesc(Phi),title('measurement mat')
subplot(2,2,3),imagesc(mat_dct_1d),title('1d dct mat')
psnr = 20*log10(255/sqrt(mean((img(:)-img_rec_1d(:)).^2)))
subplot(2,2,4),imagesc(img_rec_1d),title(strcat('1d rec img ',num2str(psnr),'dB'))

disp('over')


%************************************************************************%
function hat_x=cs_sp(y,T_Mat,m)
% y=T_Mat*x, T_Mat is n-by-m
% y - measurements
% T_Mat - combination of random matrix and sparse representation basis
% m - size of the original signal
% the sparsity is length(y)/4

n=length(y);                           % length of measurements
s=floor(n/4);                                 % sparsity                  
r_n=y;                                 % initial residuals

sig_pos_lt=[];                         % significant pos for last time iteration

for times=1:s                          % number of iterations
    
    product=abs(T_Mat'*r_n);
    [val,pos]=sort(product,'descend');
    sig_pos_cr=pos(1:s);             % significant pos for curretn iteration
    
    sig_pos=union(sig_pos_cr,sig_pos_lt);
    
    Aug_t=T_Mat(:,sig_pos);            % current selected entries of T_Mat 
    
    aug_x_cr=zeros(m,1);               
    aug_x_cr(sig_pos)=(Aug_t'*Aug_t)^(-1)*Aug_t'*y;  % temp recovered x (sparse)
    
    [val,pos]=sort(abs(aug_x_cr),'descend');
    
    hat_x=zeros(1,m);
    hat_x(pos(1:s))=aug_x_cr(pos(1:s));% recovered x with s sparsity  
    
    sig_pos_lt=pos(1:s);               % refresh the significant positions
    
    r_n=y-T_Mat*hat_x';
end
             














