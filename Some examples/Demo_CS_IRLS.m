%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the DCT basis is selected as the sparse representation dictionary
% instead of seting the whole image as a vector, I process the image in the
% fashion of column-by-column, so as to reduce the complexity.

% Author: Chengfu Huo, roy@mail.ustc.edu.cn, http://home.ustc.edu.cn/~roy
% Reference: R. Chartrand and W. Yin, ¡°Iteratively Reweighted Algorithms
% for Compressed Sensing,¡± 2008.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------ read in the image --------------
clear;clc;
img=imread('lena.bmp');     % testing image
img=double(img);
[height,width]=size(img);


%------------ form the measurement matrix and base matrix ---------------
Phi=randn(floor(height/3),width);  % only keep one third of the original data  
Phi = Phi./repmat(sqrt(sum(Phi.^2,1)),[floor(height/3),1]); % normalize each column

Psi=zeros(256,256);  % building the DCT basis (corresponding to each column)
for k=0:1:255 
    dct_1d=cos([0:1:255]'*k*pi/256);
    if k>0
        dct_1d=dct_1d-mean(dct_1d); 
    end;
    Psi(:,k+1)=dct_1d/norm(dct_1d);
end


%--------- projection ---------
y=Phi*img;          % treat each column as a independent signal


%-------- recover using omp ------------
sparse_rec_1d=zeros(height,width);            
Theta_1d=Phi*Psi;
for i=1:width
    column_rec=CS_IRLS(y(:,i),Theta_1d,height);
    sparse_rec_1d(:,i)=column_rec';           % sparse representation
end
img_rec_1d=Psi*sparse_rec_1d;          % inverse transform


%------------ show the results --------------------
figure(1)
subplot(2,2,1),imagesc(img),title('original image')
subplot(2,2,2),imagesc(Phi),title('measurement mat')
subplot(2,2,3),imagesc(Psi),title('1d dct mat')
psnr = 20*log10(255/sqrt(mean((img(:)-img_rec_1d(:)).^2)))
subplot(2,2,4),imagesc(img_rec_1d),title(strcat('1d rec img ',num2str(psnr),'dB'))

disp('over')







