function Phi = PartFourierMtx(M,N)%生成傅里叶测量矩阵
    % N = 6; M = 3;
    Phi_t = fft(eye(N,N))/sqrt(N);%Fourier matrix ,生成傅里叶正交矩阵
    RowIndex = randperm(N);
    Phi = Phi_t(RowIndex(1:M),:);%Select M rows randomly
%     normalization
    for ii = 1:N
        Phi(:,ii) = Phi(:,ii)/norm(Phi(:,ii));
    end
end