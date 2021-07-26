function [ Phi ] = CirculantMtx( M,N )  %Éú³ÉÑ­»·¾ØÕó
%CirculantMtx Summary of this function goes here
%   Generate Circulant matrix 
%   M -- RowNumber
%   N -- ColumnNumber
%   Phi -- The Circulant matrix

%% Generate a random vector
%     %(1)Gauss
%     u = randn(1,N);
    %(2)Bernoulli
    u = randi([0,1],1,N);
    u(u==0) = -1;
%% Generate Circulant matrix   
    Phi_t = toeplitz(circshift(u,[1,1]),fliplr(u(1:N)));
    Phi = Phi_t(1:M,:);
end