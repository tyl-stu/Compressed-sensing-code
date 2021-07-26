function [ Phi ] = SparseRandomMtx( M,N,d )  % Éú³ÉÏ¡ÊèËæ»ú¾ØÕó
%SparseRandomMtx Summary of this function goes here
%   Generate SparseRandom matrix 
%   M -- RowNumber
%   N -- ColumnNumber
%   d -- The number of '1' in every column,d<M 
%   Phi -- The SparseRandom matrix

%% Generate SparseRandom matrix   
    Phi = zeros(M,N);
    for ii = 1:N
        ColIdx = randperm(M);
        Phi(ColIdx(1:d),ii) = 1;
    end
end
