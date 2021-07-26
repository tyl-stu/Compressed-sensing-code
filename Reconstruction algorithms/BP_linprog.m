function [ alpha ] = BP_linprog( s,Phi )
%BP_linprog(Basis Pursuit with linprog) Summary of this function goes here
%Version: 1.0 written by jbb0523 @2016-07-21 
%Reference:Chen S S, Donoho D L, Saunders M A. Atomic decomposition by
%basis pursuit[J]. SIAM review, 2001, 43(1): 129-159.(Available at: 
%http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.37.4272&rep=rep1&type=pdf)
%   Detailed explanation goes here
%   s = Phi * alpha (alpha is a sparse vector)  
%   Given s & Phi, try to derive alpha
    [s_rows,s_columns] = size(s);  
    if s_rows<s_columns  
        s = s';%s should be a column vector  
    end 
    p = size(Phi,2);
    %according to section 3.1 of the reference
    c = ones(2*p,1);
    A = [Phi,-Phi];
    b = s;
    lb = zeros(2*p,1);
    x0 = linprog(c,[],[],A,b,lb);
    alpha = x0(1:p) - x0(p+1:2*p);
end
