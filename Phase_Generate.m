function [Phi]=Phase_Generate(P,N)
% This function aims at generating random pilot phase matrix
% Phi=exp(1i*2*pi*rand(P,N));% randomly generation
% 
% X_temp=randn(P,P);
% [U,~,~]=svd(X_temp);
% Phi=U(1:P,1:N);

TEMP=dftmtx(N);
Phi=1/sqrt(N)*TEMP(1:P,1:N);