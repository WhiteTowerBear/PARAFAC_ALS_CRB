function [SumRate_MRC,SumRate_ZF]=SumRate(M,K,N,P,Hs,Hr,Phi,var_noise,p_u,var_channel) 
% This function aims at computing MRC
for k=1:K
    R_MRC(k,1)=log2(1+p_u*(M-1)^2*var_channel^2/ (p_u*(M-1)*var_channel*(K-1)*var_channel) +var_noise);
end
SumRate_MRC=sum(R_MRC)/K;
SumRate_ZF=log2(1+p_u/va_noise);

end