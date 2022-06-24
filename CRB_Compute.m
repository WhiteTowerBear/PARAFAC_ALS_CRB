function [CRB,CRB_aa_part,CRB_bb_part]=CRB_Compute(M,K,N,P,Hs,Hr,Phi,var_noise) 
% This function aims at computing CRB of channel estimation 

%The references to theorems and equations refer to the following paper:
%
% L. Wei, C. Huang, G. C. Alexandropoulos, C. Yuen, Z. Zhang and M. Debbah, 
% "Channel Estimation for RIS-Empowered Multi-User MISO Wireless 
% Communications," in IEEE Transactions on Communications, vol. 69, 
% no. 6, pp. 4144-4157, June 2021.

%License: If you in any way use this code for research that results in 
% publications, please cite our original article listed above.

KR_Hs_Phi=kr(Hs.',Phi);
temp1=KR_Hs_Phi'*KR_Hs_Phi/var_noise;
KR_Phi_Hr=kr(Phi,Hr);
temp2=KR_Phi_Hr'*KR_Phi_Hr/var_noise;
 

% Compute the Psi_aa
Psi_aa=kron(eye(K),temp1);

% Compute the Psi_bb
Psi_bb=kron(eye(M-1),temp2);


Psi_ab=zeros(K*N,(M-1)*N);
     for m=2:M
               for k=1:K
                   noise_temp=zeros(P*M,K*P);
                    for p=1:P
                        noise_temp((m-1)*P+p,(p-1)*K+k)=1;
                    end  
                    temp_ab_pro=KR_Hs_Phi'*noise_temp*KR_Phi_Hr;     
                    Psi_ab((k-1)*N+1:k*N,(m-2)*N+1:((m-1)*N))=temp_ab_pro/var_noise;
                end
     end


CRB_aa=inv(Psi_aa-Psi_ab*inv(Psi_bb)*Psi_ab');                
CRB_bb=inv(Psi_bb-Psi_ab'*inv(Psi_aa)*Psi_ab);      
CRB_aa_part=trace(CRB_aa);
CRB_bb_part=trace(CRB_bb);
CRB_part=CRB_aa_part+CRB_bb_part;
CRB_aa_part=CRB_aa_part/(K*N);
CRB_bb_part=CRB_bb_part/((M-1)*N);


CRB=CRB_part/((K+M-1)*N);

end