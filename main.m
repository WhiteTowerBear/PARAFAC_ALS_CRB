%========================================================a
%  Author: Wei Li;
%  Date:   11/17/ 2019,  SUTD;
%  Version: V2.0 
%  Note: This is mainly for computing CRB of PARAFAC ALS based channel 
%        estimation algorithm in RIS assisted systems. 

%The references to theorems and equations refer to the following paper:
%
% L. Wei, C. Huang, G. C. Alexandropoulos, C. Yuen, Z. Zhang and M. Debbah, 
% "Channel Estimation for RIS-Empowered Multi-User MISO Wireless 
% Communications," in IEEE Transactions on Communications, vol. 69, 
% no. 6, pp. 4144-4157, June 2021.

%License: If you in any way use this code for research that results in 
% publications, please cite our original article listed above.
%%========================================================
clc
clear all 
fid=fopen('result.txt','a+');
%%========================================================
%  Set the parameters
M=16; % The number of antennas at BS
K=16; % The number of users
T=16; % The numebr of time slots
N=16; % The number of elements placed on RIS
P=16;  % The number of pilot phase matrices
var_channel=1;    % The variance of reflecting channel 

FRAME=50;
iter=30;
diff_temp=[];
snr=1;


for SNR=0:5:10  
    SNR
    mc_sum_Hs=0;
    mc_sum_Hr=0;
    mmse_sum_Hs=0;
    mmse_sum_Hr=0;
    CRB_total=0;
    CRB_aa_total=0;
    CRB_bb_total=0;
   
    for frame=1:1:FRAME
    % Generate transmitted signals X
    [X,X_inv]=Transceiver(M,T); 
    var_noise=10^(-0.1*SNR);
    %%========================================================
    %  Pass channel
    Hr= sqrt(var_channel/2)*(randn(K,N)+1i*randn(K,N)); 
    Hs=sqrt(var_channel/2)*(randn(N,M)+1i*randn(N,M));
    Hs(:,1)=1;
    % Generate phase matrix Phi
    [Phi]=Phase_Generate(P,N);
    % Received signals
    noise=zeros(K,T,P);
    rec_y=zeros(K,T,P);
    for p=1:P
        noise(:,:,p)=sqrt(var_noise/2)*(randn(K,T)+1i*randn(K,T));
        rec_y(:,:,p)=Hr*diag(Phi(p,:))*Hs*X+noise(:,:,p);
        rec_y_TEMP(:,:,p)=rec_y(:,:,p)*X_inv;
    end
%     %%========================================================
%     % Receiver
%     % Mode-1, mode-2, mode-3 unfolded forms of PARAFAC    
    for m=1:M
        for p=1:P
            for k=1:K
                Z_KP_M((p-1)*K+k,m)=rec_y_TEMP(k,m,p); % MODE2
                Z_PM_K((m-1)*P+p,k)=rec_y_TEMP(k,m,p); % MODE1
            end
        end
    end
  %%  
    % Initialization
    Hs_est=zeros(N,M,iter+1);
    Hr_est=zeros(K,N,iter+1);
    Hr_est(:,:,1)=sqrt(var_channel/2)*(randn(K,N)+1i*randn(K,N));
    
    A1=zeros(P*M,N,iter);
    A1_inv=zeros(N,P*M,iter);
    A2=zeros(K*P,N,iter);
    A2_inv=zeros(N,K*P,iter);
    
    %%========================================================
    %% PARAFAC 
    % Iteration between channels
    com_Hs=ones(iter+1,1);
    for i=2:iter        
        A2(:,:,i-1)=kr(Phi,Hr_est(:,:,i-1));
        A2_inv(:,:,i-1)=inv(A2(:,:,i-1)'*A2(:,:,i-1))*A2(:,:,i-1)';
        Hs_est(:,:,i)=A2_inv(:,:,i-1)*Z_KP_M;
        for n=1:N % for scaling ambiguity removal later:
            Hs_est(n,:,i) = Hs_est(n,:,i) / Hs_est(n,1,i);
        end
        A1(:,:,i)=kr(Hs_est(:,:,i).',Phi);
        A1_inv(:,:,i)=inv(A1(:,:,i)'*A1(:,:,i))*A1(:,:,i)';
        Hr_est(:,:,i)=(A1_inv(:,:,i)*Z_PM_K).';
        fit=0;
        for p=1:P
            fit=fit+norm(Z_KP_M-kr(Phi,Hr_est(:,:,i))*Hs_est(:,:,i),'fro')^2;
        end
        com_Hs(i,1)=fit;
        delta(i,1)=(com_Hs(i-1,1)-com_Hs(i,1))/com_Hs(i,1);
        if abs(delta(i,1))<1e-5
            break
        end
    end   
     %% LS estimation
%         A21=kr(Phi,Hr);
%         A2_inv1=inv(A21'*A21)*A21';
%         Hs_ls=A2_inv1*Z_KP_M;
%         A11=kr(Hs.',Phi);
%         A11=inv(A11'*A11)*A11';
%         Hr_ls=(A11*Z_PM_K).';
%     %%=================================================
%%    CRB part: CRB_aa_part is the CRB of Hr, and CRB_bb_part is the CRB of Hs
       [CRB_iter,CRB_aa_part_iter,CRB_bb_part_iter]=CRB_Compute(M,K,N,P,Hs,Hr,Phi,var_noise);
       CRB_aa_total=CRB_aa_total+real(CRB_aa_part_iter);
       CRB_bb_total=CRB_bb_total+real(CRB_bb_part_iter);
       CRB_total=CRB_total+real(CRB_iter);
%%========================================================%
%%
% % Simulation result
      mc_sum_Hs=mc_sum_Hs+norm(Hs(:,2:M)-Hs_est(:,2:M,i),'fro')^2/(norm(Hs(:,2:M),'fro')^2);
      mc_sum_Hr=mc_sum_Hr+norm(Hr-Hr_est(:,:,i),'fro')^2/(norm(Hr,'fro')^2);
%       mmse_sum_Hr=mmse_sum_Hr+norm(Hr-Hr_ls,'fro')^2/(norm(Hr,'fro')^2);   
%       mmse_sum_Hs=mmse_sum_Hs+norm(Hs-Hs_ls,'fro')^2/(norm(Hs,'fro')^2); 
    end
    diff_Hs(1,snr)=mc_sum_Hs/FRAME;
    diff_Hr(1,snr)=mc_sum_Hr/FRAME;
%     ls_Hr(1,snr)=mmse_sum_Hr/FRAME;
%     ls_Hs(1,snr)=mmse_sum_Hs/FRAME;
    CRB(1,snr)=CRB_total/FRAME;
    CRB_aa(1,snr)=CRB_aa_total/FRAME;
    CRB_bb(1,snr)=CRB_bb_total/FRAME;
    
    snr=snr+1;
end
fprintf(fid,'\n \n M=%d  K=%d  N=%d  T=%d  P=%d  %d Monte Carlo simulations \n',M,K,N,T,P,FRAME);
fprintf(fid,' %s\r\n ',datestr(now,31));
fprintf(fid,'NMSE of Hr is \n');
fprintf(fid, ' %d  ', diff_Hr);
fprintf(fid,'\n NMSE of Hs is \n');
fprintf(fid, ' %d  ', diff_Hs);
% fprintf(fid,'\n LS of Hr is \n');
% fprintf(fid, ' %d  ', ls_Hr);
% fprintf(fid,'\n LS of Hs is \n');
% fprintf(fid, ' %d  ', ls_Hs);
% fprintf(fid,'\n CRB is \n');
% fprintf(fid, ' %d  ', CRB);
% fprintf(fid,'\n CRB_Hr is \n');
% fprintf(fid, ' %d  ', CRB_aa);
% fprintf(fid,'\n CRB_Hs is \n');
% fprintf(fid, ' %d  ', CRB_bb);
fclose(fid);
figure
semilogy(0:5:10,diff_Hr,'-bo','linewidth',2)
hold on
semilogy(0:5:10,diff_Hs,'-ys','linewidth',2)
hold on
% semilogy(0:5:30,ls_Hr,'-g+','linewidth',2)
% hold on
% semilogy(0:5:30,ls_Hs,'-kp','linewidth',2)
% % hold on
semilogy(0:5:10,CRB_aa,'-.bo','linewidth',1)
hold on 
semilogy(0:5:10,CRB_bb,'-.ys','linewidth',1)

legend({'Hr (CRB)','Hs (CRB)'});
legend({'Hr','Hs','CRB (Hr)', 'CRB (Hs)'});
title(['M=',num2str(M),',K=',num2str(K),',T=',num2str(T), ',N=',num2str(N),',P=',num2str(P)]);
xlabel('SNR (dB)');
ylabel('NMSE');
grid on 
