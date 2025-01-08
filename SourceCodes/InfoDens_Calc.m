function [InfoDens] = InfoDens_Calc(y,h_est,x,constel,rho,s)
H_est_sqr = abs(h_est).^2 ; 
% s = s_def() ; 
%% Mismatch Decoding metric calculation
if(length(constel) == 1) % if Gaussian Codebook
   InfoDens = -s .* abs(y - h_est.*x).^2 + s.*abs(y).^2./(1+s.*rho.*H_est_sqr) + log(1+s.*rho.*H_est_sqr) ; 
else
   q_xl_log = -s.*abs(y - h_est.*x).^2 ;  
   Eq_xl_log = Expect_of_q_Xbar(y,h_est,constel,rho,s) ; 
   InfoDens = q_xl_log - Eq_xl_log ;
end


function [Eq_xl_log] = Expect_of_q_Xbar(Y_ld,H_l_est,constel,SNR_lin,s)

xbar_SNR = constel.*sqrt(SNR_lin) ; 
xbar_SNR = reshape(xbar_SNR,1,1,length(xbar_SNR)) ;  
temp = abs(repmat(Y_ld,size(xbar_SNR)) - xbar_SNR.*H_l_est).^2 ; 
mismatch_vals = exp(-temp.*s) ; % exp(-temp).^s
Eq_xl_log =  log(mean(mismatch_vals,3)); 