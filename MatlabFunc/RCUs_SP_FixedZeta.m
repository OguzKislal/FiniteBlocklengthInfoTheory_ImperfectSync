function [avg_error, s_val, debug_data_out] = RCUs_SP_FixedZeta(n,np,L, rho, constel, zetaFix, rateFix, g_list, g_hat_list, d_list, d_hat_list, tp, sigma_sq_list, nbrOfRealizations, SPMode, s_list,S_FAST_FLAG)
% Function sec3_getErrorProbability(n, rho, b, g_list ,ghat_list,sigma_sq_list, nbrOfRealizations, s_start))
% that computes the saddlepoiint approximation of the error probability given by the RCUs
if(nargin == 16)
   S_FAST_FLAG = 1 ; 
end
INDEPENDENT_CGF_MULT = 1; % will be removed. 
eps_debug = []; 
debug_data_out = []; 

eps_old = inf;
for ii = 1:length(s_list)
    s_candidate = s_list(ii);
    
    [epsilon, debug_data]= sec3_SaddlepointApprox(s_candidate, n,np,L, rho, constel, zetaFix, rateFix, g_list, g_hat_list, d_list, d_hat_list, tp, sigma_sq_list, nbrOfRealizations, SPMode,INDEPENDENT_CGF_MULT);
    eps_cur = mean(epsilon);
    eps_debug(ii)=eps_cur;
    if(S_FAST_FLAG == 1)
       if eps_cur < 1e-10
           debug_data_out = [] ; 
           debug_data_out = debug_data ; 
           break;
       end
       if eps_cur > eps_old
           eps_cur = eps_old;
           s_candidate = s_list(ii-1);        
           break;
       else
           debug_data_out = [] ; 
           debug_data_out = debug_data ; 
           eps_old = eps_cur;
       end
    end
end

if(S_FAST_FLAG == 1)
   avg_error = eps_cur;
else
   avg_error = eps_debug;
end
s_val = s_candidate;

end

function [epsilon, debug_data]= sec3_SaddlepointApprox(s, n,np,L, rho, constel, zetaFix, rateFix, g_list, g_hat_list, d_list, d_hat_list, tp, sigma_sq_list, nbrOfRealizations, SPMode,INDEPENDENT_CGF_MULT)
% that computes the quantities needed to compute the CGFs and its
% derivatives to then compute saddlepoiint approximation of the
% error probability given by the RCUs
nc = n./L ; 
nd = nc-np ; 
alpha = abs(d_list-d_hat_list)./tp; % in the article this is defined as alpha_not! 
if(SPMode == "RCUs")
   debug_data = [] ; 
   R_nats = rateFix.* log(2) ;
   for L_cnt = 1:L 
      i_s(:,:,L_cnt) = info_density(s,g_list(L_cnt,:),g_hat_list(L_cnt,:),alpha(L_cnt,:),sigma_sq_list(L_cnt,:),rho,constel,nd,nbrOfRealizations);
   end
   i_s = sum(i_s,3) ; 
   epsilon = mean(exp(-max(0, i_s - n*R_nats )),'all'); % Exact RCUs assuming 2^(n*rateFix) >> 1
elseif(SPMode == "Quasistatic")
   debug_data = [] ; 
   epsilon = nan(1, nbrOfRealizations);
%    options_Fzero = optimset('TolX',1e-4) ;
   for j = 1:nbrOfRealizations
       % Get channel, channel estimate, and effective noise:
       gTemp = g_list(:,j);
       ghatTemp = g_hat_list(:,j);
       alpha_sel = alpha(:,j) ; 
       g = gTemp .* (1-(alpha_sel >1)) ;
       ghat = ghatTemp .* (1-(alpha_sel >1)) ; 
       sigma_sq = sigma_sq_list(:,j);
       [epsilon(j)] = RCUsBound_SP_Markov(constel,rho,nd-1,nc,L,s,alpha_sel,g,ghat,rateFix) ;
   end
elseif(SPMode == "Normal")
   debug_data = [] ; 
   epsilon = nan(1, nbrOfRealizations);
   parfor j = 1:nbrOfRealizations
       % Get channel, channel estimate, and effective noise:
       gTemp = g_list(:,j);
       ghatTemp = g_hat_list(:,j);
       alpha_sel = alpha(:,j) ; 
       g = gTemp .* (1-(alpha_sel >1)) ; 
       ghat = ghatTemp .* (1-(alpha_sel >1)) ; 
       sigma_sq = sigma_sq_list(:,j);
       [epsilon(j)] = RCUsBound_Normal_Markov(constel,rho,nd-1,nc,L,s,alpha_sel,g,ghat,rateFix) ;
   end
end
end

function [i_s] = info_density(s,g,ghat,alpha,sigma_list,rho,constel,nd,N)

q_states = constel .* sqrt(rho) ; 
DATA_LOOPER = 5e1 ; 
% count = 0;
% i_s = zeros(nd,N) ;
for ii = 1 : DATA_LOOPER
   for trial = 1:N
%    parfor trial = 1:N
      if(alpha(trial) > 1)
         i_s(ii,trial) = 0 ; 
      else
         q = q_states(randi([1,length(constel)],1,nd+1));
         q_new = q(2:end) ; 
         q_old = q(1:end-1) ; 
         % [q_new.', q_old.']
         z = randcn(1,nd)  ; 
         x = (1-alpha(trial)).*q_new + alpha(trial).*q_old ; 
         y = g(trial).*x + z ;
         i_s(ii,trial) = sum(InfoDens_Calc(y,ghat(trial),q_new,constel,rho,s)) ;
      end
   end
end
% i_s = mean(i_s_temp,1) ; 
end

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

end

function [Eq_xl_log] = Expect_of_q_Xbar(Y_ld,H_l_est,constel,SNR_lin,s)

xbar_SNR = constel.*sqrt(SNR_lin) ; 
xbar_SNR = reshape(xbar_SNR,1,1,length(xbar_SNR)) ;  
temp = abs(repmat(Y_ld,size(xbar_SNR)) - xbar_SNR.*H_l_est).^2 ; 
mismatch_vals = exp(-temp.*s) ; % exp(-temp).^s
Eq_xl_log =  log(mean(mismatch_vals,3)); 

end

% function [s] = Evaluate_Opt_s(g,ghat,sigma_sq,rho,options_Fzero)
% 
% a = abs(g).^2.*rho + sigma_sq ; 
% b = rho.*abs(ghat).^2 ; 
% c = abs(g-ghat).^2.*rho+sigma_sq ; 
% 
% fun = @(s) sum( (a+b+s.*b)./((1+s.*b).^2)  - c) ;
% s = fzero(fun,[0,1.2],options_Fzero) ; 
% end

