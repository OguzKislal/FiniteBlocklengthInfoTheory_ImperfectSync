clc
clear all; 

disp('Sim Starts')

sim_number = 1; % Change this variable for the different parameter settings.


n = 288 ; % Blocklength (nb *nc in paper)
np = 0; % Number of pilot symbols to deal with pilot overhead, leave it at 0 if you don't want to deal with that
rateFix = 0.104 ; % rate bit per channel use

% SPMode = "Quasistatic" ;
sList = [0.1:0.1:0.9] ; % Parameter s, code searches the optimal s in the list.
tp = 1; % pulse shape period
dmax = 12; % maximum value of delay 
N = 10; % Oversampling rate
% constel = [1+1i, -1-1i] ; 
constel = [1+1i, -1+1i, -1-1i, 1-1i] ; % Constellation
constel = constel ./ sqrt(mean(abs(constel).^2)) ; % Power normalization

% nc = n ./ L_divBranch
JOINT_SYNCR_FLAG = 0 % Joint Synchronization
IMPERFECT_SYNC_CH_SIM = 1 % Imperfect sync. and channel estimatiopn
PERFECT_SYNC_SIM = 0  % Perfect Sync. 
PERFECT_SYNC_CH_SIM = 0  % Perfect Sync. and channel estimation. 

% 
SNR_SEARCH = 1 
NP_FAST_SIM = 1 ; % When set n_p search is stopped when no further improvement cannot be obtained. 
NP_FAST_SIM = SNR_SEARCH | NP_FAST_SIM  %
TARGET_ERROR_RATE = 1e-5; % Only useful when SNR_SEARCH is set, target error rate

% g_list = randcn(L_divBranch, nbrOfRealizations) ; 
% d_list = 3*rand(1,nbrOfRealizations); % Delay truths drawn from Unif(0,3)
if(sim_number <= 6)
   SPMode = "Quasistatic"  % Can be "RCUs", "Quasistatic"
elseif(sim_number <= 12 )
   SPMode = "RCUs" % Can be "RCUs", "Quasistatic"
   sim_number = sim_number - 6 ;
else
   SPMode = "Normal"  % Can be "RCUs", "Quasistatic"
   sim_number = sim_number - 12 ;
end

%Depending on the figure comment/comment out of the following parts: 

%%
nbrOfRealizations = 5e6 ;
JOINT_SYNCR_FLAG = 0  
if(sim_number == 1)
   L_divBranch = 1 ; 
   rho_db = [35 45]; 
elseif(sim_number == 2)
   L_divBranch = 2 ; 
   rho_db = [11 21]; 
elseif(sim_number == 3)
   L_divBranch = 4 ; 
   rho_db = [5 15]; 
elseif(sim_number == 4)
   L_divBranch = 8 ; 
   rho_db = [0 10]; 
elseif(sim_number == 5)
   L_divBranch = 12 ; 
   rho_db = [0 10]; 
elseif(sim_number == 6)
   L_divBranch = 24 ; 
   rho_db = [0 10]; 
else
   L_divBranch = 4 ;
   rho_db = [5 10] ; 
   nbrOfRealizations = 1e2 ; 
   SPMode = "Quasistatic" ; % Can be "RCUs", "Quasistatic", "Normal"
end

%%
%JOINT_SYNCR_FLAG = 1  
%if(sim_number == 1)
%   L_divBranch = 1 ; 
%   rho_db = [37 46]; 
%elseif(sim_number == 2)
%   L_divBranch = 2 ; 
%   rho_db = [15 25]; 
%elseif(sim_number == 3)
%   L_divBranch = 4 ; 
%   rho_db = [6 14]; 
%elseif(sim_number == 4)
%   L_divBranch = 8 ; 
%   rho_db = [1 9]; 
%elseif(sim_number == 5)
%   L_divBranch = 12 ; 
%   rho_db = [-1 9]; 
%elseif(sim_number == 6)
%   L_divBranch = 24 ; 
%   rho_db = [-1 9]; 
%else
%   L_divBranch = 1 ;
%   rho_db = [10 15] ; 
%   nbrOfRealizations = 1e2 ; 
%   SPMode = "RCUs" ; % Can be "RCUs", "Quasistatic", "Normal"
%end

%%
%PERFECT_SYNC_SIM = 1  
%if(sim_number == 1)
%   L_divBranch = 1 ; 
%   rho_db = [37 47]; 
%elseif(sim_number == 2)
%   L_divBranch = 2 ; 
%   rho_db = [13 23]; 
%elseif(sim_number == 3)
%   L_divBranch = 4 ; 
%   rho_db = [3 13]; 
%elseif(sim_number == 4)
%   L_divBranch = 8 ; 
%   rho_db = [0 10]; 
%elseif(sim_number == 5)
%   L_divBranch = 12 ; 
%   rho_db = [0 10]; 
%elseif(sim_number == 6)
%   L_divBranch = 24 ; 
%   rho_db = [0 10]; 
%else
%   L_divBranch = 1 ;
%   rho_db = [10 15] ; 
%   nbrOfRealizations = 1e2 ; 
%   SPMode = "RCUs" ; % Can be "RCUs", "Quasistatic", "Normal"
%end

%%
%PERFECT_SYNC_CH_SIM = 1  
%if(sim_number == 1)
%   L_divBranch = 1 ; 
%   rho_db = [30 44]; 
%elseif(sim_number == 2)
%   L_divBranch = 2 ; 
%   rho_db = [10 20]; 
%elseif(sim_number == 3)
%   L_divBranch = 4 ; 
%   rho_db = [0 11]; 
%elseif(sim_number == 4)
%   L_divBranch = 8 ; 
%   rho_db = [-3 7]; 
%elseif(sim_number == 5)
%   L_divBranch = 12 ; 
%   rho_db = [-5 5]; 
%elseif(sim_number == 6)
%   L_divBranch = 24 ; 
%   rho_db = [-5 5]; 
%else
%   L_divBranch = 1 ;
%   rho_db = [10 15] ; 
%   nbrOfRealizations = 1e2 ; 
%   SPMode = "RCUs" ; % Can be "RCUs", "Quasistatic", "Normal"
%end

if(SNR_SEARCH == 1 && length(rho_db) ~= 2)
   error('USER ERROR: Wrong SNR vector length!')
end
SNR_SEARCH_COUNT_MAX = max( ceil(log2(diff(rho_db)/ 0.05))*SNR_SEARCH,1);  % 0.05 SNR interval guaranteed


num_seqs_vec = 2:floor(log2(n/L_divBranch));
if(PERFECT_SYNC_CH_SIM == 1)
   num_seqs_vec = [] ; 
elseif(L_divBranch <= 8)
   num_seqs_vec(end) = [] ; 
end

nbrOfRealizations = L_divBranch.*round(nbrOfRealizations./L_divBranch) ; 
nbrOfRealizations = nbrOfRealizations./L_divBranch ; 

num_seqs_vec
rho_db
N 

np_vec = [] ; 
for searchCnt = 1 : SNR_SEARCH_COUNT_MAX
   rho = 10^(mean(rho_db)/10) ;
   avg_errors_perfectDelay = inf(1,length(num_seqs_vec));
   avg_errors_imperfect = inf(1,length(num_seqs_vec));
   for k = 1:length(num_seqs_vec)
      num_seqs = num_seqs_vec(k) ; 
      %Generate g_list and d_list, along with G_hat_list and D_hat_list.
      if(JOINT_SYNCR_FLAG == 1)
         [g_list,d_list,G_hat_list,D_hat_list,G_PerfSync_hat_list,avg_error_h,avg_error_d] = DataGenerator_ImperfectSync_joint(tp,dmax,rho,N,nbrOfRealizations,num_seqs,L_divBranch) ; 
         for ii = 1 : L_divBranch
            G_hat_list_new(ii,:,:) = G_hat_list(:,:,ii).' ; 
            D_hat_list_new(ii,:,:) = D_hat_list.'  ;  
            G_PerfSync_hat_list_new(ii,:,:) = G_PerfSync_hat_list(:,:,ii).' ; 
         end
         d_list = repmat(d_list,L_divBranch,1) ; 
      else
         [g_list,d_list,G_hat_list,D_hat_list,G_PerfSync_hat_list,avg_error_h,avg_error_d] = DataGenerator_ImperfectSync(tp,dmax,rho,N,nbrOfRealizations.*L_divBranch,num_seqs,L_divBranch) ; 
         g_list = reshape(g_list,L_divBranch,nbrOfRealizations) ; 
         d_list = reshape(d_list,L_divBranch,nbrOfRealizations) ; 
         for ii = 1 :size(G_hat_list,1) 
            G_hat_list_new(:,:,ii) = reshape(G_hat_list(ii,:),L_divBranch,nbrOfRealizations) ; 
            D_hat_list_new(:,:,ii) = reshape(D_hat_list(ii,:),L_divBranch,nbrOfRealizations) ; 
            G_PerfSync_hat_list_new(:,:,ii) = reshape(G_PerfSync_hat_list(ii,:),L_divBranch,nbrOfRealizations) ; 
         end
         %nbrOfRealizations = nbrOfRealizations./L_divBranch ; 
      end
      G_hat_list = G_hat_list_new ; 
      D_hat_list = D_hat_list_new ; 
      G_PerfSync_hat_list = G_PerfSync_hat_list_new ; 
      clear G_hat_list_new D_hat_list_new G_PerfSync_hat_list_new
      %[g_hat_list, d_hat_list] = synch_channel_est(g_list,d_list,np_synch,tp); % ML estimates using pilot symbols
      sigma_sq_list = ones(size(g_list)) ; % The noise has unit power.
      tic
     
   %% Get RCUs with different number of pilot symbols, and imperfect synchronization and channel estimation
   
      np = 2^(num_seqs)-1; % number of pilot symbols
      np_vec(k) = np ; 
      disp(['  Number of pilots: ' num2str(np)]);
        
      if(PERFECT_SYNC_SIM == 1)
         % Compute RCUs in case of perfect delay
         if(NP_FAST_SIM == 1 && k>2 && avg_errors_perfectDelay(k-1)*0.95 > avg_errors_perfectDelay(k-2) )
            break ;
         else
            disp('Computing RCUs with perfect delay information');
            [avg_error_perfectDelay, s_val,~] = RCUs_SP_FixedZeta(n,np,L_divBranch, ...
            rho, constel, NaN, rateFix, g_list, G_PerfSync_hat_list, d_list, d_list, tp,...
            sigma_sq_list, nbrOfRealizations,SPMode,sList);
            avg_errors_perfectDelay(k) = avg_error_perfectDelay ;
         end
      end
   
      % Compute RCUs in case of perfect Channel
      if(IMPERFECT_SYNC_CH_SIM == 1)
         if(NP_FAST_SIM == 1 && k>2 && avg_errors_imperfect(k-1)*0.95 > avg_errors_imperfect(k-2))
            break ;
         else
            disp('Computing RCUs with imperfect delay and channel information');
            [avg_error_imperfect, s_val,~] = RCUs_SP_FixedZeta(n,np,L_divBranch, ...
                 rho, constel, NaN, rateFix, g_list, G_hat_list, d_list, D_hat_list, tp,...
                 sigma_sq_list, nbrOfRealizations,SPMode,sList);
            avg_errors_imperfect(k) = avg_error_imperfect;
         end
      end
   end
   %% Compute RCUs in case of perfect delay and channel information
   avg_error_perfectDelayChannel = NaN ; 
   if(PERFECT_SYNC_CH_SIM == 1)
      disp('Computing RCUs with perfect delay and channel information');
      [g_list,d_list,G_hat_list,D_hat_list,G_PerfSync_hat_list,avg_error_h,avg_error_d] = DataGenerator_ImperfectSync_joint(tp,dmax,rho,N,nbrOfRealizations,2,L_divBranch) ; 
      d_list = repmat(d_list,L_divBranch,1) ; 
      sigma_sq_list = ones(size(g_list)) ; % The noise has unit power.
      [avg_error_perfectDelayChannel, s_val,~] = RCUs_SP_FixedZeta(n,np,L_divBranch, ...
         rho, constel, NaN, rateFix, g_list, g_list, d_list, d_list, tp,...
         sigma_sq_list, nbrOfRealizations,SPMode,sList);
   end

   if(SNR_SEARCH == 1)
      if(PERFECT_SYNC_SIM == 1)
         temp = min(avg_errors_perfectDelay) ; 
      elseif(PERFECT_SYNC_CH_SIM == 1)
         temp = avg_error_perfectDelayChannel ; 
      else
         temp = min(avg_errors_imperfect) ; 
      end
      if( 100*abs(temp- TARGET_ERROR_RATE)/TARGET_ERROR_RATE <= 0.1 )
         break; 
      elseif(temp < TARGET_ERROR_RATE)
         rho_db  = [rho_db(1), mean(rho_db)] ; 
      elseif(temp > TARGET_ERROR_RATE)
         rho_db = [mean(rho_db), rho_db(2)] ; 
      end
   end
   rho_db
   
end

save_string = ['MarkovChain_'] ; 

if(PERFECT_SYNC_CH_SIM == 1)
   save_string = [save_string, 'PerfectChSyncResult']
elseif(PERFECT_SYNC_SIM == 1)
   save_string = [save_string, 'PerfectSyncResult']
elseif(JOINT_SYNCR_FLAG == 1)
   save_string = [save_string, 'JointSyncResults_'] ; 
else
   save_string = [save_string, 'IIDSyncResults_'] ; 
end

save_string = [save_string,convertStringsToChars(SPMode),'_L',num2str(L_divBranch),'_N',num2str(N)] ;
if(SNR_SEARCH == 1)
   save_string = [save_string, '_SNR_Search'] ; 
else
   save_string = [save_string, '_SNR',num2str(rho_db)]; 
end

save_string = [save_string,'.mat']; 
% save(save_string,'nbrOfRealizations','avg_errors_perfectDelay','avg_errors_imperfect','avg_error_perfectDelayChannel','np_vec','rho_db','N','L_divBranch')
% toc
disp('Simulation ends')


% if(sim_number == 1)
%    L_divBranch = 2 ; 
%    rho_db = 0; 
%    nbrOfRealizations = 1e4;
% elseif(sim_number == 2)
%    L_divBranch = 2 ; 
%    rho_db = 5; 
%    nbrOfRealizations = 2e4;
% elseif(sim_number == 3)
%    L_divBranch = 2 ; 
%    rho_db = 10;
%    nbrOfRealizations = 1e5; 
% elseif(sim_number == 4)
%    L_divBranch = 2 ; 
%    rho_db = 15;
%    nbrOfRealizations = 1e6;
% elseif(sim_number == 5)
%    L_divBranch = 2 ; 
%    rho_db = 20;
%    nbrOfRealizations = 4e6;
% elseif(sim_number == 6)
%    L_divBranch = 2 ; 
%    rho_db = 25; 
%    nbrOfRealizations = 5e6;
% elseif(sim_number == 7)
%    L_divBranch = 4 ; 
%    rho_db = -4; 
%    nbrOfRealizations = 1e4;
% elseif(sim_number == 8)
%    L_divBranch = 4; 
%    rho_db = 0;
%    nbrOfRealizations = 1e4;
% elseif(sim_number == 9)
%    L_divBranch = 4 ; 
%    rho_db = 5;
%    nbrOfRealizations = 3e5;
% elseif(sim_number == 10)
%    L_divBranch = 4 ; 
%    rho_db = 8;
%    nbrOfRealizations = 3e6;
% elseif(sim_number == 11)
%    L_divBranch = 4 ; 
%    rho_db = 11; 
%    nbrOfRealizations = 5e6;
% elseif(sim_number == 12)
%    L_divBranch = 4 ; 
%    rho_db = 14; 
%    nbrOfRealizations = 5e6;
% elseif(sim_number == 13)
%    L_divBranch = 8 ; 
%    rho_db = 11;
%    nbrOfRealizations = 1e6;
% elseif(sim_number == 14)
%    L_divBranch = 8 ; 
%    rho_db = -5;
%    nbrOfRealizations = 1e4;
% elseif(sim_number == 15)
%    L_divBranch = 8 ; 
%    rho_db = -1;
%    nbrOfRealizations = 2e4;
% elseif(sim_number == 16)
%    L_divBranch = 8; 
%    rho_db = 2; 
%    nbrOfRealizations = 3e5;
% elseif(sim_number == 17)
%    L_divBranch = 8; 
%    rho_db =  4; 
%    nbrOfRealizations = 3e6;
% elseif(sim_number == 18)
%    L_divBranch = 8 ; 
%    rho_db = 5;
%    nbrOfRealizations = 5e6;
% elseif(sim_number == 19)
%    L_divBranch = 12; 
%    rho_db = -8;
%    nbrOfRealizations = 1e3;
% elseif(sim_number == 20)
%    L_divBranch = 12; 
%    rho_db = -5;
%    nbrOfRealizations = 1e4;
% elseif(sim_number == 21)
%    L_divBranch = 12; 
%    rho_db = -2; 
%    nbrOfRealizations = 1e4;
% elseif(sim_number == 22)
%    L_divBranch = 12; 
%    rho_db = 1; 
%    nbrOfRealizations = 1e5;
% elseif(sim_number == 23)
%    L_divBranch = 12; 
%    rho_db = 2;
%    nbrOfRealizations = 4e5;
% elseif(sim_number == 24)
%    L_divBranch = 12; 
%    rho_db = 3;
%    nbrOfRealizations = 2e6;
% elseif(sim_number == 25)
%    L_divBranch = 12; 
%    rho_db = 5;
%    nbrOfRealizations = 5e6;
% else
%    L_divBranch = 1 ;
%    rho_db = 10;
%    nbrOfRealizations = 1e3 ; 
%    SPMode = "Normal" ; % Can be "RCUs", "Quasistatic", "Normal"
% end
