function [g_list,d_list,G_hat_list,D_hat_list,G_PerfSync_hat_list,avg_error_h,avg_error_d] = DataGenerator_ImperfectSync(tp,dmax,rho,N,nbrOfRealizations,num_seqs,L_divBranch)

if(nargin == 6)
   L_divBranch = 1; 
end
ts = tp/N;
% nc = n ./ L_divBranch

g_list = randcn(1, nbrOfRealizations) ;
d_list = dmax*rand(1,nbrOfRealizations./L_divBranch); % Delay truths drawn from Unif(0,3)
d_list = repmat(d_list,L_divBranch,1) ; 
d_list = d_list(:).' ; 
G_hat_list = zeros(length(num_seqs),nbrOfRealizations);
D_hat_list = zeros(length(num_seqs),nbrOfRealizations);

N0 = 1;

h_err_list = zeros(length(num_seqs),nbrOfRealizations) ;
d_err_list = zeros(length(num_seqs),nbrOfRealizations) ;

for numm_cnt = 1:length(num_seqs)
    
    ns = num_seqs(numm_cnt);
    
    xm_p_unscaled  = mseq(2,ns,1,2);
    xm_p = sqrt(rho)*xm_p_unscaled;
    np = length(xm_p);
    
    disp(['Count number: ' num2str(numm_cnt)]);
    
    for i = 1:nbrOfRealizations
        h_true = g_list(i);
        d_true = d_list(i);
        
        % Generate X and e
        [X_real,e_real] = gen_X_e_data(d_true,ts,tp,xm_p,dmax,rho);
        
        % Generate y
        y_noiseless = h_true*X_real*e_real;
        noise = sqrt(N0/ts)*randcn(length(y_noiseless),1);
        y = y_noiseless + noise;
        %            y = y_noiseless ;
        
        % Synchronization algorithm
        
        d_set = [];
        q_max = dmax/ts;
        q_set = 0:1:round(q_max);
        
        for q = q_set
            d_set = [d_set q*ts];
            [X,~] = gen_X_e(q*ts,ts,tp,xm_p,dmax);
            
            %solve for eps_roots;
            
            % Find polynomial coefficients for N-
            % N(eps) = y3 + y2*eps + y1*eps^2
            A = (X'*y)*(y'*X);
            a = A(1);
            b = A(2);
            c = A(3);
            d = A(4);
            
            y1 = a./(ts^2) - (b+c)./(ts^2) + d/(ts^2);
            y2 = -2*a./ts + (b+c)./ts;
            y3 = a;
            
            % Find polynomial coefficients for D-
            % D(eps) = x3 + x2*eps + x1*eps^2
            
            A = X'*X;
            a = A(1);
            b = A(2);
            c = A(3);
            d = A(4);
            
            x1 = a./(ts^2) - (b+c)./(ts^2) + d/(ts^2);
            x2 = -2*a./ts + (b+c)./ts;
            x3 = a;
            
            
            % coefficients of polynomial to find roots of-
            % z1*eps^2 + z2*eps + z3
            z3 = x3*y2-y3*x2;
            z2 = 2*x3*y1-2*y3*x1;
            z1 = x2*y1-x1*y2;
            eps_roots = real(roots([z1 z2 z3]));
            
            for l = [1,2]
                eps = eps_roots(l);
                if eps>0 && eps<=ts
                    d_set = [d_set q*ts+eps];
                end
            end
        end
        % temporaryFix for floating point problem in d_set
        d_set = d_set + 1e-10 ;
        d_set(d_set>dmax) = dmax ;
        
        %
        % find objective function values associated with each possible d
        
        obj = zeros(1,length(d_set));
        
        for k = 1:length(d_set)
            delayCheck = d_set(k);
            [X,e] = gen_X_e(delayCheck,ts,tp,xm_p,dmax);
            mu = X*e;
            obj(k) = -(abs(mu'*y).^2)/(norm(mu)^2);
        end
        
        [~,ind] = min(obj);
        d_hat = d_set(ind);
        % Channel estimation
        
        [X_hat,e_hat] = gen_X_e(d_hat,ts,tp,xm_p,dmax);
        mu_hat = X_hat*e_hat;
        h_hat = (mu_hat'*y)./(mu_hat'*mu_hat);

        mu_PerfSync = X_real*e_real;
        h_PerfSync_hat = (mu_PerfSync'*y)./(mu_PerfSync'*mu_PerfSync);
        
        h_err = abs(h_true-h_hat).^2;
        d_err = abs(d_true-d_hat).^2./tp;
        h_PerfSync_err = abs(h_true-h_PerfSync_hat).^2;

        h_err_list(numm_cnt,i) = h_err ;
        h_PerfSync_err_list(numm_cnt,i) = h_PerfSync_err ;
        d_err_list(numm_cnt,i) = d_err ;
        
        G_hat_list(numm_cnt,i) = h_hat;
        G_PerfSync_hat_list(numm_cnt,i) = h_PerfSync_hat ; 
        D_hat_list(numm_cnt,i) = d_hat;
    end
    
end

avg_error_h = mean(h_err_list,2) ;
avg_error_d = mean(d_err_list,2);


end

