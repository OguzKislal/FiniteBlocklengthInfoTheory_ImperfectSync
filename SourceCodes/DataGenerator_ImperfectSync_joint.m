function [g_list,d_list,G_hat_list,D_hat_list,G_PerfSync_hat_list,avg_error_h,avg_error_d] = DataGenerator_ImperfectSync_joint(tp,dmax,rho,N,nbrOfRealizations,num_seqs,L_divBranch)

% Generating imperfectly synchronized
if(nargin == 6)
    L_divBranch = 1;
end
ts = tp/N;
% nc = n ./ L_divBranch


L = L_divBranch;

% g_list = zeros(L,nbrOfRealizations);
g_list = randcn(L, nbrOfRealizations);

d_list = dmax*rand(1,nbrOfRealizations); % Delay truths drawn from Unif(0,3)
% d_list = repmat(d_list,L_divBranch,1) ;
% d_list = d_list(:).' ;
G_hat_list = zeros(length(num_seqs),nbrOfRealizations,L);
D_hat_list = zeros(length(num_seqs),nbrOfRealizations);

G_PerfSync_hat_list = zeros(length(num_seqs),nbrOfRealizations,L);
N0 = 1;

h_err_list = zeros(length(num_seqs),nbrOfRealizations,L) ;
d_err_list = zeros(length(num_seqs),nbrOfRealizations) ;

%h_PerfSync_err_list = zeros(length(num_seqs),nbrOfRealizations,L);
for numm_cnt = 1:length(num_seqs)
    
    ns = num_seqs(numm_cnt);
    
    xm_p_unscaled  = mseq(2,ns,1,2);
%     xm_p_unscaled  = ones(2^ns-1,1)  ; 
%     xm_p_unscaled(2:2:end) = -1 ; 
    xm_p = sqrt(rho)*xm_p_unscaled;
    np = length(xm_p);
    
    disp(['Count number: ' num2str(numm_cnt)]);
    M = N*dmax./tp +(np*N) +1;
    for i = 1:nbrOfRealizations
        h_true = g_list(:,i);
        d_true = d_list(i);
        % Generate y
        Y= zeros(M,L);
%         [X_real,e_real] = gen_X_e(d_true,ts,tp,xm_p,dmax);
        [X_real,e_real] = gen_X_e_data(d_true,ts,tp,xm_p,dmax,rho);
        for l = 1:L
            % Generate X and e
            yl_noiseless = (h_true(l)*X_real)*e_real;
            noise = sqrt(N0/ts)*randcn(M,1);
            yl = yl_noiseless + noise;
            Y(:,l) = yl;
        end
        
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
            
            w1 = 0;
            w2 = 0;
            w3 = 0;
            
            for l = 1:L
                yl = Y();
                A = (X'*yl)*(yl'*X);
                a = A(1);
                b = A(2);
                c = A(3);
                d = A(4);
                
                w1 = w1 + a./(ts^2) - (b+c)./(ts^2) + d/(ts^2);
                w2 = w2 + -2*a./ts + (b+c)./ts;
                w3 = w3 + a;
            end
            
            
            
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
            z3 = x3*w2-w3*x2;
            z2 = 2*x3*w1-2*w3*x1;
            z1 = x2*w1-x1*w2;
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
            obj_num = 0;
            for l = 1:L
                yl = Y(:,l);
                obj_num = obj_num - (abs(mu'*yl).^2);
            end
            obj(k) = obj_num/(norm(mu)^2);%-(abs(mu'*y).^2)/(norm(mu)^2);
        end
        
        [~,ind] = min(obj);
        d_hat = d_set(ind);
        
        % Channel estimation
        
        [X_hat,e_hat] = gen_X_e(d_hat,ts,tp,xm_p,dmax);
        mu_hat = X_hat*e_hat;
        h_hat = zeros(L,1);
        for l = 1:L
            yl = Y(:,l);
            h_hat(l) = (mu_hat'*yl)./(mu_hat'*mu_hat);
        end
        
        
        mu_PerfSync = X_real*e_real;
        
        for l = 1:L
            yl = Y(:,l);
            h_PerfSync_hat = (mu_PerfSync'*yl)./(mu_PerfSync'*mu_PerfSync);
            h_err = abs(h_true(l)-h_hat(l)).^2;
            
            %h_PerfSync_err = abs(h_true(l)-h_PerfSync_hat).^2;
            %h_PerfSync_err_list(numm_cnt,i,l) = h_PerfSync_err ;
            
            G_hat_list(numm_cnt,i,l) = h_hat(l);
            G_PerfSync_hat_list(numm_cnt,i,l) = h_PerfSync_hat ;
            
            h_err_list(numm_cnt,i,l) = h_err ;
            
            
        end
        
        %         h_PerfSync_hat = (mu_PerfSync'*y)./(mu_PerfSync'*mu_PerfSync);
        %
%         h_err = abs(h_true-h_hat).^2;
        d_err = abs(d_true-d_hat).^2./tp;
        D_hat_list(numm_cnt,i) = d_hat;
        d_err_list(numm_cnt,i) = d_err ;
        
%         h_PerfSync_err = abs(h_true-h_PerfSync_hat).^2;
% 
%         h_err_list(numm_cnt,i) = h_err ;
%         h_PerfSync_err_list(numm_cnt,i) = h_PerfSync_err ;
%         
%         
%         G_hat_list(numm_cnt,i) = h_hat;
%         G_PerfSync_hat_list(numm_cnt,i) = h_PerfSync_hat ; 
        
   end    
end

avg_error_h = mean(mean(h_err_list,3),2);
avg_error_d = mean(d_err_list,2);


end

