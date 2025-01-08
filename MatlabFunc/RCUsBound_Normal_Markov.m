function [tail_approx] = RCUsBound_Normal_Markov(constel,rho,nd,nc,L,s,alpha,h,h_est,R)
constel_labels = unique(nchoosek([constel,constel],2), 'rows') ;
t = length(constel_labels) ; 
constel_len = length(constel) ;  

a = 1-alpha; 
b = alpha; 

R_nats = R.* log(2) ; 
M = 2.^((nc.*L).*R); 

% Generate P(s) (discrete transitition probabilities)
initProb = ones(1,t) .* 1/t ;
P = ones(t);  
for ii  = 1 : t 
   xi_vec = constel_labels(ii,:) ; 
   for jj = 1 : t
      xj_vec = constel_labels(jj,:) ; 
      if(xi_vec(2) == xj_vec(1))
         P(ii,jj) = 1/constel_len ; 
      else
         P(ii,jj) = 0 ; 
      end
   end
end

zeta_sel = 0 ; 
% if(zeta_sel <= zeta_down+1e-3)
%    disp('WARNING') ; 
%    disp(zeta_sel)
% end
for L_cnt = 1:L
   h_in = h(L_cnt) ; 
   h_est_in = h_est(L_cnt) ; 
   a_in = a(L_cnt);  
   b_in = b(L_cnt) ;
   for ii = 1:2
       xi_vec = constel_labels(ii,:) ; 
       [M_0(1,ii), M_deriv(1,ii),M_deriv2(1,ii) ] = MGF_NumericIntegral_infoDens(h_in,h_est_in,xi_vec(1),xi_vec(2),constel,rho,s,a_in,b_in,zeta_sel) ;
   end
   M_0(1,3:4) = fliplr(M_0(1,1:2)) ;
   M_deriv(1,3:4) = fliplr(M_deriv(1,1:2)) ;
   M_deriv2(1,3:4) = fliplr(M_deriv2(1,1:2)) ;

   for ii = 1 : t
       xi_vec = constel_labels(ii,:) ; 
       for jj = 1:t
          xj_vec = constel_labels(jj,:) ;      
          if(xj_vec(1) == xi_vec(2))
             varPhi(ii,jj) = M_0(1,jj)  ;
          else
             varPhi(ii,jj) = 0  ;
          end
       end
   end

   initProb_zeta = initProb .* M_0 ; 
   a_theta = varPhi(1,1,:) ; b_theta = varPhi(1,2,:) ; 
   kappa_n_zeta_vec(L_cnt) = T2Vect(log(2.^(-1-nd).*(a_theta+b_theta).^(1+nd))); 
   a_theta_prime = M_deriv(1,1,:) ; b_theta_prime = M_deriv(1,2,:); 
   a_theta_prime2 = M_deriv2(1,1,:) ; b_theta_prime2 = M_deriv2(1,2,:); 
   kappa_n_zeta_div1_vec(L_cnt) = T2Vect((1+nd) .* (a_theta_prime + b_theta_prime)./(a_theta + b_theta)) ;
   kappa_n_zeta_div2_vec(L_cnt) = T2Vect(-1.*(1+nd) .* (a_theta_prime + b_theta_prime).^2./(a_theta + b_theta).^2 + (1+nd) .* (a_theta_prime2 + b_theta_prime2)./(a_theta + b_theta)) ;  
end

kappa_n_zeta = sum(kappa_n_zeta_vec) ; 
kappa_n_zeta_div1 = sum(kappa_n_zeta_div1_vec) ; 
kappa_n_zeta_div2 = sum(kappa_n_zeta_div2_vec) ; 

mu_n_zeta = kappa_n_zeta_div1/nd ; 
sigma_n_zeta = kappa_n_zeta_div2 ./nd ;

Is = -mu_n_zeta ; 
Vs = sigma_n_zeta; 

if(M == Inf)
   tail_approx = qfunc((nd.* Is - n.*R.*log(2))./ (sqrt(nd.*Vs))) ;
else
   tail_approx = qfunc((nd.* Is - log(M-1))./ (sqrt(nd.*Vs))) ;
end




if(isnan(tail_approx))
%    disp('NAN_WARNING')
   tail_approx = exp(kappa_n_zeta + nd.*R_nats.*L) ; 
end
if(isinf(tail_approx))
   tail_approx = 1 ;
%    disp('INF_WARNING')
end
end







