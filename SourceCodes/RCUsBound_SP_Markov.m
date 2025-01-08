function [tail_approx] = RCUsBound_SP_Markov(constel,rho,nd,nc,L,s,alpha,h,h_est,R)
constel_labels = unique(nchoosek([constel,constel],2), 'rows') ;
t = length(constel_labels) ; 
constel_len = length(constel) ;  

a = 1-alpha; 
b = alpha; 

R_nats = R.* log(2) ; 

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

zeta_up = 3 ; 
zeta_down = -1.5; 
zeta_sel = 1 ; 

for ii  = 1 : 20
   rate_sel = Transform_Zeta2Rate(constel,constel_labels,rho,nd,nc,L,s,a,b,h,h_est,zeta_sel) ; 
   rateDiffdebug = abs(rate_sel - R_nats) ; 
   if(rateDiffdebug <= 1e-3)
      break ; 
   elseif(rate_sel <= R_nats)
      zeta_up = zeta_sel ; 
      zeta_sel = (zeta_down + zeta_sel) / 2 ; 
   else
      if(zeta_sel >= 1) % If zeta>=1 and we are gonna increase even more
         zeta_sel = 1; 
         break ;
      end
      zeta_down = zeta_sel ; 
      zeta_sel = (zeta_up + zeta_sel) / 2 ; 
   end
end
zeta_sel(zeta_sel>=1) = 1 ;
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

temp_n  = zeta_sel.*sqrt(nd.*sigma_n_zeta) ; 
temp_n2 = (1-zeta_sel).*sqrt(nd.*sigma_n_zeta) ; 
temp_n3  = -zeta_sel.*sqrt(nd.*sigma_n_zeta) ;
B_knot = @(u) exp((u.^2)./2).*(qfunc(u));
Psi_tilde = @(a1,a2,b_in) exp((a1).*(-kappa_n_zeta_div1 - b_in.*R_nats + kappa_n_zeta_div2./2) .* qfunc(a1.*sqrt(kappa_n_zeta_div2) - a2 .* (kappa_n_zeta_div1 + b_in*R_nats)./(kappa_n_zeta_div2))) ; 
if(zeta_sel < 0)
   tail_approx = 1 - exp(kappa_n_zeta - nd.*zeta_sel.*mu_n_zeta).* (B_knot(temp_n3) - B_knot(temp_n2)) ; 
   tail_approx(isnan(tail_approx)) = 1 ; % To address rare numerical inconsistencies.
   tail_approx(isinf(tail_approx)) = 1 ; 
elseif(zeta_sel >= 1)
   tail_approx = exp(kappa_n_zeta + nd.*R_nats.*L).* ( Psi_tilde(1,1,nd) + Psi_tilde(0,-1,nd) ) ; 
else
   tail_approx = exp(kappa_n_zeta - nd.*zeta_sel.*mu_n_zeta)  .* (B_knot(temp_n)  + B_knot(temp_n2)  )  ; 
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

function [rate] = Transform_Zeta2Rate(constel,constel_labels,rho,nd,nc,L,s,a_in,b_in,h_in,h_est_in,zeta)
% To do extend it for QPSK. 
kappa_n_zeta_div1 = zeros(1,L) ;
for L_looper = 1 : L 
   h = h_in(L_looper);  
   h_est = h_est_in(L_looper) ;
   a = a_in(L_looper) ; 
   b = b_in(L_looper) ; 
   for ii = 1:2
      [x_vec] = constel_labels(ii,:) ; 
      xi = x_vec(1) ; xj = x_vec(2) ; 
      xk_i = sqrt(rho).*constel ;
      xk_j = sqrt(rho).*constel ;
      xi = sqrt(rho).*xi ; 
      xj = sqrt(rho).*xj ; 
      mu_ij = h.*a.*xj + h.*b.*xi ; 
      mu_ijj = mu_ij - h_est.*xj ; 
      mu_ijk = mu_ij - h_est.*xk_i; 
      mu_ijk = mu_ijk(:); 
      funcExpectDeriv = @(zr,zi) 0.5.^(zeta) .* exp(-1.*zi.^2-1.*zr.^2 + s.*zeta.*abs(mu_ijj+i.*zi+zr).^2) .*... 
                      (exp(-s.*abs(mu_ijk(1) + i.*zi +zr).^2) + exp(-s.*abs(mu_ijk(2) + i.*zi +zr).^2) ).^(zeta) ... 
                      .* (-0.220636 + 0.31831.*s.*abs(mu_ijj + i.*zi + zr ).^2 + 0.31831.* log(exp(-s.*abs(mu_ijk(1) + i.*zi +zr ).^2) + exp(-s.*abs(mu_ijk(2) + i.*zi +zr ).^2) ) ) ;
      funcExpect = @(zr,zi) exp(zeta.*s.*abs(mu_ijj+zr+ 1i.*zi).^2) .* (0.5 ...
           .* ( exp(-s.*abs(mu_ijk(1)+zr +1i.*zi).^2) + exp(-s.*abs(mu_ijk(2)+zr +1i.*zi).^2)) ).^(zeta) ... 
            .* sqrt(2)/sqrt(2.*pi) .* exp(-0.5 .* (zr.*sqrt(2)).^2 ) .* sqrt(2)/sqrt(2.*pi) .* exp(-0.5 .* (zi.*sqrt(2)).^2 ) ;
      M_0(1,ii) = simp2D(funcExpect,-5,5,-5,5,30,30) ; 
      M_deriv(1,ii) = simp2D(funcExpectDeriv,-5,5,-5,5,30,30) ;  
   %    M_0(1,ii) = quad2d(funcExpect,-5,5,-5,5,'AbsTol',0,'RelTol',1e-5,'MaxFunEvals', 100) ; 
   %    M_deriv(1,ii) = quad2d(funcExpectDeriv,-5,5,-5,5,'AbsTol',0,'RelTol',1e-5,'MaxFunEvals', 100) ;  
   end
   % M_deriv(1,3:4,zeta_cnt) = fliplr(M_deriv(1,1:2,zeta_cnt)) ;
   M_0(1,3:4) = fliplr(M_0(1,1:2)) ;

   t = length(constel_labels) ; 
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

   a_theta = varPhi(1,1) ; b_theta = varPhi(1,2) ; 
   a_theta_prime = M_deriv(1,1) ; b_theta_prime = M_deriv(1,2); 

   kappa_n_zeta_div1(L_looper) = (1+nd) .* (a_theta_prime + b_theta_prime)./(a_theta + b_theta) ;
end
mu_n_zeta = sum(kappa_n_zeta_div1)/nd ; 
rate = -mu_n_zeta.*nd./(nc.*L) ; 
end






