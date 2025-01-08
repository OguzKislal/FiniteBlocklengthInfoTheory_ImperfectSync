function [err_val,failFlag] = RCUsBound_SP_Approx_Quasistatic(h_vect,h_est_vect,rho,sigmaSqr,s,n,nd,L,R,InfoDens)

err_val = 1 ;
failFlag = 0 ;
if(nargin == 9)
   InfoDens = -1 ; 
end
nc = n./L ; 

zeta_vectLen = 10 ; 

beta_A = s .* (rho .* abs(h_vect-h_est_vect).^2 + sigmaSqr) ; 
beta_B = s ./ (1+s.*rho.*abs(h_est_vect).^2) .* (rho.* abs(h_vect).^2 + sigmaSqr) ; 
v = s.^2 .* abs(  rho.*abs(h_vect).^2 + sigmaSqr - conj(h_vect) .* h_est_vect .*rho ).^2 ./ ( beta_A .* beta_B .* (1+s.*rho.*abs(h_est_vect).^2) ) ; 
zeta_term1 = sqrt((beta_B - beta_A).^2 + 4.*beta_A.*beta_B.*(1-v) ) ; 
zeta_down =  -1.*(zeta_term1 + beta_A - beta_B ) ./ (2.*beta_A.*beta_B.*(1-v)) ; 
zeta_up =  (zeta_term1 - beta_A + beta_B ) ./ (2.*beta_A.*beta_B.*(1-v)) ; 
zeta_up_smallest = min(zeta_up) ; 
zeta_down_biggest = max(zeta_down) ;

CGF = @(zeta) sum(-zeta.*log(1+s.*rho.*abs(h_est_vect).^2) - log(1+(beta_B-beta_A).*zeta-beta_A.*beta_B.*(1-v).*zeta.^2),1);
CGF_d1 = @(zeta) sum(-log(1+s.*rho.*abs(h_est_vect).^2) - ((beta_B-beta_A) - 2.*beta_A.*beta_B.*(1-v).*zeta) ./ (1+(beta_B - beta_A).*zeta - beta_A.*beta_B.*(1-v).*(zeta.^2)),1) ; 
CGF_d2 = @(zeta) sum((((beta_B-beta_A) - 2.*beta_A.*beta_B.*(1-v).*zeta) ./ (1+(beta_B - beta_A).*zeta - beta_A.*beta_B.*(1-v).*(zeta.^2))).^2 + (2.*beta_A.*beta_B.*(1-v))./(1+(beta_B - beta_A).*zeta - beta_A.*beta_B.*(1-v).*zeta.^2),1)  ;

% rateVect_zeta = (-K_dif1).*nd./(n) ; 
R_nats = R.* log(2) ; 

% zeta_stepsize = (zeta_up_smallest - zeta_down_biggest)/zeta_vectLen ; 
zeta_interval = linspace(zeta_down_biggest,zeta_up_smallest,zeta_vectLen) ; 
% zeta_interval = (zeta_down_biggest+zeta_stepsize : zeta_stepsize : zeta_up_smallest-zeta_stepsize) ;
if(zeta_down_biggest >= zeta_up_smallest)
   failFlag =1 ; 
   err_val = 1 ; 
   return ; 
end
for ii = 1 : 15
   K_dif1 = CGF_d1(zeta_interval) ; 
%    K_dif1 = CGF_d1(zeta_interval); 
   rateVect_zeta = -K_dif1.*nd/n ; 
   [rateDiff_debug, zeta_pos ]= min(abs(rateVect_zeta - R_nats)) ; 
   if(rateDiff_debug <= 1e-3)
      break ; 
   end
   zeta_pos(zeta_pos<= 1) = 2 ;
   zeta_pos(zeta_pos >= length(zeta_interval)) = length(zeta_interval)-1  ;
   zeta_down = zeta_interval(zeta_pos-1)  ;
   zeta_up = zeta_interval(zeta_pos+1) ; 
   zeta_interval = linspace(zeta_down,zeta_up,zeta_vectLen) ; 
end
zeta = zeta_interval(zeta_pos) ; 
 
if(zeta <= 1 && zeta >= 0)
   K = CGF(zeta) ;
   K_dif1 = CGF_d1(zeta) ; 
   K_dif2 = CGF_d2(zeta) ; 
   RCUs_term1exp = nd.*(K - zeta.*K_dif1) + nd.*zeta.^2/2.*K_dif2  ;
   RCUs_term1 = exp(RCUs_term1exp + log(qfunc(zeta.*sqrt(nd.*K_dif2)))) ; 
   RCUs_term2exp = nd.*(K - zeta.*K_dif1) + nd.*(1-zeta).^2/2.*K_dif2  ;
   RCUs_term2 = exp(RCUs_term2exp + log(qfunc((1-zeta).*sqrt(nd.*K_dif2)))) ; 
   err_val = RCUs_term1 + RCUs_term2 ; 
elseif(zeta > 1) % R < R_cr
   R_critical = -CGF_d1(1).*nd./n; 
   K_critical = CGF(1) ;
   K_dif1 = CGF_d1(zeta) ;
   Kd2_critical = CGF_d2(1) ; 
   RCUs_term1exp = nd.*(K_critical - K_dif1) + nd.*( R_critical - R_nats + (Kd2_critical)/2 ) ;
   RCUs_term1 = exp(RCUs_term1exp + log(qfunc(sqrt(nd.*Kd2_critical)+ nd.*(R_critical - R_nats)./(sqrt(nd.*Kd2_critical))))) ; 
   RCUs_term2exp = nd.*(K_critical - K_dif1) ;
   RCUs_term2 = exp(RCUs_term2exp + log(qfunc(-nd.*(R_critical - R_nats)./(sqrt(nd.*Kd2_critical))))) ; 
   err_val = RCUs_term1 + RCUs_term2 ;
elseif(zeta< 0)
   K = CGF(zeta) ;
   K_dif1 = CGF_d1(zeta) ; 
   K_dif2 = CGF_d2(zeta) ; 
   RCUs_term1exp = nd.*(K - zeta.*K_dif1) + nd.*(-zeta).^2/2.*K_dif2  ;
   RCUs_term1 = exp(RCUs_term1exp + log(qfunc(-zeta.*sqrt(nd.*K_dif2)))) ; 
   RCUs_term2exp = nd.*(K - zeta.*K_dif1) + nd.*(1-zeta).^2/2.*K_dif2  ;
   RCUs_term2 = exp(RCUs_term2exp + log(qfunc((1-zeta).*sqrt(nd.*K_dif2)))) ; 
   err_val = 1-(RCUs_term1 - RCUs_term2) ; 
end
if isnan(err_val)
   err_val = 1 ;
   failFlag = 1 ; 
end
err_val(err_val>1) = 1 ; 
err_val(err_val<0) = 0 ; 

%Alternative algorithm to find zeta (Although seems nice it is slower).
% CGFd1_R = @(zeta) CGF_d1(zeta).*nd./n + R_nats;
% error_zeta = 1;
% iter = 0;
% small_tmp = eps;
% while error_zeta>1e-6
%     small_tmp=small_tmp*10;
%     iter = iter+1;
%     zeta_interval = [zeta_down_biggest+small_tmp, zeta_up_smallest-small_tmp]; %[zeta_low zeta_high];
%     if isinf(CGFd1_R(zeta_interval(1))) || isinf(CGFd1_R(zeta_interval(2)))
%         continue
%     elseif sign(CGFd1_R(zeta_interval(1)))==-1 &&  sign(CGFd1_R(zeta_interval(2))) ==-1
%         continue
%     elseif sign(CGFd1_R(zeta_interval(1)))==1 &&  sign(CGFd1_R(zeta_interval(2))) ==1
%         continue
%     end
%     %options = optimset('Display','iter'); % show iterations
%     [zeta2, error_zeta] = fzero(CGFd1_R,zeta_interval);
% end
