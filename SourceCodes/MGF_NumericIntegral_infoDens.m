function [MGF_numericApprox,MGF_numericApprox_deriv1,MGF_numericApprox_deriv2] = MGF_NumericIntegral_infoDens(h,h_est,xi,xj,constel,rho,s,alpha_not,alpha,zeta)

xk_i = sqrt(rho).*constel ;
xi = sqrt(rho).*xi ; 
xj = sqrt(rho).*xj ; 
mu_ij = h.*alpha_not.*xj + h.*alpha.*xi ; 
mu_ijj = mu_ij - h_est.*xj ; 
mu_ijk = mu_ij - h_est.*xk_i; 
mu_ijk = mu_ijk(:); 

% if(xi == xj)
%    a = h.*sqrt(rho) - h_est.*sqrt(rho) ; 
%    a_tilde = h.*sqrt(rho) + h_est.*sqrt(rho) ; 
% else
%    a = h.*sqrt(rho).*(alpha_not-alpha) - h_est.*sqrt(rho); 
%    a_tilde =  h.*sqrt(rho).*(alpha_not-alpha) + h_est.*sqrt(rho); 
% end
% sigma2 = 8.*rho.*abs(h_est).^2 ; 
if(length(constel) == 2)
%    funcExpect_temp = @(z) 1/sqrt(2.*pi.*sigma2) .*exp(-(-log(exp(-(z))-1)./s - (abs(a_tilde).^2 - abs(a).^2) ).^2 ./ (2.*sigma2)).* ...
%        1./(s.*(1-exp((z)))).* exp(-zeta.*z) ; 
   funcExpect = @(zr,zi) exp(zeta.*s.*abs(mu_ijj+zr+ 1i.*zi).^2) .* (0.5 ...
        .* ( exp(-s.*abs(mu_ijk(1)+zr +1i.*zi).^2) + exp(-s.*abs(mu_ijk(2)+zr +1i.*zi).^2)) ).^(zeta) ... 
         .* sqrt(2)/sqrt(2.*pi) .* exp(-0.5 .* (zr.*sqrt(2)).^2 ) .* sqrt(2)/sqrt(2.*pi) .* exp(-0.5 .* (zi.*sqrt(2)).^2 ) ;
   funcExpectDeriv = @(zr,zi) 0.5.^(zeta) .* exp(-1.*zi.^2-1.*zr.^2 + s.*zeta.*abs(mu_ijj+i.*zi+zr).^2) .*... 
                (exp(-s.*abs(mu_ijk(1) + i.*zi +zr).^2) + exp(-s.*abs(mu_ijk(2) + i.*zi +zr).^2) ).^(zeta) ... 
                .* (-0.220636 + 0.31831.*s.*abs(mu_ijj + i.*zi + zr ).^2 + 0.31831.* log(exp(-s.*abs(mu_ijk(1) + i.*zi +zr ).^2) + exp(-s.*abs(mu_ijk(2) + i.*zi +zr ).^2) ) ) ; 
   funcExpectDeriv2 = @(zr,zi) 1./pi .* 0.5.^(zeta) .* exp(-1.*zi.^2-1.*zr.^2+s.*zeta.*abs(mu_ijj + i.*zi +zr).^2) .* ... 
                     (exp(-s.*abs(mu_ijk(1) + i.*zi +zr).^2) + exp(-s.*abs(mu_ijk(2) + i.*zi +zr).^2)).^zeta .* (0.480453 +s.^2.*abs(mu_ijj + i.*zi +zr).^4- ...
                     1.38629 .*  log(exp(-s .* abs(mu_ijk(1) + i.*zi +zr).^2) + exp(-s.*abs(mu_ijk(2) + i.*zi +zr).^2) ) + ...
                     log(exp(-s .* abs(mu_ijk(1) + i.*zi +zr).^2) + exp(-s.*abs(mu_ijk(2) + i.*zi +zr).^2) ).^2 + s.*abs(mu_ijj + i.*zi +zr).^2 .* ...
                     (-1.3862 + 2.*log(exp(-s .* abs(mu_ijk(1) + i.*zi +zr).^2) + exp(-s.*abs(mu_ijk(2) + i.*zi +zr).^2) ) ) ) ; 
%    MGF_numericApprox_deriv1 = quad2d(funcExpectDeriv,-5,5,-5,5,'AbsTol',0,'RelTol',1e-5) ; 
%    MGF_numericApprox_deriv2 = quad2d(funcExpectDeriv2,-5,5,-5,5,'AbsTol',0,'RelTol',1e-5) ; 
   MGF_numericApprox_deriv1 = simp2D(funcExpectDeriv,-5,5,-5,5,30,30) ; 
   MGF_numericApprox_deriv2 = simp2D(funcExpectDeriv2,-5,5,-5,5,30,30) ;
else
   funcExpect = @(zr,zi) exp(zeta.*s.*abs(mu_ijj+zr+ 1i.*zi).^2) .* (0.25 ...
        .* ( (exp(-s.*abs(mu_ijk(1)+zr +1i.*zi).^2) + exp(-s.*abs(mu_ijk(2)+zr +1i.*zi).^2)) ...
         +  exp(-s.*abs(mu_ijk(3)+zr +1i.*zi).^2) +exp(-s.*abs(mu_ijk(4)+zr +1i.*zi).^2) )  ).^(zeta) ... 
         .* sqrt(2)/sqrt(2.*pi) .* exp(-0.5 .* (zr.*sqrt(2)).^2 ) .* sqrt(2)/sqrt(2.*pi) .* exp(-0.5 .* (zi.*sqrt(2)).^2 ) ;
   MGF_numericApprox_deriv1 = NaN ;
   MGF_numericApprox_deriv2 = NaN ;
end
% MGF_numericApprox = quad2d(funcExpect,-5,5,-5,5,'AbsTol',0,'RelTol',1e-5,'MaxFunEvals', 100) ; 
tic
MGF_numericApprox = simp2D(funcExpect,-5,5,-5,5,30,30) ; 


debug = 1; 

