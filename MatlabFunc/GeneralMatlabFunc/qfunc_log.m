function [logOfQfuncOut] = qfunc_log(x,t)

   if(x >= 45)
      logOfQfuncOut = -Inf ; 
      return ;
   end
   if( x<=t)
      t= x;  
   end
   phi_t = exp(-t.^2/2) ./ (sqrt(2.*pi)) ;  
   q_t = qfunc(t) ;
   temp_term = sqrt(2/pi) - (phi_t ./ q_t) ; 

   term1_log =  1./(2.*t) .* temp_term .* x.^2 ; 
   if(isnan(term1_log)) 
      term1_log = sqrt(2/pi).*1/(2.*t) - 1/2 .* x.^2 ; 
   end
   term2_log = (1./t .*log(2.*q_t) - 1/2.*temp_term) .*x  ; 
   if(isnan(term2_log))
      term2_log = (1./t .*log(2.*q_t) - 1/2.*(sqrt(2/pi) - t) )  .*x ; 
   end
   
   logOfQfuncOut = -log(2) + term1_log + term2_log ; 
   
end