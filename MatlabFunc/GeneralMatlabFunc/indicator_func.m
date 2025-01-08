function [x,posInv] = indicator_func(x,IndicUpLim)
   posInv = find(x > IndicUpLim) ; 
   x(posInv) = []; 
   
end
