function Complex_Gaussian_out = randcn(row,col,R)
if(nargin  == 2 )
   Complex_Gaussian_out = 1/sqrt(2).*(randn(row,col) + 1i.*randn(row,col)) ; 
else
   Complex_Gaussian_out = 1/sqrt(2).*(randn(col,row)*R + 1i.*randn(col,row)*R) ; 
   Complex_Gaussian_out = Complex_Gaussian_out.' ; 
end

end