function Mat_out = T2Mat(Mat_in)

[size1, size2, size3 ]= size(Mat_in) ; 

if(size3 == 1)
   Mat_out = Mat_in ; 
else
   if(size2 == 1)
      Mat_out = reshape(Mat_in,size1,size3) ; 
   else
      Mat_out = reshape(Mat_in,size2,size3) ; 
   end
end

end