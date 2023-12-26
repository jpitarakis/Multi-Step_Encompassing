function lrv = andrews_lrv(e)
   [t,~] = size(e);

   a  =  e(2:t)\e(1:t-1); 
   k = 0.8; 
   c_left = 1.1447*(((4*(a^2)*t)/(((1+a)^2)*((1-a)^2)))^(1/3));
   c_right = 1.1447*(((4*(k^2)*t)/(((1+k)^2)*((1-k)^2)))^(1/3));
   veco = [c_left;c_right];
   l  =  min(veco);
   l = fix(l);
   lrv = e'*e/t;
   
   for i=1:l
       w = (1-(i/(1+l)));
       lrv = lrv + 2*e(1:t-i)'*e(1+i:t)*w/t;
   end
   
   %i = 1;
   %do until i>l; 
   %    w = (1-i/(l+1)); 
   %     lrv = lrv+2*e(1:t-i)'*e(1+i:t)*w/t;
   %     i = i+1;
   % endo;
   
   

