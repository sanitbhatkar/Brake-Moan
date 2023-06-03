
function op = lugre_static_fn(fc,fs,sigma2,vs,omegam)
  
 A = omegam/vs;
 B = fc + (fs-fc)*exp(-A^2);
 C = sigma2*omegam;
 
 op = B+C;