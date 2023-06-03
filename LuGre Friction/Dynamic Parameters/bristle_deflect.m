function op = bristle_deflect(t,x0,omega_m,g)
  
##  x0
##  omega_m
##  g
  op = omega_m - abs(omega_m)*(x0(1)/g);
  
  
  