
function op = lugre_static_fn1(X,Raw_pos)
  
##----------------------- LuGre Formula Parameters -----------------------------
 
fc = X(1,1);
fs = X(1,2);
sigma2 = X(1,3);
vs = X(1,4);

omegam = Raw_pos(:,1);
f_measured = Raw_pos(:,2);

[sm,~] = size(omegam);


##------------------------------ LuGre Formula ---------------------------------

err_fn = 0;

for i = 1:sm

 A = omegam(i)/vs;
 B = fc + (fs-fc)*exp(-A^2);
 C = sigma2*omegam(i);
 f_estimated = B+C;
 
 
## ------------------------ objective Function Definition ----------------------
 
 % Error for each population
 
 err_itr(i) = abs(f_measured(i) - f_estimated);
 err_fn = (err_itr(i))^2 + err_fn;
 
end

err_fn = 0.5*err_fn;
##err_fn = max(err_itr);

op = err_fn;
