
function op = lugre_dynamic_fn(time,f_measured,dt,X,Xs,omega_max,t0,t1,t3,tmax)



%% ----------------- Input of static parameters START --------------------------

Tc = Xs(1,1);                               % Sliding Torque
Ts = Xs(1,2);                               % Static Torque
omega_s = Xs(1,3);                          % Sliding speed
sigma_2 = Xs(1,4);                          % LuGre parameters
% X = [10,10];                              % Dynamic parameter matrix

%% ------------------- Input of static parameters END --------------------------




%% -------------------- Dynamic prameter calculation START ----------------------

[sz,~] = size(time);
##g = NaN(sz,1);
sigma_0 = X(1,1);
sigma_1 = X(1,2);

%% --------------------------------- Equation 1 ---------------------------------

##for i=1:sz
##
##g(i) = Tc + (Ts - Tc)*exp(-(omega_m(i)/omega_s)^2);
##g(i) = g(i)/sigma_0;
##
##end
##
##% Initial Defletion at t = 0
##x0 = [0];
z = zeros(1,sz);

%% --------------------------------- Equation 2 ---------------------------------

% Bristle deflection calculation

for i = 2:sz
  
dum = omg_fn(time(i),omega_max,t0,t1,t3,tmax,Xs,X);
omega_m = dum(1);
g = dum(2);
  
##[~,dum] = ode45(@bristle_deflect,[0:dt:time(i)],x0,omega_m(i),g(i));
##z(i) = dum(end,1);
##zdot(i) = omega_m(i) - abs(omega_m(i))*(z(i)/g(i));

[~,dum1] = ode45(@bristle_deflect,[time(i-1):dt/15:time(i)],z(i-1),omega_m,g);
z(i) = dum1(end,1);
zdot(i) = omega_m - abs(omega_m)*(z(i)/g);

end



%% --------------------------------- Equation 3 ---------------------------------
err_fn = 0;

for i =  1:sz
  
dum = omg_fn(time(i),omega_max,t0,t1,t3,tmax,Xs,X);
omega_m = dum(1);

f_estimated(i) = sigma_0*z(i) + sigma_1*zdot(i) + sigma_2*omega_m;

%%  ------------------------ objective Function Definition ----------------------
 
 % Error for each population
 
 err_itr(i) = abs(f_measured(i) - f_estimated(i));
 err_fn = (err_itr(i))^2 + err_fn;
 
end

err_fn = 0.5*err_fn;

op = err_fn;

