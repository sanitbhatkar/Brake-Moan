
function op = lugre_dynamic_fn_rk(time,f_measured,omega_max,t0,t1,t3,tmax,Xs,X)

%% ----------------- Input of static parameters START --------------------------

Tc = Xs(1,1);                               % Sliding Torque
Ts = Xs(1,2);                               % Static Torque
omega_s = Xs(1,3);                          % Sliding speed
sigma_2 = Xs(1,4);                          % LuGre parameters
% X = [10,10];                              % Dynamic parameter matrix

%% ------------------- Input of static parameters END --------------------------




%% -------------------- Dynamic prameter calculation START ----------------------

% g = NaN(sz,1);
sigma_0 = X(1,1);
sigma_1 = X(1,2);

%% --------------------------------- Equation 1 ---------------------------------

% for i=1:sz
% 
% g(i) = Tc + (Ts - Tc)*exp(-(omega_m(i)/omega_s)^2);
% g(i) = g(i)/sigma_0;
% 
% end

% Initial Defletion at t = 0

[sz,~] = size(time);
z = zeros(1,sz);
zdot = zeros(1,sz);
f_estimated = zeros(1,sz);
err_itr = zeros(1,sz);

%% --------------------------------- Equation 2 ---------------------------------

% Bristle deflection calculation

% for i = 2:sz
%   
% [~,dum] = ode45(@bristle_deflect,[0:dt:time(i)],x0,omega_m(i),g(i));
% z(i) = dum(end,1);
% zdot(i) = omega_m(i) - abs(omega_m(i))*(z(i)/g(i));
% 
% end

% for i = 2:sz
%   
% [~,dum] = ode45(@bristle_deflect,[time(i-1):dt/5:time(i)],z(i-1),omega_m(i),g(i));
% z(i) = dum(end,1);
% zdot(i) = omega_m(i) - abs(omega_m(i))*(z(i)/g(i));
% 
% end



for i = 2:sz

z(i) = RK_fun(z(i-1),time(i-1),time(i),omega_max,t0,t1,t3,tmax,Xs,X);
dum = omg_fn(time(i),omega_max,t0,t1,t3,tmax,Xs,X);
omega_m = dum(1);
g = dum(2);
zdot(i) = omega_m - abs(omega_m)*(z(i)/g);

end


%% --------------------------------- Equation 3 ---------------------------------
err_fn = 0;

for i =  1:sz

dum = omg_fn(time(i),omega_max,t0,t1,t3,tmax,Xs,X);
omega_m = dum(1);
f_estimated(i) = (sigma_0*z(i) + sigma_1*zdot(i))*sign(omega_m) + sigma_2*omega_m;

%%  ------------------------ objective Function Definition ----------------------
 
 % Error for each population
 
 err_itr(i) = abs(f_measured(i) - f_estimated(i));
 err_fn = (err_itr(i))^2 + err_fn;
 
end

err_fn = 0.5*err_fn;

op = err_fn;

