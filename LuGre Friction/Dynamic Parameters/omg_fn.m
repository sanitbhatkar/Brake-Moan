function op = omg_fn(time,omega_max,t0,t1,t3,tmax,Xs,X)

% omega_max = 0.0265;                               % Max Rotor Velocity
% 
% t0 = 0;
% t1 = 60;
% t2 = 120;
% t3 = 180;
% tmax = 240;

omega_m_local = NaN(1,1);

%% ----------------- Input of static parameters START --------------------------

Tc = Xs(1,1);                               % Sliding Torque
Ts = Xs(1,2);                               % Static Torque
omega_s = Xs(1,3);                          % Sliding speed
sigma_2 = Xs(1,4);                          % LuGre parameters
% X = [10,10];                              % Dynamic parameter matrix

%% ------------------- Input of static parameters END --------------------------




%% ------------------------ Input Velocity Cycle START -------------------------
 
  if (time <= t1)
    
    omega_m_local = omega_max*(time/t1);
    
  elseif ( time <= t3 && time > t1)
    
    omega_m_local = omega_max - 2*omega_max*(time-t1)/(t3-t1);
    
%   elseif ( time <= t3 && time > t2)
%     
%     omega_m_local = -omega_max*(time-t2)/(t3-t2);
    
  else 
    
    omega_m_local = -omega_max + omega_max*(time-t3)/(tmax-t3);
    
  end
   
%% ------------------------- Input Velocity Cycle END --------------------------




%% --------------------- LuGre parameter calculation START ----------------------


g = NaN(1,1);
sigma_0 = X(1,1);

%% --------------------------------- Equation 1 ---------------------------------


g = Tc + (Ts - Tc)*exp(-(omega_m_local/omega_s)^2);
g = g/sigma_0;

%% ---------------------- LuGre parameter calculation END -----------------------


op = [omega_m_local,g];


