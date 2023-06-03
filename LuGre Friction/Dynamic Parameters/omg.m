clc
clearvars

omega_max = 0.0265;                               % Max Rotor Velocity

t0 = 0;
t1 = 60;
t2 = 120;
t3 = 180;
tmax = 240;
dt = (tmax-t0)/100;


%% ----------------- Input of static parameters START --------------------------

Tc = 4.352;                             % Sliding Torque
Ts = 8.802;                             % Static Torque
omega_s = 0.0613;                       % Sliding speed
sigma_2 = 6.416;                        % LuGre parameters
X = [10,10];                            % Dynamic parameter matrix

%% ------------------- Input of static parameters END --------------------------



time = 0:dt:tmax;
[~,sz] = size(time);
omega_m = NaN(sz,1);

%% ------------------------ Input Velocity Cycle START -------------------------


for i = 1:sz
  
  
  if (time(i) <= t1)
    
    omega_m(i) = omega_max*(time(i)/t1);
    
  elseif ( time(i) <= t3 && time(i) > t1)
    
    omega_m(i) = omega_max - 2*omega_max*(time(i)-t1)/(t3-t1);
    
%   elseif ( time(i) <= t3 && time(i) > t2)
%     
%     omega_m(i) = -omega_max*(time(i)-t2)/(t3-t2);
    
  else 
    
    omega_m(i) = -omega_max + omega_max*(time(i)-t3)/(tmax-t3);
    
   end
   
end

figure(1)
plot(time,omega_m)
grid on
% set(gca,'fontsize', 18)
xlabel('Time (s)')
ylabel('Velocity (r/s)')
title('Velocity Curve')

%%-------------------------- Input Velocity Cycle END --------------------------




%%--------------------- LuGre parameter calculation START ----------------------

[sz,~] = size(time);
g = NaN(sz,1);
sigma_0 = X(1,1);

%%--------------------------------- Equation 1 ---------------------------------

for i=1:sz

g(i) = Tc + (Ts - Tc)*exp(-(omega_m(i)/omega_s)^2);
g(i) = g(i)/sigma_0;

end    

%%---------------------- LuGre parameter calculation END -----------------------

