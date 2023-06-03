clc
clearvars


##--------------------------- User inputs START --------------------------------

% Heuristic algorithm inputs

D = 2;                        % Number of variables
Np = 5;                       % Population size
T = 10;                        % Max iterations
Pc = 0.8;                     % Cross-over probability            
F = 0.85;                     % Mutation probability


%Bounds

l_sigma_0 = 0;                        % Lower bound
u_sigma_0 = 5000;                     % Uppper bound

l_sigma_1 = 0;                        % Lugre Viscous Damping Lower Bound
u_sigma_1 = 100;                      % Lugre Viscous Damping Upper Bound

##--------------------------- User inputs END ----------------------------------




%% ----------------- Input of static parameters START --------------------------

Tc = 4.352;                             % Sliding Torque
Ts = 8.802;                             % Static Torque
omega_s = 0.0613;                       % Sliding speed
sigma_2 = 6.416;                        % LuGre parameters
X = [10,10];                            % Dynamic parameter matrix

%% ------------------- Input of static parameters END --------------------------




##------------------------ Excel Raw Data Import START -------------------------

pkg load io
Raw_hyst = xlsread ('Raw_hyst.xlsx');
% Raw = xlsread ('Raw.xls', 'Sheet1', 'A1:B73');

[sz,~] = size(Raw_hyst);


 % Initialization

time = Raw_hyst(:,1);                      % Time 
omega_m = Raw_hyst(:,2);                   % Motor velocity
f_measured = Raw_hyst(:,3);                % Friction Torque

dt = Raw_hyst(2,1) - Raw_hyst(1,1);        % Time step
 
## % Plot
## 
## figure(1)
## 
## plot(omega_m,Torque,'*')
## grid on
## set(gca, "fontsize", 18)
## xlabel('Velocity (r/s)')
## ylabel('Torque (N-m)')
## title('Hysterisis Curve')
 
##------------------------ Excel Raw Data Import END ---------------------------




##------------------------ Matrix Initialization START -------------------------
%% Initial matrix

lb = NaN(1,D);
ub = NaN(1,D);

lb(1:D-1) = l_sigma_0*ones(1,D-1);             % Lower bound matrix sigma_0
ub(1:D-1) = u_sigma_0*ones(1,D-1);             % Upper bound matrix sigma_0

lb(1,D) = l_sigma_1;                           % Lower bound matrix sigma_1
ub(1,D) = u_sigma_1;                           % Upeer bound matrix sigma_1

% Solution vector
f = NaN(Np,1);
fu = NaN(Np,1);
U = NaN(Np,D);

% Population vector
P = NaN(Np,D);
% Formula guarantees that the population remains within bounds
P(:,1:D) = repmat(lb,Np,1)+repmat(ub-lb,Np,1).*rand(Np,D);  


for i = 1:Np
    f(i) = lugre_dynamic_fn(Tc,Ts,omega_s,sigma_2,time,omega_m,f_measured,dt,P(i,:));
end

##f;
fbest_itr = NaN(1,T+1);
fbest_itr(1) = min(f);

##------------------------ Matrix Initialization END ---------------------------




##----------------------- Optimization : DE START ------------------------------

%% Iteration loop

for t = 1:T
 
    for i = 1:Np
    
        %% Mutation
        
        cand = [1:i-1 i+1:Np];
        idx = cand(randperm(Np-1,3));
        
        X1 = P(idx(1),:);
        X2 = P(idx(2),:);
        X3 = P(idx(3),:);
        
        V = X1 + F*(X2-X3);
        
        %% Crossover
        
        % Uniform crossover
        
        delta = randi(D,1);
        
        for k = 1:D
            
            if (rand <= Pc || k == delta)
                U(i,k)= V(k);
            else
                U(i,k)=P(i,k);
            end
        end     

    end
    
    %% Bounding and greedy selection
    
    for k = 1:Np
        
        U(k,:)= max(U(k,:),lb) ;        % Lower bound
        U(k,:)= min(U(k,:),ub) ;        % Upper bound
        
##        fu(k) = fun(U(k,:));
        fu(k) = lugre_dynamic_fn(Tc,Ts,omega_s,sigma_2,time,omega_m,f_measured,dt,U(k,:));
        
        if fu(k) < f(k)
            
            f(k) = fu(k);
            P(k,:) = U(k,:);
        end
        
    end
             
    % Best value after every iteration
    fbest_itr(t+1) = min(f);
    
end 


[fbest,ind] = min(f);
Pbest = P(ind,:)

##----------------------- Optimization : DE END --------------------------------




##---------------------- Estimated Value Calculation START ---------------------

g = NaN(sz,1);
sigma_0 = Pbest(1,1);
sigma_1 = Pbest(1,2);

##--------------------------------- Equation 1 ---------------------------------

for i=1:sz

g(i) = Tc + (Ts - Tc)*exp(-(omega_m(i)/omega_s)^2);
g(i) = g(i)/sigma_0;

end

% Initial Defletion at t = 0
x0 = [0];
z = zeros(1,sz);

##--------------------------------- Equation 2 ---------------------------------

% Bristle deflection calculation

for i = 2:sz
  
[~,dum] = ode45(@bristle_deflect,[0:dt:time(i)],x0,omega_m(i),g(i));
z(i) = dum(end,1);
zdot(i) = omega_m(i) - abs(omega_m(i))*(z(i)/g(i));

end

##for i = 2:sz
##  
##[~,dum] = ode45(@bristle_deflect,[time(i-1):dt/5:time(i)],z(i-1),omega_m(i),g(i));
##z(i) = dum(end,1);
##zdot(i) = omega_m(i) - abs(omega_m(i))*(z(i)/g(i));
##
##end

##--------------------------------- Equation 3 ---------------------------------
err_fn = 0;

for i =  1:sz

f_estimated(i) = sigma_0*z(i) + sigma_1*zdot(i) + sigma_2*omega_m(i);

end

##---------------------- Estimated Value Calculation END -----------------------




## --------------------------------- Results -----------------------------------
 
 figure(3)
 
 plot(omega_m,f_measured,'*')
 hold on
 plot(omega_m,f_estimated,'-r')
 hold off
 grid on
 set(gca, "fontsize", 18)
 xlabel('Velocity (r/s)')
 ylabel('Friction Torque (N-m)')
 title('Hysterisis Curve')
 legend('Measured Torque','Estimated Torque')
 

 figure (4)
 semilogy(0:T,fbest_itr*1e-3)
 grid on
 set(gca, "fontsize", 18)
 xlabel('Number of Iterations')
 ylabel('|T_{Measured} - T_{Estimated}|')
 title('Estimation Error')