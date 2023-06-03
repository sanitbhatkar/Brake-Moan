clc
clearvars

tic

N_repeat = 20;

for i = 1:N_repeat

##--------------------------- User inputs START --------------------------------

% Heuristic algorithm inputs

D = 4;                    % Number of variables
Np = 50;                   % Population size
T = 200;                   % Max iterations
Pc = 0.8;                 % Cross-over probability            
F = 0.85;                 % Mutation probability


%Bounds

l_torque = 0;           % Lower bound
u_torque = 10;          % Uppper bound

l_sigma2 = 0;           % Lugre Viscous Damping Lower Bound
u_sigma2 = 10;           % Lugre Viscous Damping Upper Bound

l_vs = 0;               % Sliding Velocity Lower Bound
u_vs = 0.1;             % Sliding Velocity Upper Bound

##--------------------------- User inputs END ----------------------------------




##------------------------ Excel Raw Data Import START -------------------------

pkg load io
Raw = xlsread ('Raw.xlsx');
% Raw = xlsread ('Raw.xls', 'Sheet1', 'A1:B73');

[sz,~] = size(Raw);

itr = 0;

for i = 1:sz
  
  if Raw(i,1) < 0
    
    itr = itr + 1;
    
  endif
  
 end 
 
 % Initialization
 
 Raw_neg = zeros(itr-1,2);
 Raw_pos = zeros(sz-itr,2);
 
 Raw_neg(:,1) = Raw(1:itr-1,1);
 Raw_neg(:,2) = Raw(1:itr-1,2);
 
 Raw_pos(:,1) = Raw(itr+1:sz,1);
 Raw_pos(:,2) = Raw(itr+1:sz,2);
 
## % Plot
## 
## figure(1)
## 
## plot(Raw_pos(:,1),Raw_pos(:,2),'*')
## grid on
## set(gca, "fontsize", 18)
## xlabel('Velocity (r/s)')
## ylabel('Torque (N-m)')
## title('Stribeck Curve')
 
## figure(2)
## 
## plot(Raw_neg(:,1),Raw_neg(:,2),'*')
## grid on
## set(gca, "fontsize", 18)
## xlabel('Velocity (r/s)')
## ylabel('Torque (N-m)')
## title('Stribeck Curve')
 
## figure(3)
##  plot(Raw(:,1),Raw(:,2),'*')
## grid on
## set(gca, "fontsize", 18)
## xlabel('Velocity (r/s)')
## ylabel('Torque (N-m)')
## title('Stribeck Curve')


##------------------------ Excel Raw Data Import END ---------------------------




##------------------------ Matrix Initialization START -------------------------
%% Initial matrix

lb = NaN(1,D);
ub = NaN(1,D);

lb(1:D-2) = l_torque*ones(1,D-2);             % Lower bound matrix
ub(1:D-2) = u_torque*ones(1,D-2);             % Upper bound matrix

lb(1,D) = l_vs;                               % Sliding Velocity Lower Bound
ub(1,D) = u_vs;                               % Sliding Velocity Upper Bound

lb(1,D-1) = l_sigma2;                                % Lugre Viscous Damping Lower Bound
ub(1,D-1) = u_sigma2;                                % Lugre Viscous Damping Upper Bound



% Solution vector
f = NaN(Np,1);
fu = NaN(Np,1);
U = NaN(Np,D);

% Population vector
P = NaN(Np,D);
% Formula guarantees that the population remains within bounds
P(:,1:D) = repmat(lb,Np,1)+repmat(ub-lb,Np,1).*rand(Np,D);  


##fc = P(:,1);
##fs = P(:,2);
##sigma2 = P(:,3);
##vs = P(:,4);

for i = 1:Np
    f(i) = lugre_static_fn1(P(i,:),Raw_pos);
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
        fu(k) = lugre_static_fn1(U(k,:),Raw_pos); 
        
        if fu(k) < f(k)
            
            f(k) = fu(k);
            P(k,:) = U(k,:);
        end
        
    end
             
    % Best value after every iteration
    fbest_itr(t+1) = min(f);
    
end 


[fbest,ind] = min(f);
Pbest = P(ind,:);

##----------------------- Optimization : DE END --------------------------------




##---------------------- Estimated Value Calculation START ---------------------

fc = Pbest(1,1);
fs = Pbest(1,2);
sigma2 = Pbest(1,3);
vs = Pbest(1,4);

omegam = Raw_pos(:,1);
f_measured = Raw_pos(:,2);

[sm,~] = size(omegam);

for i = 1:sm

 A = omegam(i)/vs;
 B = fc + (fs-fc)*exp(-A^2);
 C = sigma2*omegam(i);
 f_estimated(i) = B+C;
 
end


## --------------------------------- Results -----------------------------------
 
## figure(1)
## 
## plot(Raw_pos(:,1),Raw_pos(:,2),'*')
## hold on
## plot(omegam,f_estimated,'-r')
## hold off
## grid on
## set(gca, "fontsize", 18)
## xlabel('Velocity (r/s)')
## ylabel('Friction Torque (N-m)')
## title('Stribeck Curve')
 

 figure (2)
 semilogy(0:T,fbest_itr*1e-3)
 hold on
 grid on
 set(gca, "fontsize", 18)
 xlabel('Number of Iterations')
 ylabel('|T_{Measured} - T_{Estimated}|')
 title('Estimation Error')

end

hold off

 toc

