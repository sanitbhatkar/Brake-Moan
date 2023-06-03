clc
clearvars



%% User inputs

D = 4;     % Number of variables
Np = 10;    % Population size
T = 5;     % Max iterations
omegam = 0.51691;  % Radial Velocity
f_measured = 7.5691; % Measured friction torque

%Bounds
l_torque = 0;     % Lower bound
u_torque = 10;     % Uppper bound

l_vs = 0;
u_vs = 0.1;

%% Initial matrix

lb = NaN(1,D);
ub = NaN(1,D);

lb(1:D-1) = l_torque*ones(1,D-1);      % Lower bound matrix
ub(1:D-1) = u_torque*ones(1,D-1);       % Upper bound matrix

lb(1,D) = l_vs;
ub(1,D) = u_vs;

% Solution vector
f = NaN(Np,1);

% Population vector
P = NaN(Np,D);
P(:,1:D) = repmat(lb,Np,1)+repmat(ub-lb,Np,1).*rand(Np,D);  % Formula guarantees that the population remains within bounds

fc = P(:,1);
fs = P(:,2);
sigma2 = P(:,3);
vs = P(:,4);

for i = 1:Np
    f(i) = abs(f_measured - lugre_static_fn(fc(i),fs(i),sigma2(i),vs(i),omegam));
end

f;
fbest_itr = NaN(1,T+1);
fbest_itr(1) = min(f);

%% Iteration loop

for j = 1:T
 
    for i = 1:Np
    
        %% Teacher Phase

        Xmean = mean(P);    % Mean population vector

        [~,ind]= min(f);    % Index of best solution
        Xbest = P(ind,:);   % Best solution vector X

        TF = randi([1,2]);  % Randomly selected teaching parameter
        Xnew = P(i,:)+rand(1,D).*(Xbest-Xmean*TF);

        % Bounding
        Xnew = max(lb,Xnew);    % Bounding lower bound
        Xnew = min(ub,Xnew);   % Bounding upper bound

        % New solution calculation
##        fnew = fun(Xnew)
        fnew = abs(f_measured - lugre_static_fn(Xnew(1),Xnew(2),Xnew(3),Xnew(4),omegam));

        % Greedy selection
        if (fnew < f(i))
            P(i,:) = Xnew;
            f(i) = fnew;
        end

        %% Learner Phase

        % Random partner selection
        ind = randi([1,Np]);

        % To ensure that unique partner is selected
        while ind == i
            ind = randi([1,Np]);
        end

        if f(i) < f(ind)
            Xnew = P(i,:)+rand(1,D).*(P(i,:)-P(ind,:));
        else
            Xnew = P(i,:)-rand(1,D).*(P(i,:)-P(ind,:));
        end
        
        % Bounding
        Xnew = max(lb,Xnew);    % Bounding lower bound
        Xnew = min(ub,Xnew);    % Bounding upper bound

        % New solution calculation
##        fnew = fun(Xnew)
        fnew = abs(f_measured - lugre_static_fn(Xnew(1),Xnew(2),Xnew(3),Xnew(4),omegam));

        % Greedy selection
        if (fnew < f(i))
            P(i,:) = Xnew;
            f(i) = fnew;
        end
    

    end
    
    % Best value after every iteration
    fbest_itr(j+1) = min(f);
    
end    

[fbest,ind] = min(f);
Pbest = P(ind,:)

figure (1)
semilogy(0:T,fbest_itr)
grid on






