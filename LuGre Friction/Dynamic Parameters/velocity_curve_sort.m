clc
clearvars

%% -------------------------- Excel Inut START ---------------------------------
% pkg load io
filename = 'Raw_hyst.xlsx';
Raw_hyst = xlsread(filename);

[sz,~] = size(Raw_hyst);

[omega_max,idx1] = max(Raw_hyst(:,2));
[omega_min,idx2]  = min(Raw_hyst(:,2));

t0 = Raw_hyst(1,1)
t1 = Raw_hyst(idx1,1)
t3 = Raw_hyst(idx2,1)
tmax = Raw_hyst(end,1)

%% -------------------------- Excel Inut END -----------------------------------




%% ----------------- Input of static parameters START --------------------------

Tc = 4.352;                             % Sliding Torque
Ts = 8.802;                             % Static Torque
omega_s = 0.0613;                       % Sliding speed
sigma_2 = 6.416;                        % LuGre parameters

Xs = [Tc,Ts,omega_s,sigma_2];           % Static parameters
X = [10,10];                            % Dynamic parameter matrix

%% ------------------- Input of static parameters END --------------------------



