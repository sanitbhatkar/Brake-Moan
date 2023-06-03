clc
clearvars

% Solve y'(t)=-2y(t) with y0=3, 4th order Runge Kutta
y0 = 3;                  % Initial Condition
h=0.2;                   % Time step
t = 0:h:2;               % t goes from 0 to 2 seconds.
yexact = 3*exp(-2*t)     % Exact solution (in general we won't know this
ystar = zeros(size(t));  % Preallocate array (good coding practice)

ystar(1) = y0;           % Initial condition gives solution at t=0.
for i=1:(length(t)-1)
  
  k1 = -2*ystar(i)  % Approx for y gives approx for deriv
  y1 = ystar(i)+k1*h/2;      % Intermediate value (using k1)
  
  k2 = -2*y1        % Approx deriv at intermediate value.
  y2 = ystar(i)+k2*h/2;      % Intermediate value (using k2)
  
  k3 = -2*y2        % Another approx deriv at intermediate value.
  y3 = ystar(i)+k3*h;        % Endpoint value (using k3)
  
  k4 = -2*y3        % Approx deriv at endpoint value.
  
  ystar(i+1) = ystar(i) + (k1+2*k2+2*k3+k4)*h/6; % Approx soln
end
plot(t,yexact,t,ystar,'*');
legend('Exact','Approximate');