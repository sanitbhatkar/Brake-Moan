function op = RK_fun(y0,t_start,t_end,omega_max,t0,t1,t3,tmax,Xs,X)


y0;                       % Initial Condition
t_start;                  % Interval Start
t_end;                    % Interval End
% omega_m;                  % Value of omega at t = t_start
% g;                        % LuGre parameter


h = (t_end-t_start)/25;            % Time step
t = t_start:h:t_end;               % t goes from t_start to t_end seconds.

ystar = zeros(size(t));         % Preallocate array (good coding practice)

ystar(1) = y0;                  % Initial condition gives solution at t=t_start.


for i=1:(length(t)-1)
  
  dum = omg_fn(t(i),omega_max,t0,t1,t3,tmax,Xs,X);
  omega_m = dum(1);
  g = dum(2);
  
  k1 = omega_m - abs(omega_m)*(ystar(i)/g);     % Approx for y gives approx for deriv
  y1 = ystar(i)+k1*h/2;                         % Intermediate value (using k1)
  
  k2 = -2*y1;                                    % Approx deriv at intermediate value.
  y2 = ystar(i)+k2*h/2;                         % Intermediate value (using k2)
  
  k3 = -2*y2;                                    % Another approx deriv at intermediate value.
  y3 = ystar(i)+k3*h;                           % Endpoint value (using k3)
  
  k4 = -2*y3 ;                                   % Approx deriv at endpoint value.
  
  ystar(i+1) = ystar(i) + (k1+2*k2+2*k3+k4)*h/6; % Approx soln
  
end


% plot(t,yexact,t,ystar,'*');
% legend('Exact','Approximate');

op = ystar(end);