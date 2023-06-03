clc
clearvars

alpha = 0:pi/50:2*pi;
theta = pi/6;

Rx = 0.03;
Ry = 0.01;

Cx = 0;
Cy = 0;

[~,sz] = size(alpha);

for i = 1:sz
  
  x(i) = Rx*cos(alpha(i))*cos(theta) - Ry*sin(alpha(i))*sin(theta) + Cx;
  y(i) = Rx*cos(alpha(i))*sin(theta) + Ry*sin(alpha(i))*cos(theta) + Cy;
  
end

figure(1)
plot(x,y*1000)
grid on
ylabel('Frtiction Torque T_{f}(t)')
xlabel('Velocity (r/s)')
set(gca, 'fontsize', 18)
title('Pre-sliding Hysterisis')
  


