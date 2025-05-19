% Jonathan Kocur
% Computing for Physicists
% HW 5 Problem 2A

% Shooting Method

% Here we are defining the empty arrays for each of the y values
x = [];
y1 = [];
y2 = [];
y3 = [];
y4 = [];

% Here we are defining the initial conditions for both shots of the
% shooting method and the x values
y1(1) = 0;
y2(1) = 0;
y3(1) = 0;
y4(1) = 0;

x(1) = 0;
dx = 0.001;
N = 1/dx;

% This for loop calculates the shooting method where we use runge-kutta
% twice with different initial conditions, then a final function is found
% with the upper boudary by finding lamda
for n = 2:N
    x(n) = x(n-1) + dx;
    xx = x(n-1);

    K1 = dx * y2(n-1);
    L1 = dx * (-4 * (y1(n-1) - xx));
    K2 = dx * (y2(n-1) + L1);
    L2 = dx * (-4 * ((y1(n-1) + K1) - (xx + dx)));
    
    y1(n) = y1(n-1) + 0.5 * (K1 + K2);
    y2(n) = y2(n-1) + 0.5 * (L1 + L2);
    y3(n) = y3(n-1) + 0.5 * (K1 + K2);
    y4(n) = y4(n-1) + 0.5 * (L1 + L2);
end

y1(1) = 0;
y2(1) = 1;
x(1) = 0;

% This for loop updates the values of y1 and y2 for the second shot of the
% shooting method
for n = 2:N
    x(n) = x(n-1) + dx;
    xx = x(n-1);

    K1 = dx * y2(n-1);
    L1 = dx * (-4 * (y1(n-1) - xx));
    K2 = dx * (y2(n-1) + L1);
    L2 = dx * (-4 * ((y1(n-1) + K1) - (xx + dx)));
    
    y1(n) = y1(n-1) + 0.5 * (K1 + K2);
    y2(n) = y2(n-1) + 0.5 * (L1 + L2);
end

% Here lamda is calculated and the analytical solution of the function is
% caluculated for comparison with the shooting method
lamda = ((2 - y3(N))/(y1(N) - y3(N)));
y_analytical = sin(2*x)/sin(2) + x;

% This loop calculates the final function for y with the calculated value
% of lamda and y1 and y3
for n = 1:N
    y_final(n) = lamda * y1(n) + (1 - lamda) * y3(n);
end

figure(1)
plot(x,y1)
hold on
plot(x,y3)
hold on
plot(x,y_final)
hold on
plot(x,y_analytical)
legend('y1','y3','y_final','y_analytical')
hold off

