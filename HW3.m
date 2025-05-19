% Jonathan Kocur
% 2/20/23
% HW 3

clear all
close all

% Here is the defined data set that will be interpolated
x = [1,2,3,4,5];
y = [3.4,4.6,1.8,6.3,1.9];
N = length(x); 
n = N-1; 
h = (x(N)-x(1))/n;

size(x);
N = size(x,2);

% Here all of our 0 vectors and arrays are established. These will be used
% to be filled in for solving our SLE
K = zeros(N-2,N-2);
dx = zeros(1,N-1);
dy = zeros(1,N-1);
k = zeros(1,N-2);
b = zeros(1,N-2);

% Here is where we are first going through and calculating our dx and dy
% and dy values
for i = 1:N-1
    dx(i) = x(i+1) - x(i);
    dy(i) = y(i+1) - y(i);
end

% This for loop fills in the matrix for the cubic spline. First the
% diagonal elements are filled. Then, if statements are used to fill in the
% remaining spots. We use the if statements because the pattern is
% different in the first and last rows.
for i = 1:N-2    
    K(i,i) = 2*(dx(i)+dx(i+1));
    if i == 1
        K(i,i+1) = dx(i+1);
    elseif i == N-2
        K(i,i-1) = dx(i);
    else
        K(i,i+1) = dx(i+1);
        K(i,i-1) = dx(i);
    end    
    b(i) = 6*((dy(i+1) / dx(i+1)) - (dy(i) / dx(i)));
end

% Next we will use gaussian elimination to for the k values which will then
% be used to find the piecewise equations to fit the data
lambda = zeros(3);
Knew = K;
bnew = b;
for i=1:1:N-2
    for j=i+1:1:3
        lambda(j,i)=Knew(j,i)/Knew(i,i);
        for k=i:1:N-2
            Knew(j,k)=Knew(j,k)-lambda(j,i)*Knew(i,k);
        end
        bnew(j) = bnew(j) - lambda(j,i)*bnew(i);
    end
end

% Here is where our k values are found and our k vector is defined. This
% will be used to find our C and D coefficients and create the piecewise
% functions that will fit our data set
k3 = bnew(3) / Knew(3,3);
k2 = (bnew(2) - (Knew(2,3) * k3)) / Knew(3,3);
k1 = (bnew(1) - (Knew(1,2) * k2)) / Knew(1,1);
k = [k1,k2,k3]';
knew = [0; k; 0];


% This loop is used to calculate the coefficients that will be used in the 
% interpolated function 
for i = 1:n
    D(i) = y(i);
    B(i) = knew(i)/2;
    A(i) = (knew(i+1)-knew(i))/(6*h);
    C(i) = (y(i+1)-y(i))/h-h/6*(2*knew(i)+knew(i+1));
end

% Here is where the x plotting information is defined. This is what will
% determine how many x points will be used to run through the piecewise
% functions to them be plotted with the spline
r = 100;
hh = h/r;
xx = x(1):hh:x(N);

% This nested for loop is used to calculate the function values for the
% cubic spline using the calculated coefficients
for i = 1:n
    for j = r*(i-1) + 1:r*i
        f(j) = A(i)*(xx(j)-x(i))^3 + B(i)*(xx(j)-x(i))^2 + C(i)*(xx(j)-x(i)) + D(i);
    end
end

f(r*n+1) = y(N);

% Here is the plot of the data set and the interpolated line from the
% cubic spline
plot(x, y,'bo')
hold on
plot(xx,f)
hold off
xlabel('x')
ylabel('y')
title('Interpolated Data Set')
legend('Data Points','Interpolated Line')