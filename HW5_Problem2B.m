% Jonathan Kocur
% Computing for Physicists
% HW 5 Problem 2B

% Here we are initializing the amount of points and preallocating the
% arrays of x and y
m = 5;
x = zeros(1,m+1);
dx = 1/m;
y = zeros(1,m+1);

% Here we are creating the A matrix and b vector and filling in the initial
% conditions
b = zeros(1,m+1);
b(1) = 0;
b(m+1) = 2;
A = zeros(m+1,m+1);
A(1,1) = 1;
A(m+1,m+1) = 1;


% Here we are creating the xvector based on the amount of points we want
% with the given interval from the initial condition, then filling the A
% matrix and b vector with the values from the system of equations
for n = 1:m
    x(n+1) = x(n) + dx;
end
for n = 2:m
    A(n,n-1) = 1;
    A(n,n) = 4*dx^2 - 2;
    A(n,n+1) = 1;
    b(n) = 4*x(n)*dx^2;
end
b = b';


% Here we are performing gaussian elimination to solve the matrix A with
% vector b which will give us the resulting y points
N = m+1;
Anew = A;
Bnew = b;
lambda = zeros(N);
for i=1:1:N
    for j=i+1:1:N
        lambda(j,i)=Anew(j,i)/Anew(i,i);
        for k=i:1:N
            Anew(j,k)=Anew(j,k)-lambda(j,i)*Anew(i,k);
        end
        Bnew(j) = Bnew(j) - lambda(j,i)*Bnew(i);
    end
end
U = Anew;
Bnew2 = Bnew;
for i = m+1:-1:1
    for k = m+1:-1:i+1
        Bnew2(i) = Bnew2(i) - y(k)*U(i,k);
    end
    y(i) = Bnew2(i)/U(i,i);
end

y_analytical = sin(2*x)/sin(2) + x;
plot(x,y)
hold on
plot(x,y_analytical)
legend('y','y_analytical')
hold off