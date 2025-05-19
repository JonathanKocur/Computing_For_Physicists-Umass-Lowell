% Jonathan Kocur
% 2/5/23
% Computing For Physacists HW1

% 1)

% This function manually multiplies matrices by creating two sample
% matrices and creating a for loop that goes through each matrix element
% and multiplies and sums them.

A = [6,-2,1;-2,7,2;1,2,-5];
B = [1,1,1;1,1,1;1,1,1];

if size(A) ~= size(B)
    disp('Matrices are not compatible')
end

N = 3;
product = zeros(N,N);

for i=1:N
    for j=1:N
        product(i,j) = 0;
        for k=1:N
            product(i,j) = product(i,j) + A(i,k) * B(k,j);
        end
    end
end

MatrixMulTest = A * B;
MatrixMul=product

% 2)

% This function uses the gauss-seidel algorithm to solve a linear system of
% equations in matrix form. It starts with an initial guess for the x value
% and updates the given equations then retests until the x is found within
% a given amount of error.

x = zeros(3,1);
b = [11,5,-1];
n = size(x,1);
nmax = 8;
iter = 0;

while iter < nmax
    x_old = x;    
    for i = 1:n        
        guess = 0;              
        for j = 1:i-1
            guess = guess + A(i,j) * x(j);
        end        
        for j=i+1:n
            guess = guess + A(i,j) * x_old(j);
        end               
        x(i) = (1/A(i,i)) * (b(i) - guess);       
    end
    iter = iter + 1;   
end

GaussSeidelOutput = x

% 3)

% This function uses for loops that scan through any given matrix, it
% starts with the first value of the maxtrix and replaces it if the next
% value is larger.

maxval = A(1,1);

for i = 1:length(A)
    for j = 1:length(A)
        if A(i,j) > maxval
            maxval = A(i,j);
            indices = [i,j];
        end
    end
end

MaxOutput = [maxval,indices]
