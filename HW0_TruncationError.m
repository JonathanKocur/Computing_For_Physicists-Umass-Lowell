% 1/18/2023
% HW 0 Truncation Error

% First we're going to clear the workspace
clear all
close all

% Now we'll define the input variables, the given xvalue and the number of
% iterations N in the loop 
x=0.69;
N=10;

% Here we will call the function, calculate the difference between the
% taylor expansion and the actual sine
a=JonMan(N,x);
diff = abs(sin(x)-a);

% Finally we will plot the difference as a function of N
Nindex = 0:N;
plot(Nindex,diff)
xlabel('N value')
ylabel('Error')
title('Taylor Expansion Error by N value')


% Now we define the function that performs the taylor expansion of sine
function returnvalue=JonMan(N,x)
    y = zeros(1,N+1);
    y(1) = 0;

    for i=0:N-1
        y(i+2) = y(i+1) + (-1).^i/factorial(2*i+1)*x.^(2*i+1);        
    end

    returnvalue = y;
end