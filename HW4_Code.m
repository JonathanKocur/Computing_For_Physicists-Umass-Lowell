% Jonathan Kocur
% Computing For Physacists
% HW 4 Code

% Here the 5 coding problems will be evaluated. Below are each of the
% function calls to solve each of the problems.
ProbA13PT = GQuadrature1(-1,1,3);
ProbA15PT = GQuadrature1(-1,1,5);
ProbA2 = GQuadrature1(0,10,3);
ProbA3 = GQuadrature2(0,1,3);
ProbB1 = GQuadrature3(3,6,2);
ProbB23PT = GQuadrature4(0.5,1,3);
ProbB25PT = GQuadrature4(0.5,1,5);
ProbB3 = GQuadrature5(-1,1,3);
ProbB4 = GQuadrature6(-1,1,3);
ProbB5 = GQuadrature7(-1,1,-1,1,3);

% These functions represent the calculation of the gaussian quadrature,
% where an if statement is used to check what the amount of points used in
% the sum, then the integral is evaluated and summed

function ret = GQuadrature1(a,b,n) % Problems A1 - A2
    ret = 0;
    x = zeros(1,10);
    A = zeros(1,10);

    if n == 3
        A(1) = 5/9;
        A(2) = 8/9;
        A(3) = 5/9;
        x(1) = sqrt(3/5);
        x(2) = 0;
        x(3) = -sqrt(3/5);
    end

    if n == 5
        A(1) = 128/225;
        A(2) = (322+13*sqrt(70))/900;
        A(3) = (322+13*sqrt(70))/900;
        A(4) = (322-13*sqrt(70))/900;
        A(5) = (322-13*sqrt(70))/900;

        x(1) = 0;
        x(2) = 1/3*sqrt(5-2*sqrt(10/7));
        x(3) = -1/3*sqrt(5-2*sqrt(10/7));
        x(4) = 1/3*sqrt(5+2*sqrt(10/7));
        x(5) = -1/3*sqrt(5+2*sqrt(10/7));
    end

    % This for loop is a sum of the A values and the function values for
    % each polynomial degree
    for i = 1:n
        ret = ret + A(i) * A1Function((b-a)/2*x(i) + (b+a)/2);
    end
    ret = ret * (b-a)/2;
end

function ret = GQuadrature2(a,b,n) % Problems A3 - A4
    ret = 0;
    x = zeros(1,10);
    A = zeros(1,10);
    
    if n == 3
        A(1) = 5/9;
        A(2) = 8/9;
        A(3) = 5/9;
        x(1) = sqrt(3/5);
        x(2) = 0;
        x(3) = -sqrt(3/5);
    end

    % This for loop is a sum of the A values and the function values for
    % each polynomial degree
    for i = 1:n
        ret = ret + A(i) * A3Function((b-a)/2*x(i) + (b+a)/2);
    end
    ret = ret * (b-a)/2;
end

function ret = GQuadrature3(a,b,n) % Problem B1
    ret = 0;
    x = zeros(1,10);
    A = zeros(1,10);
    
    % This if statement declares the values of x1 and x2 based on the
    % polynomial degree, then the function can be properly evaluated
    if n == 2
        A(1) = 1;
        A(2) = 1;
        x(1) = sqrt(1/3);
        x(2) = -sqrt(1/3);
    end

    % This for loop is a sum of the A values and the function values for
    % each polynomial degree
    for i = 1:n
        ret = ret + A(i) * B1Function((b-a)/2*x(i) + (b+a)/2);
    end
    ret = ret * (b-a)/2;
end

function ret = GQuadrature4(a,b,n) % Problem B2
    ret = 0;
    x = zeros(1,10);
    A = zeros(1,10);

    if n == 3
        A(1) = 5/9;
        A(2) = 8/9;
        A(3) = 5/9;
        x(1) = sqrt(3/5);
        x(2) = 0;
        x(3) = -sqrt(3/5);
    end

    if n == 5
        A(1) = 128/225;
        A(2) = (322+13*sqrt(70))/900;
        A(3) = (322+13*sqrt(70))/900;
        A(4) = (322-13*sqrt(70))/900;
        A(5) = (322-13*sqrt(70))/900;

        x(1) = 0;
        x(2) = 1/3*sqrt(5-2*sqrt(10/7));
        x(3) = -1/3*sqrt(5-2*sqrt(10/7));
        x(4) = 1/3*sqrt(5+2*sqrt(10/7));
        x(5) = -1/3*sqrt(5+2*sqrt(10/7));
    end

    % This for loop is a sum of the A values and the function values for
    % each polynomial degree
    for i = 1:n
        ret = ret + A(i) * B2Function((b-a)/2*x(i) + (b+a)/2);
    end
    ret = ret * (b-a)/2;
end

function ret = GQuadrature5(a,b,n) % Problem B3
    ret = 0;
    x = zeros(1,10);
    A = zeros(1,10);

    if n == 3
        A(1) = 0.804914;
        A(2) = 0.813128;
        A(3) = 0.813128;
        x(1) = 0.524648;
        x(2) = 1.650680;
        x(3) = -1.650680;
    end

    % This for loop is a sum of the A values and the function values for
    % each polynomial degree
    for i = 1:n
        ret = ret + A(i) * B3Function((b-a)/2*x(i) + (b+a)/2);
    end
    ret = ret * (b-a)/2;
end

function ret = GQuadrature6(a,b,n) % Problem B4
    ret = 0;
    x = zeros(1,10);
    A = zeros(1,10);

    if n == 3
        A(1) = 5/9;
        A(2) = 8/9;
        A(3) = 5/9;
        x(1) = sqrt(3/5);
        x(2) = 0;
        x(3) = -sqrt(3/5);
    end

    % This for loop is a sum of the A values and the function values for
    % each polynomial degree
    for i = 1:n
        ret = ret + A(i) * B4Function((b-a)/2*x(i) + (b+a)/2);
    end
    ret = ret * (b-a)/2;
end

function ret = GQuadrature7(a,b,c,d,n) % Problem B5
    ret = 0;
    x = zeros(1,10);
    y = zeros(1,10);
    A = zeros(1,10);
    gx = zeros(1,10);

    if n == 3
        A(1) = 5/9;
        A(2) = 8/9;
        A(3) = 5/9;
        x(1) = sqrt(3/5);
        x(2) = 0;
        x(3) = -sqrt(3/5);
        y(1) = sqrt(3/5);
        y(2) = 0;
        y(3) = -sqrt(3/5);
    end

    % This for loop is a sum of the A values and the function values for
    % each polynomial degree
    for j = 1:n
        for i = 1:n
            gx(i+1) = (gx(i) + A(j) * B5Function(x(i),y(j)))*((d-c)/2);
        end
    end

    for j = 1:n
        for i = 1:n
            ret = ((b-a)/2)*(ret + A(i)*gx(i));
        end
    end
end


% These functions represent the users function to be tested above using the
% gaussian quadrature
function out = A1Function(x)
    out = exp(-x^2);
end
function out = A3Function(x)
    out = exp(-x^2 + x - 2);
end
function out = B1Function(x)
    out = x^3 + 2*x^2 + 1;
end
function out = B2Function(x)
    out = sin(sqrt(x))*sqrt(x);
end
function out = B3Function(x)
    out = x^3 + x^2 + x + 1;
end
function out = B4Function(x)
    out = (1 - x^2)^2;
end
function out = B5Function(x,y)
    out = x^2*y^2;
end
