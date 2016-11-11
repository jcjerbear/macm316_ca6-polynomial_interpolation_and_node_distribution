% Assignment: MACM 316 Computing Assignment 6
% Title: Polynomial interpolation and node distribution
% Description: Examine how the locations of the nodes x0,x1,...,Xn affect
%              the accuracy and robustness of polynomial interpolation
% Author: Jerry Chen
% File name: ca06.m


% Part1: Polynomial interpolation with equally spaced nodes
% Define f(x)
f = @(x) exp(sin(5*x+0.5));

% Parameters
N = 100;
errLinear = zeros(0,N);
n = 1:N;

for i=1:N
    i;
    x = linspace(-1,1,i);
    errLinear(i) = polyinterperr(f,x);
end

figure(1)
semilogy(n,errLinear)
title('Error of Polynomial Interpolation At Equidistant Nodes','fontsize',12)
xlabel('Node','fontsize',10)
ylabel('log_{10}Error','fontsize',10)
grid on


% Part2: Compute nodes using greedy algorithm
% Parameters
M = 10000;
xTilda = linspace(-1,1,M);
m = 1:N+1;
result = ones(1,N);
for i=1:N
    result(1)=-1;
    for j=2:i+1
        Vx = ones(1,M);
        for k=1:M
            for l=1:j-1
                Vx(k) = Vx(k)*abs(xTilda(k)-result(l));
            end
        end
        [argvalue,argmax] = max(Vx);
        result(j) = xTilda(argmax);
        Xsort = sort(result);
    end
end
%{
result = ones(1,N+1);
result(1) = -1;
v = ones(1,M);
m = 1:M;

for j=1:N+1
    j;
    for k=1:M
        k;
        v(k) = v(k)*abs(xTilda(k)-result(j));
        %v = v.*value;
    end
    [argvalue,argmax] = max(v);
    result(j+1) = xTilda(argmax);
end
resultSorted = sort(result);
%}
figure(2)
plot(n,x,m,Xsort)
axis([0 100 -1 1])
title('Linear and Nonlinear Nodes','fontsize',14)
xlabel('Node','fontsize',12)
ylabel('Value of Node','fontsize',12)
grid on


% Part3: Polynomial interpolation with greedy algorithm nodes
errGreedy = ones(1,N);
for i=1:N
    i;
    x = result(1:i);
    errGreedy(i) = polyinterperr(f,x);
end

figure(3)
semilogy(n,errGreedy)
title('Error of Polynomial Interpolation At Greedy Nodes','fontsize',12)
xlabel('Iteration','fontsize',10)
ylabel('log_{10}Error','fontsize',10)
grid on


% Comparison of both errors
figure(4)
semilogy(n,errLinear,n,errGreedy)
title('Error of Polynomial Interpolations','fontsize',12)
xlabel('Node','fontsize',10)
ylabel('log_{10}Error','fontsize',10)
grid on


    