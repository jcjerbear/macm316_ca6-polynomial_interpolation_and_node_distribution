function err = polyinterperr(f,x)

% Inputs:
% x - the nodes for interpolation
% f - the function to interpolate

% Output:
% err - the maximum error max_x | f(x) - P(x) | over [-1,1]
% Note: err is computed on a grid of 10,000 equispaced points

% Convert x to a column vector
[r,s] = size(x);
if s>r
    x = x';
end

% construct fine grid for evaluating error on
m =10000;
grid = (linspace(-1,1,m))';

y = f(x); % sample f at the nodes x
w=baryweights(x); % compute the barycentric weights w
u=baryinterp(x,w,y,grid); % compute the interpolating polynomial P(x)

err = max(abs(f(grid) - u));
end

% function to compute barycentric weights w
function w = baryweights(x)

n = length(x);
w = zeros(n,1);
for i = 1:n
    X = x-x(i)*ones(n,1);
    X = X([1:i-1 i+1:n],1);
    w(i) = 1/prod(X);
end

end

% function to evaluate P(x) on a grid using the barycentric formula
function u = baryinterp(x,w,y,grid)

n = length(x);
m = length(grid);
u = zeros(m,1);

for i = 1:m
    
    diff = grid(i).*ones(n,1)-x;
    l = sum(diff==0);
    
    if l==0
        z = w./diff;
        u(i) = (z'*y)/sum(z);
    else
        u(i) = y(find(diff==0));
    end
end

end


