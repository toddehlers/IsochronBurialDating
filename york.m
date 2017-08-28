
function [b,sigmab,a,sigmaa,diag] = york(data);

% york.m
%
% Slope of straight-line fit according to York(1966) using
% inverse-square-error weighting.  
%
% Syntax: [b,sigmab,a,sigmaa,diag] = york(data);
%
% Input argument data must have fields data.x, data.y, data,dx, data.dy
% All must be equal size vectors
%
% Output arguments are 
% b, slope; sigmab, error in slope
% a, intercept; sigmaa, error in intercept
% diag, diagnostics structure
%   has fields: diag.numits.
%
% york.m from Greg Balco 2007



n = length(data.x);

if n < 2;
    error('york.m -- less than 2 data -- stopping');
end;

if n < 3;
    disp('york.m -- only 2 data -- line is fully determined');
    % If only two data, slope and intercept are determined
    % So fitting scheme is not necessary...short-circuit calculation
    % Apply formula
    b = (data.y(2)-data.y(1))./(data.x(2)-data.x(1));
    a = data.y(1) - b.*data.x(1);
    % propagate errors
    dbdx = (data.y(2)-data.y(1))./((data.x(2)-data.x(1)).^2);
    dbdy = 1/(data.x(2)-data.x(1));
    sigmab = sqrt( (dbdx.*data.dx(1)).^2 + (dbdx.*data.dx(2)).^2 + ...
        (dbdy.*data.dy(1)).^2 + (dbdy.*data.dy(2)).^2 );
    dadx1 = -b + data.x(1).*dbdx;
    dadx2 = data.x(1).*dbdx;
    dady1 = 1 + data.x(1).*dbdy;
    dady2 = -data.x(1).*dbdy;
    sigmaa = sqrt( (dadx1.*data.dx(1)).^2 + (dadx2.*data.dx(2)).^2 + ...
        (dady1.*data.dy(1)).^2 + (dady2.*data.dy(2)).^2 );
    diag.numits = 0;
    return;
end;

% Precalculate some terms
xbar = mean(data.x);
ybar = mean(data.y);
u = data.x - xbar;
v = data.y - ybar;
 
sumv2  = sum(v.^2);
sumu2 = sum(u.^2);
sumuv = sum(u.*v);

% Initial guess for b estimated from 'major axis' fit
% York (1966) Equation 1

b = ( sumv2 - sumu2 + ((sumv2-sumu2)^2 + 4.*(sum(u.*v).^2)).^(0.5) )/(2 * sumuv);

% Weights for the measurements:
% Inverse squared error weighting.

wx = 1./data.dx.^2;
wy = 1./data.dy.^2;

% Obtain the solution of the cubic equation slope by iterating.
% b calculated above serves as the initial guess.
% Calculate the important root as a function of b; replace with new
% value of b; repeat. 

tol = b.*1e-5; % Tolerance 1e-5 of slope
its = 0; % Iteration counter

while 1;
    % the term w is a function of the slope b; calculate it now
    % York p. 1081, Eqn. 14.5
    w = (wx.*wy)./((b.^2).*wy + wx);

    % Now calculate the new slope
    % York p. 1084, halfway down
    % Various parameters
    
    alpha = (2.*sum( (w.^2).*u.*v./wx) )./(3.*sum( (w.^2).*(u.^2)./wx ));
    beta = (sum( (w.^2).*(v.^2)./wx ) - sum(w.*(u.^2)))./(3.*sum( (w.^2).*(u.^2)./wx ));
    gamma = -sum(w.*u.*v)./sum( (w.^2) .* (u.^2) ./wx );
    cosphi = (alpha.^3 - (3.*alpha.*beta./2) + (gamma./2))./((alpha.^2 - beta).^(3/2));
    phi = acos(cosphi);
    
    % Now calculate the new solution
    b3 = alpha + 2.*sqrt(alpha.^2 - beta).*(cos((phi + 2.*2.*pi)./3));
    % Check for convergence
    if abs(b3-b) < abs(tol);
        b = b3; % Use the current value
        break; % bail
    end;
    b = b3; % Replace with new slope
    its = its + 1; % Increment iteration counter
end;


% Uncertainty in b. York p. 1084 bottom. 
sigma2b = (1/(n-2)).*(sum(w.*((b3.*u - v).^2)))./(sum(w.*(u.^2)));
sigmab = sqrt(sigma2b);

% Intercept a. York Equation 18. 

a = ybar - b.*xbar;

% Uncertainty in a. York p. 1085 top. 

sigma2a = sigma2b.*( sum(w.*(data.x.^2)) )./sum(w);
sigmaa = sqrt(sigma2a);

% out assignments

diag.numits = its;