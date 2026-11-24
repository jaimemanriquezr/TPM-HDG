function [p, dp, ddp] = jacobi_eval(n,a,b,x)
%-------------------------------------------------------------------------
%Description: Evaluates Jacobi polynomial of degree n and parameters (a,b)
%           and its first and second derivatives at x.
%Input:     - n: degree of the polynomial
%           - a: parameter (alpha)
%           - b: parameter (beta)
%           - x: point (can be vector) for evaluation
%Output:    - p: evaluation of the Jacobi polynomial
%           - dp: evaluation of the first derivative of the Jacobi
%           polynomial
%           - ddp: evaluation of the second derivative of the Jacobi
%           polynomial
%-------------------------------------------------------------------------
p = jacobi_pol(n, a, b, x);
switch n
    case 0
        dp = zeros(size(x));
        ddp = zeros(size(x));
    case 1
        dp = .5*(n+a+b+1)*jacobi_pol(n-1,a+1,b+1,x);
        ddp = zeros(size(x));
    otherwise
        dp = .5*(n+a+b+1)*jacobi_pol(n-1,a+1,b+1,x);
        ddp = .25*(n+a+b+1)*(n+a+b+2)*jacobi_pol(n-2,a+2,b+2,x);
end
end

function P = jacobi_pol(n, alpha, beta, t)
P0 = ones(size(t));
P1 = 0.5*(2*(alpha + 1) + (alpha + beta + 2)*(t - 1));

switch n
    case 0
        P = P0;
    case 1
        P = P1;
    otherwise
        for k = 1:n-1
            ck = 2*k + alpha + beta;

            A10 = (ck + 1) * (alpha^2 - beta^2);
            A11 = ck * (ck + 1) * (ck + 2);
            A0 = 2 * (k + alpha) * (k + beta) * (ck + 2);
            L = (2*(k + 1)*((ck + 1) - k) * ck);
            
            P = ((A10 + t*A11).*P1 - A0*P0) / L;
            P0 = P1;
            P1 = P;
        end
end
end
