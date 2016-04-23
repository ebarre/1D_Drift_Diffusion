function [ J, p ] = current1D_noRecomb(x, V,  pl, pr )
%This gives the current density and charge distribution for a 1D voltage
%distribuion
%assuming boundary condition for leftmost side (i.e. set left to 0 V)
%Inputs
% V: array of voltages
% pl: boundary condition for hole concentration at the left contact
% pr: boundary condition for hole concentration at the right contact

global Dp mup q0

%step size
nx = length(x);
dx = x(2) -x(1);
E(1:nx-1) = -(V(2:nx) - V(1:nx-1))/dx; %E(n) is actually E(n+1/2)
%figure(3)
%plot(E)
J = zeros(1,nx);

A = mup.*E/2 - Dp/dx ;
B = mup.*E/2 + Dp/dx ;
G = sparse(nx);
for n = 1:nx
    if n == 1
        G(n,n) = 1; %Boundary condition for leftmost segment
    elseif n == nx
        G(n,n) = 1;
    else
        G(n,n+1) = A(n);
        G(n,n) = B(n) - A(n-1);
        G(n,n-1) = -B(n-1);
    end
end
condest(G)
BC = zeros(1,nx);
%Ohmic Boundary conditions
BC(1) = pl;      % p0 at left end
BC(end) = pr;    % p0 at right end
%sometimes check the condition number of G if having trouble converging
%condest(G)
p = G\BC';
gradp = (p(2:nx)-p(1:nx-1))/dx;
avgp = (p(2:nx)+p(1:nx-1))/2;
J = q0.*(mup.*avgp'.*E-Dp.*gradp');
end

