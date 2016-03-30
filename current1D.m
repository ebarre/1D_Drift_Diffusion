function [ J, p ] = current1D( V, pb, pt )
%This gives the current density and charge distribution for a 1D voltage
%distribuion assuming boundary condition for leftmost side (i.e. set left 
%to 0 V) and right (set righthandside to some voltage Vapp)

%Inputs
% V: array of voltages
% pb: boundary condition for concentration of holes in the bottom contact
% pt: boundary condition for concentration of holes in the top contact

global kT dx tau q0 mup

%mobility of holes
Dp = kT*mup;                %m^2/s^2

%step size
nx = length(V);
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
        G(n,n) = B(n) - A(n-1) - dx/tau;
        G(n,n-1) = -B(n-1);
    end
end

BC = zeros(nx,1);
%Ohmic Boundary conditions
BC(1) = pb;      % p0 at bottom (0 V)
BC(end) = pt;    % p0 at top (0.2 V)

RC = ones(nx,1)*(-dx*pb/tau);
RC(1) = 0;
RC(end) = 0;

p = G\(BC+RC);

gradp = (p(2:nx)-p(1:nx-1))/dx;
avgp = (p(2:nx)+p(1:nx-1))/2;
J = q0.*(mup.*avgp'.*E-Dp.*gradp');
 
end

