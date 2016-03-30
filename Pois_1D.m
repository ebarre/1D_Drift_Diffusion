function [ V ] = Pois_1D(rho, BCb, BCt, x)
%rho - charge distribution
%BCb - bottom contact voltage (x=0)
%BCt - top contact voltage
%x   - mesh to analyze (must have even step size)

global eps0  

nx = length(x);
dx = x(2)-x(1);         %update later for arbitrary meshes
G = sparse(nx,nx);      % matrix for calculating V distribution
B = zeros(1,nx);        % boundary conditions

%BCb = 0;                %bottom contact voltage
%BCt = 1;                %top contact voltage

for n=1:nx
    if n==1
        G(n,n) =1;
        B(n) = BCb;
    elseif n==nx 
        G(n,n)  = 1;
        B(n) = BCt;
    else
        G(n,n-1) =1;
        G(n,n)  = -2;
        G(n,n+1) =1;
    end
    
end

V =  G\((dx^2/eps0)*rho' + B');


end

