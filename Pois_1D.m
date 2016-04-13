function [ V ] = Pois_1D( x, rho, BCb, BCt)
%x   - mesh to analyze (must have even step size)
%rho - charge distribution
%BCb - bottom contact voltage (x=0)
%BCt - top contact voltage


global eps0  

nx = length(x);
dx = x(2)-x(1) ;        %update later for arbitrary meshes
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

%check to make sure ends of rho are set to 0 (essential for setting
%boundary conditions)
if(rho(1)|| rho(end))
    rho(1) =0;
    rho(end) = 0;
end
    
%sometimes check the condition number of G if having trouble converging
%condest(G);
V =  G\((dx^2/eps0)*rho' + B');


end

