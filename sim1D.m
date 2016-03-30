%For a simple 1D equation with no generation and Dirichlet
%Boundary conditions calculates the current/charge distribution/voltage

%takes a while to converge when using ridiculous input but eventually does
close all

global eps0 k T kT q0 tau mup m hbar dx
%constants
eps0 = 8.854e-12;           % F/m  permittivity of free space
hbar = 6.682119e-16;        % eVs
q0 =1.609e-19;              % C
k = 8.61734318e-5;          % eV/K Boltzmann constant
T = 300;                    % Temperature in Kelvin
kT = k*T;                   % eV
tau = 1e-7;                 % Estimate a hole lifetime of about 10 ns
mup = 100*(1/100)^2;        %cm^2/Vs*(1m/100cm)^2 
m0 = 5.68562975e-12;        % eV*s^2/m^2 electron mass
m = 0.5*m0;                 % Assuming an effective mass (Look for value)

%choice of variables
nx = 5000;                  % steps in x direction
dx = 1e-12;                  %spatial step is defined as 1nm
x = linspace(0, dx*nx, nx); %m


%initial guess for carrier density is linear and only considering holes
Ef = 0.5;                   % eV setting a fermi level rel to valence
Vleft = 0;
Vright = 0.2; 

%The density of holes for left and right contact based on parabolic 
%assumption (+degeneracy of K & K') + MaxwellBoltzmann Approx of
%Fermi-level
pl = 2*(m*kT/(2*pi*hbar^2))^(3/2)*exp(-(Ef)/kT); 
pr = 2*(m*kT/(2*pi*hbar^2))^(3/2)*exp(-(Ef)/kT); 
%p0 = 1e12*(1e2)^3;          %1e12 carriers /cm^3 (1e2cm/1m)^3        
p = pl*ones(1,nx);            % the free carrier density
%p = (1e8*x-50).^2*10*(1e2)^3 +p0; %weird initial guess 
%rho = 100*exp(-((x-mean(x))/(100*dx)).^2);
rho = q0*p; %need to think about this one a bit more
rho(1) = 0; %i define these so that my boundary conditions are correct
rho(nx) = 0; %i define these so that my boundary conditions are correct

maxIter =1e3;
V = Pois_1D(rho, Vleft, Vright, x);
%plotting V& p 1 of 3
%{
%figure(4)
%plot(x, p, 'Color', [0.5,0.5,0])
%ylim([pl pr])
%hold on
figure(3)
plot(x, V, 'Color', [0.5,0.5,0])
ylim([0, 0.2])
hold on
%}
for i=1:maxIter
    Vnew = Pois_1D(rho, Vleft, Vright, x);
    %rarely vary by greater than 10^14
    %all(abs(V-Vnew)<1e-14)
    %pause(1)
    V = Vnew;
    %p0 = 2.*(m*kT/(2*pi*hbar^2))^(3/2).*exp(-(Ef-V)./kT);
    [J, pnew] = current1D(V, pl, pr);
    
    %plotting V& p 2 of 3
    %{
    figure(4)
    color = (1/1000)*floor(1000*i/maxIter);
    plot(x,pnew-p', 'Color', [0.5,0.5,color])
    figure(3)
    plot(x, Vnew, 'Color', [0.5,0.5,color])
    pause(1) 
    %}
    if(all(abs(pnew-p')< pl/1e10))
        i
        break
    end
    p=pnew';
    rho = q0*p; 
    rho(1) = 0;
    rho(nx) = 0;
    
    if(i==maxIter)
        fprintf('Does not converge. \n')
    end
end

%plotting V& p 3 of 3
figure(3)
hold off
figure(4)
hold off
%}

figure(1)
[hAx,hLine1,hLine2] = plotyy(x*1e6, rho, x*1e6, V, 'semilogy', 'plot');
xlabel('x-axis (\mum)')
ylabel(hAx(1),'Charge density (C/m^3)') % left y-axis
ylabel(hAx(2),'Potential (V)') % right y-axis
%ylim(hAx(1), [0 1e-10])

figure(2) 
semilogy(x*1e6, p, 'r', 'linewidth', 2);
%ylim([pr pl])
xlabel('x-axis (\mum)')
ylabel('Hole concentration /m^3')
%ylim([(4/5)*p0 1.2*p0])