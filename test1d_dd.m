%Now to compare the results if I vary the recombination rate


close all

global eps0 k T kT q0  mup m hbar Dp
%constants
%%
eps0 = 8.854e-12;           % F/m  permittivity of free space
hbar = 6.682119e-16;        % eVs
q0 =1.609e-19;              % C
k = 8.61734318e-5;          % eV/K Boltzmann constant
T = 300;                    % Temperature in Kelvin
kT = k*T;                   % eV
lt = [10^(-9), 10^(-8), 10^(-7), 10^(-6)];    % s estimate a hole lifetime 
mup = 100*(1/100)^2;        % cm^2/Vs*(1m/100cm)^2 
m0 = 5.68562975e-12;        % eV*s^2/m^2 electron mass
m = 0.5*m0;                 % Assuming an effective mass (Look for value)
Dp = kT*mup;                % m^2/s^2
%Diffusion Length (Assuming that I need  a step size smaller than the
%diffusion length)
L_Diff = sqrt(Dp*lt);

%%
%initial guess for carrier density is linear and only considering holes
Ef = 0.5;                   % eV setting a fermi level rel to valence
V_x1 = 0;
V_xN = 1; 

%The density of holes for left and right contact based on parabolic 
%assumption (+degeneracy of K & K') + MaxwellBoltzmann Approx of
%Fermi-level
p0 = 2*(m*kT/(2*pi*hbar^2))^(3/2)*exp(-(Ef)/kT);
pl =  p0;
pr = p0; 
%conductivity should be
sigma = p0*q0*mup;
res = 1/sigma; %resistivity
%R = V/I = res * L /A 
%we are dealing with J in A/m^2
% I =  J*A 
% R = V/(J*A) = res * L / A
% res *L = V/J
%Use current code with no recombination so J is constant

%% %choice of variables
nx = 1000;              % steps in x direction
xL = 5e-6;             
x = linspace(0,xL, nx);
dx = x(2)-x(1);

%Debye length: (represents largest mesh size - see Computation Electronics
%by Vasileska and Goodnick) 
Ld = sqrt(eps*(kT*q0)/(p0*q0^2));  % meters (kT*q0 is in Joules )

p = p0*ones(1,nx);            % the free carrier density
%p = (1e8*x-50).^2*10*(1e2)^3 +p0; %weird initial guess 
%rho = 100*exp(-((x-mean(x))/(100*dx)).^2);
rho = q0*p; 
rho(1) = 0;         %Defined in order to implement boundary conditions
rho(nx) = 0;        %Defined in order to implement boundary conditions

%%
maxIter =1e3;
V = Pois_1D(x, rho, V_x1, V_xN);

for nt = 1:length(lt)
for i=1:maxIter
    Vnew = Pois_1D(x, rho, V_x1, V_xN);
    %rarely vary by greater than 10^14
    %all(abs(V-Vnew)<1e-14)
    %pause(1)
    V = Vnew;
    %p0 = 2.*(m*kT/(2*pi*hbar^2))^(3/2).*exp(-(Ef-V)./kT);
    [J, pnew] = current1D_lifetimetest(x, V, pl, pr, lt(nt));
    
    if(all(abs(pnew-p')< p0/1e10))
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
color = nt/length(lt);
figure(1)
plot(x*1e6, V,'Color', [0, 0, color], 'linewidth', 2)
xlabel('x-axis (\mum)')
ylabel('Potential (V)') % right y-axis
set(gca, 'fontsize', 14);
hold on

figure(2) 
semilogy(x*1e6, p,'Color', [0, 0, color], 'linewidth', 2);
xlabel('x-axis (\mum)')
ylabel('Hole concentration /m^3')
set(gca, 'fontsize', 14);
hold on
pause(1)
end

%%


%%
%{
% Actually looking at the current density and resistivity
%plot current density to check
figure(3)
plot(x(1:end-1)*1e6, J, 'linewidth', 2);
xlabel('x-axis (\mum)')
ylabel('Current Density A/m^2')

%Okay, so the current looks mostly constant which is expected, let's check
%a bit more carefully
Javg = mean(J);
figure(4)
plot(x(1:end-1)*1e6, J-Javg, 'linewidth', 2);
xlabel('x-axis (\mum)')
ylabel('Deviation from average current density A/m^2')
%Deviation smaller than 10^-12
resCode = abs ( V_xN/ (Javg*xL) ) ;
res;
res_pavg = 1/(mean(p)*q0*mup);
Deviation = abs(res-resCode);
%}

