% ELEC 4700 Assignment 2 Part 1
% Andrew Branicki 100973961
% February 24, 2019

global C
    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;                      % metres (32.1740 ft) per s²
    C.m_n = 0.26*C.m_0;                 % effective mass of electrons
    
% Declare boundaries of the canvas
% For plots, use a ratio of 3/2 for L/W
% Similar to EIPA from February 6, 2019
nx = 75; ny = 50; dx = 1; dy = 1;
G = sparse(nx*ny, ny*nx);
V = zeros(nx*ny,1);
map = @(i,j) j + (i - 1)*ny;

%% Part a --------------------------------------------------------

% Start setting up boundary conditions
% X direction
for i=1:nx
    % Y direction
    for j=1:ny
        n = map(i,j);
        nxm = map(i-1,j);
        nxp = map(i+1,j);
        nym = map(i,j-1);
        nyp = map(i,j+1);
        
        % nx = 0 (in this case, V = Vo)
        if i == 1
            G(n,:) = 0;
            G(n,n) = 1;
            V(n) = 1;
        % nx = L (in this case, V = 0)
        elseif i == nx
            G(n,:) = 0;
            G(n,n) = 1;
            V(n) = 0;
        elseif (j == 1 || j == ny)
            G(n,:) = 0;
            G(n,n) = -3;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nyp) = 1;
        else
            G(n,:) = 0;
            G(n,n) = -4;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1;
        end
    end
end

% Solve for F using G\V = F
F = G\V;

% Set up a surf plot
surfs_up = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = map(i,j);
        surfs_up(i,j) = F(n);
    end
end

%Plot
figure(1)
surf(surfs_up)
% It's almost spring so this feels appropriate
colormap summer
shading interp
colorbar

title('Electrostatic Potential in Rectangular Region L/W = 3/2; V = V0 at x = 0 and V = 0 at x = L')
xlabel('Width')
ylabel('Length')
zlabel('Voltage')

%% Part b --------------------------------------------------------

% Start setting up boundary conditions
% X direction
for i=1:nx
    % Y direction
    for j=1:ny
        n = map(i,j);
        nxm = map(i-1,j);
        nxp = map(i+1,j);
        nym = map(i,j-1);
        nyp = map(i,j+1);
        
        % nx = 0 (in this case, V = Vo)
        if i == 1
            G(n,:) = 0;
            G(n,n) = 1;
            V(n) = 1;
        % nx = L (in this case, V = 0)
        elseif i == nx
            G(n,:) = 0;
            G(n,n) = 1;
            V(n) = 1;
        elseif j == 1
            G(n,:) = 0;
            G(n,n) = 1;
            V(n) = 0;
        elseif j == ny
            G(n,:) = 0;
            G(n,n) = 1;
            V(n) = 0;
        else
            G(n,:) = 0;
            G(n,n) = -4;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1;
        end
    end
end

% Solve for F using G\V = F
F = G\V;

% Set up a surf plot
surfs_up = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = map(i,j);
        surfs_up(i,j) = F(n);
    end
end

%Plot
figure(2)
surf(surfs_up)
% It's almost spring so this feels appropriate
colormap summer
shading flat
colorbar

title('Electrostatic Potential in Rectangular Region L/W = 3/2; V = V0 at x = 0 and V = 0 at x = L')
xlabel('Width')
ylabel('Length')
zlabel('Voltage')