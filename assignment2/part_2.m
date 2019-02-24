% ELEC 4700 Assignment 2 Part 2
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
L = 75; W = 50; dx = 1; dy = 1;
G = sparse(L*W, W*L);
V = zeros(L*W,1);
cond = zeros(L,W);
map = @(i,j) j + (i - 1)*W;

% Conductivity values
sigma_o = 1;        % outside boxes 
sigma_in = 1e-2;    % inside boxes

% set up the regions for the boxes
region_x = L/3;
region_y = W/3;
middle_x = L/2;
middle_y = W/2;

% Start setting up boundary conditions
% Need to set up the conductivity matrix
% So that the conductivity is sigma_in in the areas
% where the block is present
% X direction
for i=1:L
    % Y direction
    for j=1:W
        n = map(i,j);
        nxm = map(i-1,j);
        nxp = map(i+1,j);
        nym = map(i,j-1);
        nyp = map(i,j+1);
        
        % nx = 0 (in this case, V = Vo)
        % we also don't have boxes on the left and right edges so don't
        % need to check anything :)
        if i == 1
            cond(i,j) = sigma_o;
            G(n,:) = 0;
            G(n,n) = 1;
            V(n) = 1;
        % nx = L (in this case, V = 0)
        elseif i == L
            cond(i,j) = sigma_o;
            G(n,:) = 0;
            G(n,n) = 1;
            V(n) = 0;
        % Y edges of canvas. Need to check if the specific X coordinate at
        % the top and bottom edges are inside the box or not
        elseif (j == 1 || j == W)
            G(n,:) = 0;
            % IN THIS CASE WE ARE INSIDE THE BOX!!!
            if (i > (middle_x - region_x/2) && i < (middle_x + region_x/2))
                cond(i,j) = sigma_in;
                G(n,n) = -3*sigma_in;
                G(n,nxm) = sigma_in;
                G(n,nxp) = sigma_in;
                G(n,nyp) = sigma_in;
            % IN THIS CASE WE ARE NOT INSIDE THE BOX!!!
            else
                cond(i,j) = sigma_o;
                G(n,n) = -3*sigma_o;
                G(n,nxm) = sigma_o;
                G(n,nxp) = sigma_o;
                G(n,nyp) = sigma_o;
            end
        % THIS IS THE GENERAL CASE, WHERE THE POINT ISN'T EDGY.    
        else
            G(n,:) = 0;
            G(n,n) = -4;
            % Check if inside the box again (except here we check x and y)
            if ((i > (middle_x - region_x/2) && i < (middle_x + region_x/2)) && (j > (middle_y + region_y/2) || j < (middle_y - region_y/2)))  
                cond(i,j) = sigma_in;
                G(n,nxm) = sigma_in;
                G(n,nxp) = sigma_in;
                G(n,nym) = sigma_in;
                G(n,nyp) = sigma_in;
            else % NOT INSIDE BOX!
                cond(i,j) = sigma_o;
                G(n,nxm) = sigma_o;
                G(n,nxp) = sigma_o;
                G(n,nym) = sigma_o;
                G(n,nyp) = sigma_o;
            end
        end
    end
end

% Solve for F using G\V = F
F = G\V;

% Set up a surf plot
surfs_up = zeros(L,W);
for i = 1:L
    for j = 1:W
        n = map(i,j);
        surfs_up(i,j) = F(n);
    end
end

% TIME TO PLOT COOL THINGS!
[EX, EY] = gradient(surfs_up);
J = cond.*gradient(surfs_up);

% CONDUCTIVITY PLOT
figure(3)
surf(cond)
% It's almost spring so this feels appropriate
colormap summer
colorbar
title('Conductivity Plot')
xlabel('Width')
ylabel('Length')
zlabel('Sigma')

% VOLTAGE PLOT
figure(4)
surf(surfs_up)
colormap summer
colorbar
title('Voltage Plot')
xlabel('Width')
ylabel('Length')
zlabel('Voltage')

% EX PLOT
figure(5)
surf(EX)
colormap summer
colorbar
title('X Direction Electric Field Plot')
xlabel('Width')
ylabel('Length')
zlabel('Electric Field in the X Direction')

% EY (lmao) PLOT
figure(6)
surf(EY)
colormap summer
colorbar
title('Y Direction Electric Field Plot')
xlabel('Width')
ylabel('Length')
zlabel('Electric Field in the Y Direction')

% CURRENT DENSITY PLOT
figure(7)
surf(J)
colormap summer
colorbar
title('Current Density Plot')
xlabel('Width')
ylabel('Length')
zlabel('Current Density')