% ELEC 4700 Assignment 2 Part 2 MESH DENSITY
% Andrew Branicki 100973961
% February 24, 2019
clear;
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
    
% make a new variable that will change the mesh density and LOOP IT!
for mesh_density = 1:10
    
    % Declare boundaries of the canvas
    % For plots, use a ratio of 3/2 for L/W
    % Similar to EIPA from February 6, 2019
    % SMALLER DIMENSIONS TO MAKE THIS CODE RUNNABLE...
    L = 15; W = 10; dx = 1; dy = 1;
    G = sparse(L*W, W*L);
    V = zeros(L*W,1);
    cond = zeros(L,W);

    % Conductivity values
    sigma_o = 1;        % outside boxes 
    sigma_in = 1e-2;    % inside boxes

    % set up the regions for the boxes
    region_x = L/3;
    region_y = W/3;
    middle_x = L/2;
    middle_y = W/2;

    % This is how many points we will actually need since our density changes
    points_numX = L*mesh_density;
    points_numY = W*mesh_density;
    map = @(i,j) j + (i - 1)*points_numY;
    mesh_density_sqr = mesh_density^2;

    % Start setting up boundary conditions
    % Need to set up the conductivity matrix
    % So that the conductivity is sigma_in in the areas
    % where the block is present
    % X direction
    for i=1:points_numX
        % Y direction
        for j=1:points_numY
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
                % Now we have a new spacing, similar to the previous PA
                G(n,n) = 1/mesh_density_sqr;
                V(n) = 1;
            % nx = L (in this case, V = 0)
            elseif i == points_numX
                cond(i,j) = sigma_o;
                G(n,:) = 0;
                G(n,n) = 1/mesh_density_sqr;
                V(n) = 0;
            % Y edges of canvas. Need to check if the specific X coordinate at
            % the top and bottom edges are inside the box or not
            elseif (j == 1 || j == points_numY)
                G(n,:) = 0;
                % IN THIS CASE WE ARE INSIDE THE BOX!!!
                if (i/mesh_density > (middle_x - region_x/2) && i/mesh_density < (middle_x + region_x/2))
                    cond(i,j) = sigma_in;
                    G(n,n) = -3*sigma_in/mesh_density_sqr;
                    G(n,nxm) = sigma_in/mesh_density_sqr;
                    G(n,nxp) = sigma_in/mesh_density_sqr;
                    G(n,nyp) = sigma_in/mesh_density_sqr;
                % IN THIS CASE WE ARE NOT INSIDE THE BOX!!!
                else
                    cond(i,j) = sigma_o;
                    G(n,n) = -3*sigma_o/mesh_density_sqr;
                    G(n,nxm) = sigma_o/mesh_density_sqr;
                    G(n,nxp) = sigma_o/mesh_density_sqr;
                    G(n,nyp) = sigma_o/mesh_density_sqr;
                end
            % THIS IS THE GENERAL CASE, WHERE THE POINT ISN'T EDGY.    
            else
                G(n,:) = 0;
                G(n,n) = -4/mesh_density_sqr;
                % Check if inside the box again (except here we check x and y)
                if ((i/mesh_density > (middle_x - region_x/2) && i/mesh_density < (middle_x + region_x/2)) && (j/mesh_density > (middle_y + region_y/2) || j/mesh_density < (middle_y - region_y/2)))  
                    cond(i,j) = sigma_in;
                    G(n,nxm) = sigma_in/mesh_density_sqr;
                    G(n,nxp) = sigma_in/mesh_density_sqr;
                    G(n,nym) = sigma_in/mesh_density_sqr;
                    G(n,nyp) = sigma_in/mesh_density_sqr;
                else % NOT INSIDE BOX!
                    cond(i,j) = sigma_o;
                    G(n,nxm) = sigma_o/mesh_density_sqr;
                    G(n,nxp) = sigma_o/mesh_density_sqr;
                    G(n,nym) = sigma_o/mesh_density_sqr;
                    G(n,nyp) = sigma_o/mesh_density_sqr;
                end
            end
        end
    end

    % Solve for F using G\V = F
    F = G\V;

    % Set up a surf plot
    surfs_up = zeros(points_numX,points_numY);
    for i = 1:points_numX
        for j = 1:points_numY
            n = map(i,j);
            surfs_up(i,j) = F(n);
        end
    end
    
    J = cond.*gradient(surfs_up); % gradient of surfs up is E field
    % J is now an ARRAY WITH A WHOLE BUNCH OF POINTS!
    % since we can only plot ONE current point for this specific
    % mesh density, we're gonna have to take the average.
    J_average = sum(sum(J))/(points_numX*points_numY);
    % And since J is current density, we will need to divide by area
    % to get the actual current.
    CURRENT(mesh_density) = J_average/(L*W);
end
% Once we leave this for loop, we now have 10 CURRENT VALUES
% for 10 different mesh densities! Woo!

figure(8)
plot(1:10,CURRENT)
title('Effect of Mesh Density on Current')
xlabel('Number of Points per Mesh')
ylabel('Average Current (A)')
grid on