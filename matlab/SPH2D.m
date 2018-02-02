
% A minimal 2D SPH simulation in MATLAB.
% Follows the method of "Particle-Based Fluid Simulation for Interactive Applications" [Müller et al. 2003]

clear;
set(gcf,'Renderer','OpenGL');


% Simulation settings.
% It's tough to get something stable and looks right.
nodeResolution = 22;	% Determines how many particles to use.
dt = 1e-3;
gravity = -1e4;
mass = 10;
hRad = 3;			% The radius of influence.
gas_k = 1500;		% Gas constant
rho_rest = 1000;	% Rest density
mu = 100;			% Viscosity
damper = .5;		% Damping for wall collisions.

numNodes = nodeResolution^2;

% Kernels and kernel derivative pieces stored as constants for efficiency.
poly6 = 315/(65*pi*hRad^9);
Dspike = -45/(pi*hRad^6);
DDvisc = 45/(pi*hRad^6);

% Changing this requires finding all new parameters. Best to leave it I think.
windowSize = 100;

% Set up particles in space
[X, Y] = meshgrid(linspace(0,1,nodeResolution), linspace(0,1,nodeResolution));
fieldNodes = [X(:) Y(:)];
fieldNodes = fieldNodes + rand(numNodes,2)*.01;	% Wiggle the particles a bit
fieldNodes = fieldNodes .* .7*windowSize;
fieldNodes = fieldNodes+.15*windowSize;

% initial velocity set to 0
fieldNodesV = zeros(numNodes,2);

frames = 5000;
for t=1:frames
	rho = zeros(numNodes, 1);
	pressure = zeros(numNodes, 1);
	forces = zeros(numNodes, 2);
	
	neighbors = rangesearch(fieldNodes, fieldNodes, hRad);
	
	% This loop computes density and pressure for each node
	for i=1:numNodes
		% get the local neighborhood from the neighbors data structure
		hood = neighbors(i); hood = hood{1};
				
		% Difference in position between r_j and r.
		difference = [fieldNodes(hood,1)-fieldNodes(i,1) fieldNodes(hood,2)-fieldNodes(i,2)];

		% in general we want to sum mass_j, but in this sim, mass is const.
		rho(i) = mass * poly6*sum((hRad^2 - sum(difference.^2,2)).^3);
		
		% pressure = c^2 * (density_computed - density_reference)
		pressure(i) = gas_k*(rho(i) - rho_rest);
	end
	
	% This loop computes the internal forces for each node, using the densities and pressures we just computed.
	for i=1:numNodes
		% get the local neighborhood from the neighbors data structure
		hood = neighbors(i); hood = hood{1}; hood = hood(hood ~= i);
		if(numel(hood) == 0)
			continue;
		end

		difference = [fieldNodes(hood,1)-fieldNodes(i,1) fieldNodes(hood,2)-fieldNodes(i,2)];
		dists = sqrt(sum(difference.^2,2));
		
		fPressure = -mass*(difference./dists) .* ((pressure(i) + pressure(hood)) ./ (rho(hood)*2)) .* (Dspike*(hRad-dists).^2);
		fPressure = sum(fPressure,1);
		
		fViscosity = mu*mass*((fieldNodesV(hood,:) - fieldNodesV(i,:))./rho(hood)) .* (DDvisc*(hRad-dists));
		fViscosity = sum(fViscosity,1);
		
		forces(i,:) = fPressure + fViscosity;	
	end
	
	forces(:,2) = forces(:,2) + gravity;

	accel = forces ./ rho;

	% Integrate with euler integration
	fieldNodesV = fieldNodesV + dt*accel; % velocity = \int{accel}
	oldFieldNodes = fieldNodes;
	fieldNodes = fieldNodes + dt*fieldNodesV; % position = \int{velocity}

	% Simple wall collisions.
	% Collect nodes that are hitting the sides or top/bottom
	toFlipX = fieldNodes(:,1) < 0 | fieldNodes(:,1) > windowSize;
	toFlipY = fieldNodes(:,2) < 0 | fieldNodes(:,2) > windowSize;
	
	% Reflect velocity component that is perpendicular to wall:	
	fieldNodesV(toFlipX,1) = -damper*fieldNodesV(toFlipX,1);
	fieldNodesV(toFlipY,2) = -damper*fieldNodesV(toFlipY,2);
	
	% Then push out of wall:
	fieldNodes(fieldNodes(:,1) < 0, 1) = .1;
	fieldNodes(fieldNodes(:,1) > windowSize, 1) = windowSize-.1;
	fieldNodes(fieldNodes(:,2) < 0, 2) = .1;
	fieldNodes(fieldNodes(:,2) > windowSize, 2) = windowSize-.1;
		
	% Display
	hold on
	cla
	axis([0 windowSize 0 windowSize])
	scatter(fieldNodes(:,1), fieldNodes(:,2),[],'filled');
% 	quiver(fieldNodes(:,1), fieldNodes(:,2), fieldNodesV(:,1),fieldNodesV(:,2));
	drawnow
		
end

