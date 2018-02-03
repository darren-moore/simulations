
% A 2D meshfree non-linear elastisity simulation for static and dynamic problems.

clear;
rng(2);

figure(1);
hold on
cla
axis square
axis([0 .8 0 .8]);

%% Simulation settings
dt = 1e-3;
gravity = -10;

E = 1e5;	% Younge's Modulus.
nu = .3;	% Poisson's ratio.

% Feel free to set this to 0 for a linear solve.
max_newton_iterations = 3;

% We simulate a bar pinned to the wall.
% The following controls how many particles the bar is on the Y-axis
pointDensity = 5;

% Setting the integration grid density equal to pointDensity should give nice results.
% Values that are less than perhaps (pointDensity/2) will cause the simulation to fail.
integrationGridDensity = pointDensity;

dynamic = true;	% dynamic vs static simulation
randomShift = .01; % how much random shift to apply to points for irregular distribution

% The radius of influence for a particle.
% This is computed automatically in shapeCompute, but this value will be 
% used if you uncomment a line in that function.
hRad = .1;

optionsQuad = optimoptions(@quadprog, 'Display', 'off');
warning('off','optim:quadprog:HessianNotSym'); % turned off since error was 1e-15. *just* above eps

optionsUnc = optimoptions(@fminunc, 'Display', 'iter','OutputFcn', @outfun,...
	'MaxIterations', max_newton_iterations ,...
	'SpecifyObjectiveGradient', true, 'Algorithm', 'quasi-newton');


% Lame coefficients. A function of E and nu.
mu = E / (2*(1+nu));
lam = (E*nu) / ((1+nu)*(1-2*nu));


%% Set up the simulation

nodesX = linspace(0, .4, pointDensity*4);
nodesY = linspace(.5, .6, pointDensity);
[X, Y] = meshgrid(nodesX, nodesY);
fieldNodes = [X(:) Y(:)];
numNodes = length(fieldNodes);
boundNodes0 = (fieldNodes(:,1) == 0);

% add some random noise
shifter = (2*rand(numNodes,2)-1) * randomShift;
fieldNodes = fieldNodes+shifter;



% create background integration grid
integY = linspace(.5, .6, integrationGridDensity);
gaussSquareRad = (integY(2) - integY(1))/2;
integX = linspace(0, .4, ceil(.4/(gaussSquareRad * 2)) + 1); % overshoot
currentLen = .4+2*gaussSquareRad;
actualLen = (ceil(.4/(gaussSquareRad * 2)) + 1) * (gaussSquareRad * 2);
integX = integX * actualLen/currentLen - (actualLen - currentLen)/2; % and adjust
[X, Y] = meshgrid(integX, integY);
gaussNodes = [X(:) Y(:)];	% By "gaussNode", I mean the centre of an integration grid cell.
numGauss = length(gaussNodes);



% We should remember the initial configuration.
initConfig = fieldNodes;

% We will store previous displacements for time integration.
disp_0 = zeros(numNodes*2,1);
disp_00 = zeros(numNodes*2,1);

% We can draw our initial configuration and integration grid before starting the simulation:
drawGrid(gaussNodes, gaussSquareRad*2);
scatter(fieldNodes(:,1), fieldNodes(:,2),5000/numNodes);
drawnow


% Boundary conditions: "Pin" nodes on the very left to a wall.
boundNodes = zeros(numNodes*2,1);
boundNodes(1:2:end) = boundNodes0;
boundNodes(2:2:end) = boundNodes0;

% The degrees of freedom are those nodes that *aren't* pinned.
dof = ~logical(boundNodes0);
dof2 = ~logical(boundNodes);

% Setup external forces. Just gravity for now.
extForces = zeros(numNodes*2,1);
extForces(1:2:end) = ones(numNodes,1) * 0;
extForces(2:2:end) = ones(numNodes,1) * gravity;

% Compute shape functions.
% We store the shape functions in a cell data structure.
% One cell per integration grid cell.
shapeCell = {numGauss};
shapeDxCell = {numGauss};
shapeDyCell = {numGauss};
hRadCell = {numGauss};
for i=1:numGauss
	% Using 4 point gauss quadrature.
	quadEvalPts = [-sqrt(3)/3 sqrt(3)/3; sqrt(3)/3 sqrt(3)/3; sqrt(3)/3 -sqrt(3)/3; -sqrt(3)/3 -sqrt(3)/3];
	quadEvalPts = quadEvalPts * gaussSquareRad; % scale to square space
	quadEvalPts = repmat(gaussNodes(i,:),4,1) + quadEvalPts; % translate to the centre of cell
	[shapeCell{i}, shapeDxCell{i}, shapeDyCell{i}, hRadCell{i}] = shapeCompute(quadEvalPts, fieldNodes, hRad);
end

% Compute stiffness and mass matrices
[K, M] = assembleKM( fieldNodes, gaussNodes, hRad, E, nu, gaussSquareRad, shapeCell, shapeDxCell, shapeDyCell, hRadCell);

%% Dynamics
if(dynamic)
	for t=1:1000
		fun = @(x) neohookeanEnergy(x,initConfig,gaussNodes, hRad,gaussSquareRad,mu,lam, shapeCell, shapeDxCell, shapeDyCell, hRadCell);
		fUnc = @(x) energyWrapper(fun,x,M, extForces, dt, dof, dof2, disp_0, disp_00);

		% Implicit euler for the linear solve.
		A = M + (dt^2 * K);
		B =  M*(2*disp_0 - disp_00) + (extForces)*(dt^2);
		
		disp_linearGuess = quadprog(A(dof2,dof2), -B(dof2), [],[],[],[],[],[],[], optionsQuad);
		if(max_newton_iterations > 0)
			disp_final = fminunc(fUnc, disp_linearGuess, optionsUnc);
		else
			disp_final = disp_linearGuess;
		end

		% "upgrade" displacement x with boundary condition nodes of displacement zero.
		disp_1 = zeros(numNodes*2,1);
		disp_1(dof2) = disp_final;

		% Un-interleave the displacement
		disp_columns = zeros(numNodes,2);
		disp_columns(:,1) = disp_1(1:2:end);
		disp_columns(:,2) = disp_1(2:2:end);

		fieldNodes = initConfig + disp_columns;
		disp_magnitude = sqrt(disp_columns(:,1).^2 + disp_columns(:,2).^2);
		
		% Display, with colors representing displacement.
		cla
		scatter(fieldNodes(:,1),fieldNodes(:,2),15000/numNodes,disp_magnitude,'filled');
		drawnow
	
		disp_00 = disp_0;
		disp_0 = disp_1;
	end
end

%% Statics
if(~dynamic)

	fun = @(x) neohookeanEnergy(x,fieldNodes,gaussNodes, hRad,gaussSquareRad,mu,lam, shapeCell, shapeDxCell, shapeDyCell, hRadCell);
	fUnc = @(x) staticEnergyWrapper(fun,x,M,extForces, dt, dof, dof2);

	% Compute displacements with boundary conditions
	disp_linearGuess = quadprog(K(dof2,dof2), -extForces(dof2), [],[],[],[],[],[],[], optionsQuad);
	if(max_newton_iterations > 0)
		disp_final = fminunc(fUnc, disp_linearGuess, optionsUnc);
	else
		disp_final = disp_linearGuess;
	end
	
	% Re-add zeros to displacement
	disp_full = zeros(numNodes*2,1);
	disp_full(dof2) = disp_final;

	% Un-interleave the displacement
	disp_columns = zeros(numNodes,2);
	disp_columns(:,1) = disp_full(1:2:end);
	disp_columns(:,2) = disp_full(2:2:end);

	fieldNodes = initConfig + disp_columns;
	disp_magnitude = sqrt(disp_columns(:,1).^2 + disp_columns(:,2).^2);

	% Display, with colors representing displacement.
	cla
	scatter(fieldNodes(:,1),fieldNodes(:,2),15000/numNodes,disp_magnitude,'filled');
	drawnow
end


% An output function for fminunc. Not utilized right now.
% I used to use it to visualize the optimization between iterations hence
% the commented out "drawnow"
function stop = outfun(x,optimValues,state)
% 	drawnow
	stop = false;
end
