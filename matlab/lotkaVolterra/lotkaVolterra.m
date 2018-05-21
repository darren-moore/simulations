% Lotka-Volterra predator-prey equations
% https://en.wikipedia.org/wiki/Lotka–Volterra_equations

% Parameters
alpha = 2/3;
beta = 4/3;
gamma = 1;
delta = 1;

% initial conditions, where prey is x(1) and predators is x(2)
x = [1 .7];

dt = .01;

cla
hold on
preyGraph = animatedline('Color', 'g', 'LineWidth',2);
predGraph = animatedline('Color', 'r', 'LineWidth',2);
MAXTIME = 1000;
axis([0,MAXTIME, 0,inf]);

% Comment in and out to experiment with different integrators.
for t=1:MAXTIME
% 	x = forwardEuler(x, dt, alpha, beta, gamma, delta);
% 	x = backwardEuler(x, dt, alpha, beta, gamma, delta);
	x = RK(x, dt, alpha, beta, gamma, delta);
	
	addpoints(preyGraph,t,x(1));
	addpoints(predGraph,t,x(2));
	drawnow
end

% The Lotka-Volterra system.
% Note: the derivative doesn't depend on time!
function dxdt = LV(x,alpha, beta, gamma, delta)
	dxdt(1) = alpha*x(1) - beta*x(1)*x(2);
	dxdt(2) = delta*x(1)*x(2) - gamma*x(2);
end

% Lotka-Volterra Jacobian
function J = DLV(x, alpha, beta, gamma, delta)
	J(1,1) = alpha - beta*x(2);
	J(1,2) = -beta*x(1);
	J(2,1) = delta*x(2);
	J(2,2) = delta*x(1) - gamma;
end


function x_np1 = forwardEuler(x_n, dt, alpha, beta, gamma, delta)
	dxdt_np1 = LV(x_n, alpha, beta, gamma, delta);
	x_np1 = x_n + dt*dxdt_np1;
end


function x_np1 = backwardEuler(x_n, dt, alpha, beta, gamma, delta)
	% Have a system of non-linear ODE.
	% So, restate as a root-finding problem and use Newton's method.
	% Just use 5 iterations
	x = x_n';
	for i=1:5
		J = DLV(x, alpha, beta, gamma, delta);
		fx = LV(x, alpha, beta, gamma, delta);
		Fx = x - x_n - dt*fx';		
		DF = eye(2) - dt*J;
		x = x - inv(DF)*Fx;
	end
	
	x_np1 = x;
end

% 4th order Runge-Kutta: https://en.wikipedia.org/wiki/Runge–Kutta_methods
% Since LV system doesn't vary in time, no need to pass it...
function x_np1 = RK(x_n, dt, alpha, beta, gamma, delta)
	k1 = dt*LV(x_n, alpha, beta, gamma, delta);
	k2 = dt*LV(x_n + k1/2, alpha, beta, gamma, delta);
	k3 = dt*LV(x_n + k2/2, alpha, beta, gamma, delta);
	k4 = dt*LV(x_n + k3, alpha, beta, gamma, delta);
	
	x_np1 = x_n + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
end

