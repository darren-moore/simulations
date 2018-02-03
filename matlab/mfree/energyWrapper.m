function [ energy, energyGrad ] = energyWrapper( energyFunc, x, M, extForces, dt, dof, dof2, disp_0, disp_00 )
% energyWrapper A wrapper for an energy function that allows for reduced degrees of freedom. (dynamics problems)

%	"upgrade" displacement x with boundary condition nodes of displacement zero.
	u = zeros(numel(extForces), 1);
	u(dof2) = x;

	% Compute energy and internal forces
	[energyPre, f_int] = energyFunc(u);
	
	% Implicit euler time integration.
	energy = u'*M*(u - 2*disp_0 + disp_00) - u'*extForces*dt^2 + dt^2*energyPre;
	gradPre = M*(u - 2*disp_0 + disp_00) - extForces*dt^2 + dt^2*(f_int);
	
%	"downgrade" gradient by removing irrelevent degrees of freedom
	energyGrad = gradPre(dof2);
	
end

