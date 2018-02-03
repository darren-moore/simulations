function [ energy, energyGrad ] = staticEnergyWrapper( energyFunc, x, M, extForces, dt, dof, dof2 )
% staticEnergyWrapper A wrapper for an energy function that allows for reduced degrees of freedom. (statics problems)

%	"upgrade" displacement x with boundary condition nodes of displacement zero.
	u = zeros(numel(extForces), 1);
	u(dof2) = x;

	% Compute energy and internal forces
	[energyPre, f_int] = energyFunc(u);
	
	energy = energyPre - (extForces)'*u;
	gradPre = f_int - extForces;

%	"downgrade" gradient by removing irrelevent degrees of freedom
	energyGrad = gradPre(dof2);
	
end

