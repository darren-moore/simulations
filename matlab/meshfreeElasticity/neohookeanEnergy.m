function [energy, f_int] = neohookeanEnergy(u, fieldNodes, gaussNodes, hRad, gaussSquareRad, mu, lam, shapeFunG, shapeDxG, shapeDyG, hRadCell)
% neohookeanEnergy Computes neohookean energy and energy derivative(internal forces) for a given displacement u.

numNodes = size(fieldNodes, 1);
numGauss = size(gaussNodes, 1);

energy = 0;
f_int = zeros(numNodes*2,1);

u0 = zeros(numNodes,2);
u0(:,1) = u(1:2:end);
u0(:,2) = u(2:2:end);
refConfig = fieldNodes;

% for each background cell
for i=1:numGauss
	quadEvalPts = [-sqrt(3)/3 sqrt(3)/3; sqrt(3)/3 sqrt(3)/3; sqrt(3)/3 -sqrt(3)/3; -sqrt(3)/3 -sqrt(3)/3];
	quadEvalPts = quadEvalPts * gaussSquareRad; % scale to square space
	quadEvalPts = repmat(gaussNodes(i,:),4,1) + quadEvalPts; % translate to the centre of cell
	
	shapeDx = shapeDxG{i};
	shapeDy = shapeDyG{i};
	
	for quadPt=1:size(quadEvalPts,1)
		theRads = hRadCell{i};
		quadEvalNeighbors = rangesearch(refConfig,quadEvalPts(quadPt,:),theRads(quadPt));
		toCheck = quadEvalNeighbors{1};
		toCheck = sort(toCheck);
		
		die_W = zeros(2);
		
		% compute disp. grad. at the quad point
		for beta=1:numel(toCheck)
			k = toCheck(beta);
			shapeD = [shapeDx(quadPt,k) shapeDy(quadPt,k)];
			die_W = die_W + u([2*k-1 2*k]) * shapeD;
		end
		
		% Add and eye and we get the defo. grad.
		F = die_W + eye(2);
		F_inv = inv(F);
		
		% Some invariants
		I_1 = trace(F'*F);
		I_3 = det(F)^2;

		% first P-K stress.
		PK_F = mu*(F - F_inv') + (lam/2)*log(I_3)*F_inv';
		
		% Sum internal forces
		for alpha=1:numel(toCheck)
			j = toCheck(alpha);
			f_int([2*j - 1, 2*j]) = f_int([2*j - 1, 2*j]) + (gaussSquareRad*2).^2 * PK_F * [shapeDx(quadPt,j) ; shapeDy(quadPt,j)];
		end
		
		% Compute neohookean energy...
		energy_piece = (lam/2)*(I_1 - log(I_3) - 2) + (mu/8)*log(I_3)^2;
		
		% ...multiply by volume and add.
		energy = energy + (gaussSquareRad*2).^2 * energy_piece;
	end
end


% Finally, sum the nodal energies.
energy = sum(energy);

% These were used to visualize the optimization between iterations.
% I enjoy seeing fminunc in action on static simulations especially.
% Feel free to uncomment to see this.

% uMag = sqrt(u0(:,1).^2 + u0(:,2).^2);
% fieldNodes = refConfig + u0;
% cla
% scatter(fieldNodes(:,1), fieldNodes(:,2), 15000/numNodes, uMag,'filled');
% drawnow

end

