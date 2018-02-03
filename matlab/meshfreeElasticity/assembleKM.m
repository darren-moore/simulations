function [K, M] = assembleKM( fieldNodes, gaussNodes, hRad, E, nu, gaussSquareRad, shapeCell, shapeDxCell, shapeDyCell, hRadCell)
% assembleKM Assembles stiffness and mass matrices.

	numNodes = size(fieldNodes, 1);
	numGauss = size(gaussNodes, 1);
	
	% Material matrix
	D = E/(1-nu^2) * [1 nu 0 ; nu 1 0 ; 0 0 (1-nu)/2];
	
	K = zeros(numNodes*2, numNodes*2);
	M = zeros(numNodes*2, numNodes*2);
	
	refConfig = fieldNodes;
	
	% for each background cell
	for i=1:numGauss
		quadEvalPts = [-sqrt(3)/3 sqrt(3)/3; sqrt(3)/3 sqrt(3)/3; sqrt(3)/3 -sqrt(3)/3; -sqrt(3)/3 -sqrt(3)/3];
		quadEvalPts = quadEvalPts * gaussSquareRad; % scale to square space
		quadEvalPts = repmat(gaussNodes(i,:),4,1) + quadEvalPts; % translate to the centre of cell
		numQuadPts = 4;
		
		shapeFunG = shapeCell{i};
		shapeDxG = shapeDxCell{i};
		shapeDyG = shapeDyCell{i};
		hRad_custom =  hRadCell{i};
		
		for quadPt=1:numQuadPts
			quadEvalNeighbors = rangesearch(refConfig,quadEvalPts(quadPt,:),hRad_custom(quadPt));
			toCheck = quadEvalNeighbors{1};

			% assemble the gradient operator, B
			B = zeros(3,2*numNodes);
			for idx=1:numel(toCheck)
				j = toCheck(idx);
				die_W = [shapeDxG(quadPt,j) shapeDyG(quadPt,j)];
				B(1,2*j-1)	= die_W(1);
				B(2,2*j)	= die_W(2);
				B(3,2*j-1)	= die_W(2);
				B(3,2*j)	= die_W(1);
			end

			% assemble K_ij and M_ij. Weights for 4 gauss points are all 1.
			for alpha=1:numel(toCheck)
				j = toCheck(alpha);

				for beta=1:numel(toCheck)
					k = toCheck(beta);

					% Slice the differential operator
					B_j = B(:,[2*j-1 2*j]);
					B_k = B(:,[2*k-1 2*k]);

					K_ij_piece = B_j'*D*B_k * (gaussSquareRad*2).^2; % Jacobian is the square length
					K([2*j-1 2*j], [2*k-1 2*k]) = K([2*j-1 2*j], [2*k-1 2*k]) + K_ij_piece;

					phi_j = shapeFunG(quadPt,j);
					phi_k = shapeFunG(quadPt,k);
					M_ij_piece = eye(2) * phi_j*phi_k*(gaussSquareRad*2).^2;
					M([2*j-1 2*j], [2*k-1 2*k]) = M([2*j-1 2*j], [2*k-1 2*k]) + M_ij_piece;
				end
			end
		end
	end
end

