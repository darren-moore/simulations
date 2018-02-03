function [shapeFun, shapeDx, shapeDy hRad_custom] = shapeCompute(targetNodes, samplePoints, hRad)
% Computes Moving Least Squares shape functions, shape function derivatives, and a suitable neighborhood for particles.

	targetLen = size(targetNodes,1);
	sampleLen = size(samplePoints,1);
	
	% For a point's neighborhood, we ensure it contains at least 'm' points.
	m = 5;
	
	Mdl = KDTreeSearcher(samplePoints);
	hRad_custom = zeros(targetLen,1);
	theHood = cell(targetLen,1);
	for i=1:targetLen
		idx = knnsearch(Mdl,targetNodes(i,:),'K',m); % 10, 1.3; 5, 1.7
		theRad = max(pdist2(samplePoints(idx,:), targetNodes(i,:)));
		hRad_custom(i) = 1.3*theRad;
		
		% If you want to use a provided radius of influence, uncomment.
% 		hRad_custom(i) = hRad;
		
		this = rangesearch(samplePoints,targetNodes(i,:),hRad_custom(i));
		theHood{i} = this{1};
	end
	
	shapeFun = zeros(targetLen,sampleLen);
	M = zeros(3,3,targetLen);
	B = zeros(3,sampleLen,targetLen);
	
	M_inv = zeros(3,3,targetLen);
	M_Dx = zeros(3,3,targetLen);
	M_Dy = zeros(3,3,targetLen);
	
	% For shape derivatives
	shapeDx = zeros(targetLen,sampleLen);
	shapeDy = zeros(targetLen,sampleLen);
	
	% Compute moment matrices
	for i=1:targetLen
		toCheck = theHood{i};
		if(numel(toCheck)==0)
			continue;
		end
		hRad = hRad_custom(i);
		difference = samplePoints(toCheck,:) - targetNodes(i,:);
		
		% Sum
		for k=1:numel(toCheck)
			j = toCheck(k);
			p = [1 samplePoints(j,:)];
			M(:,:,i) = M(:,:,i) + kernelMLS(difference(k,:),hRad,0) * (p'*p);
			B(:,j,i) = kernelMLS(difference(k,:),hRad,0) * p';
			
			kernelD = kernelMLS(difference(k,:),hRad,1);
			M_Dx(:,:,i) = M_Dx(:,:,i) +  kernelD(1)*(p'*p);
			M_Dy(:,:,i) = M_Dy(:,:,i) +  kernelD(2)*(p'*p);
			
		end
		M_inv(:,:,i) = inv(M(:,:,i));
		p = [1 targetNodes(i,:)];
		
		shapeFun(i,:) = p*M_inv(:,:,i)*B(:,:,i);
		
		% Compute shape function derivatives
		for k=1:numel(toCheck)
			j = toCheck(k);
			p = [1 targetNodes(i,:)]';
			w = kernelMLS(difference(k,:),hRad,0);
			dw = kernelMLS(difference(k,:),hRad,1);
			pdx = [0 1 0];
			pdy = [0 0 1];
			
			Mdx = M_Dx(:,:,i);
			Mdy = M_Dy(:,:,i);

			dM_invdx = -M_inv(:,:,i) * Mdx * M_inv(:,:,i);
			dM_invdy = -M_inv(:,:,i) * Mdy * M_inv(:,:,i);
			
			p2 = [1 samplePoints(j,:)]';
			
			shapeDx(i,j) = dw(1)*p'*M_inv(:,:,i)*p2 + w*p'*dM_invdx*p2 + w*pdx*M_inv(:,:,i)*p2;
			shapeDy(i,j) = dw(2)*p'*M_inv(:,:,i)*p2 + w*p'*dM_invdy*p2 + w*pdy*M_inv(:,:,i)*p2;
		end
	end

	% Ensure it's a partition of unity.
	shapeFun = shapeFun ./ sum(shapeFun,2);
	
end
	