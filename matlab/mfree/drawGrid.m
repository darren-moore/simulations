function drawGrid( nodes, diam )
%drawGrid Draws squares of given diameter around each given node
	hold on
	for i=1:size(nodes,1)
		c = nodes(i,:);
		r = diam/2;
		rectangle('Position', [c(1)-r c(2)-r diam diam]);
	end
end

