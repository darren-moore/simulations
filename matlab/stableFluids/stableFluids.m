% A 2D fluid simulation with periodic boundary conditions.
% Based on Stable Fluids by Jos Stam. http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/ns.pdf

cla
hold on
axis square

% Initialization
global gridres;
gridres = 128;
axis([.5 gridres+.5 .5 gridres+.5])
[ex,why] = meshgrid(1:gridres,1:gridres);	% Used for rendering.

c_temp = 10;
dt = 1/gridres;		% dt is constant for simplicity

% Set initial velocity field.
v_x = zeros(gridres);
v_y = zeros(gridres);

% Set initial temperatures.
q = zeros(gridres);

% An interesting pattern...
q(floor(gridres*1/8):floor(gridres*1/2),floor(gridres*3/8):floor(gridres*5/8)+1) = -1;
q(floor(gridres*1/2)+1:floor(gridres*7/8)+1,floor(gridres*3/8):floor(gridres*5/8)+1) = 1;
q(floor(gridres*1/2)+1:floor(gridres*7/8)+1,floor(gridres*1/8):floor(gridres*2/8)) = -1;
q(floor(gridres*1/8):floor(gridres*1/2),floor(gridres*6/8)+1:floor(gridres*7/8)+1) = 1;

MAXTIME = 1000;
for t=1:MAXTIME
	v_x_old = v_x;
	v_y_old = v_y;
	
	% Add forces due to temperature. Hot rises. Cold sinks.
	v_y = v_y + dt*c_temp*q;
	
	% Advect vector field and heat with semi-lagrangian technique.
	% Integration is explicit euler.
	v_x = advect_vectorized(v_x, v_x_old, v_y_old, dt);
	v_y = advect_vectorized(v_y, v_x_old, v_y_old, dt);
	q = advect_vectorized(q, v_x_old, v_y_old, dt);

	% Move into frequency domain
	fd_x = fftshift(fft2(v_x));
	fd_y = fftshift(fft2(v_y));

	% Project out the divergence.
	[fd_x, fd_y] = project_out_divergence_vectorized(fd_x, fd_y);

	% Return to spacial domain
	v_x = ifft2(ifftshift(fd_x),'symmetric');
	v_y = ifft2(ifftshift(fd_y),'symmetric');
	
	% Visualize pressure and velocity field
	cla;
	imagesc(q');
% 	quiver(ex',why',v_x, v_y, 'AutoScale', 'on', 'AutoScaleFactor', .9, 'Color','black', 'LineWidth', 1, 'ShowArrowHead', 'off');
	drawnow
end

% Advection with semi-lagrangian method. Loopy edition.
function q_new = advect(q, v_x, v_y, dt)
	global gridres;
	q_new = zeros(gridres);
	for i=1:gridres
		for j=1:gridres	
			X_curr = ind2pos(i,j);
			X_prev = X_curr - dt*[v_x(i,j) v_y(i,j)];	% Explicit euler
			q_new(i,j) = interpolate_value(X_prev, q);	% Use bilinear interpolation.
		end
	end
end

% Projects a frequency domain signal into a divergence-free space. Loopy edition.
function [fd_x, fd_y] = project_out_divergence(fd_x, fd_y)
	global gridres;
	for i=1:gridres
		for j=1:gridres
			v_wave = ind2pos(i-.5,j-.5) - [.5 .5];
			nrm = norm([-v_wave(2) v_wave(1)]);
			if nrm == 0
				continue;
			end
			v_wave_perp = [-v_wave(2) v_wave(1)]/nrm;
			c = fd_x(i,j)*v_wave_perp(1) + fd_y(i,j)*v_wave_perp(2);
			fd_x(i,j) = c*v_wave_perp(1);
			fd_y(i,j) = c*v_wave_perp(2);
		end
	end
end

% Converts a grid index to a position in in the unit square.
function pos = ind2pos(i,j)
	global gridres;
	pos = ([i j] - [.5 .5])/gridres;
end

% A helper function for bilinear interpolation.
% Computes grid index to the bottom left, and bilinear interpolation coeffs.
function [i, j, alpha, beta] = pos2ind(pos)
	global gridres;
	pos = pos*gridres + [.5 .5];
	
	i = floor(pos(1));
	j = floor(pos(2));
	alpha = pos(1) - i;
	beta = pos(2) - j;
end

% Bilinearly interpolate velocity from a given position
function q_itp = interpolate_value(pos, q)
	global gridres;
	% Find interpolation coefficients.
	[i, j, alpha_x, beta_x] = pos2ind(pos);

	wrap = @(x) (1 + mod(x-1, gridres));
	ip1 = wrap(i+1);
	i = wrap(i);
	jp1 = wrap(j+1);
	j = wrap(j);
		
	% Bilinear interpolation.
	q_itp = [1-alpha_x alpha_x]*[q(i,j) q(i, jp1);q(ip1,j) q(ip1,jp1)]*[1-beta_x ; beta_x];
end

% Advection with semi-lagrangian method. Vectorized edition.
function q_new = advect_vectorized(q, v_x, v_y, dt)
	global gridres;
	% Repeat values to achieve periodic boundary conditions.
	[xx, yy] = meshgrid(linspace(.5/gridres-1,2-.5/gridres,gridres*3), linspace(.5/gridres-1,2-.5/gridres,gridres*3));
	q = repmat(q, 3);
	
	[indList_x, indList_y] = meshgrid(1:gridres, 1:gridres);
	indList = [indList_x(:) indList_y(:)];
	X_curr = (indList - [.5 .5])/gridres;
	X_prev = X_curr - dt*[v_x(sub2ind(size(v_x),indList(:,1), indList(:,2))) ...
						  v_y(sub2ind(size(v_y),indList(:,1), indList(:,2)))];		% Explicit euler
	q_new = interp2(xx, yy, rot90(fliplr(q)), X_prev(:,1), X_prev(:,2), 'bilinear');
	q_new = rot90(fliplr(reshape(q_new, gridres, gridres)));
end

% Projects a frequency domain signal into a divergence-free space. Vectorized edition.
function [fd_x, fd_y] = project_out_divergence_vectorized(fd_x, fd_y)
	global gridres;
	[indList_x, indList_y] = meshgrid(1:gridres, 1:gridres);
	indList = [indList_x(:) indList_y(:)];
	v_wave = (indList - [1 1])/gridres - [.5 .5];
	v_wave_perp = [-v_wave(:,2) v_wave(:,1)];
	norms = sqrt(v_wave_perp(:,1).^2 + v_wave_perp(:,2).^2);
	goodInd = norms ~= 0;
	v_wave_perp(goodInd,:) = v_wave_perp(goodInd,:) ./ norms(goodInd);
	c = fd_x(sub2ind(size(fd_x),indList(:,1), indList(:,2))) .* v_wave_perp(:, 1) + ...
	    fd_y(sub2ind(size(fd_y),indList(:,1), indList(:,2))) .* v_wave_perp(:, 2);
	fd_x = rot90(fliplr(reshape(c.*v_wave_perp(:,1), gridres, gridres)));
	fd_y = rot90(fliplr(reshape(c.*v_wave_perp(:,2), gridres, gridres)));
end
