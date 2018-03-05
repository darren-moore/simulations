pico-8 cartridge // http://www.pico-8.com
version 16
__lua__

-- the 1d heat equation is:
-- u_t = k*u_xx

-- this cart solves it on the
-- real line given initial conditions.
-- we use the fundamental solution
-- aka the heat kernel.

function _init()
	poke(0x5f2d, 1)
	simulating = false
	
	u = {}
	init = {}
	
	for x=-64,64 do
		u[x] = 0
		init[x] = 0
	end
	
	past_mouse_x = 0
	past_mouse_y = 0
	kernel = {}
	k = 1
	t = .001
	
	-- no boundary conditions
	-- so lets scale to make nicer
	scaler = kernel_integral(.001)
end

function _update60()
	
	if(btn(4)) then
		-- reset
		simulating = false
		t = .001
		for x=-64,64 do
			init[x] = 0
		end
	end
	
	if(simulating) then
		sim_update()
	else
		no_sim_update()
	end
		
	past_mouse_x = stat(32)
	past_mouse_y = stat(33)
		
end

function sim_update()
	-- we advance only a few
	-- points per frame for speed
	-- it also looks cooler
	
	-- precompute kernel per time
	kernel = compute_kernel(t)
	for i=1,10 do
		x = flr(rnd(128)) - 64
		u[x] = convolve(x,kernel,init)*1/scaler
	end
	
	t+=.002

end

function no_sim_update()
 if(stat(34) == 1) then
		start_x = past_mouse_x
		end_x = stat(32)
		if(past_mouse_x > stat(32)) then
			start_x = stat(32)
			end_x = past_mouse_x
		end
		
		for i=start_x,end_x do
			y = past_mouse_y
			init[i - 64] = y/64 - 1
		end
		
	end
			
	if(btn(5)) then
		simulating = true
		kernel = compute_kernel(t)
		-- compute all points first
		for x=-64,64 do
		 u[x] = convolve(x,kernel,init)*1/scaler
	 end
	end
	
end

function _draw()
	cls(7)
	if(simulating) then
		for x=-64,64 do
			rect(x+64, 64, x+64, 64+u[x]*64, cmap(u[x]))
		end
	else
		for x=-64,64 do
			rect(x+64, 64, x+64, 64+init[x]*64, cmap(init[x]))
		end
	end
	
	rect(0,64,128,64, 6)
	circfill(stat(32), stat(33), 1,0)
	print("draw some initial conditions!", 6, 1)
	print("❎ to simulate", 2,122)
	print("🅾️ to reset", 83, 122)
end



-->8

function kernel_integral(t)
-- a little hack to make
-- the solution smoother
	total = 0
	for y=-64,64 do
		val = heat_kernel((0 - y)/64,t)
		total += 1/128 * val
	end

	return(total)

end

function convolve(x,f,g)
	total = 0
	for y=-64,64 do
		val = f[x-y]*g[y]
		total += 1/128 * val
	end

	return(total)

end


function compute_kernel(t)

	kernel = {}
	for y=-128,128 do
		kernel[y] = heat_kernel(y/64,t)
	end
	
	return(kernel)
end

function heat_kernel(x,t)
	c = 1/sqrt(4*3.14*t*k)
	power = -(x^2)/(4*k*t)
	-- lua can't do negative powers?
	return(c*(1/(2.71^-power)))
	
end
-->8
function cmap(x)
	x = abs(x)
	if(x > .6 ) then
		return 8
	elseif(x > .5) then
		return 9
	elseif(x > .4) then
		return 10
	elseif(x > .3) then
		return 11
	elseif(x > .2) then
		return 12
	elseif(x > .1) then
		return 13
	elseif(x > .05) then
		return 1
	else
		return 0
	end
end
