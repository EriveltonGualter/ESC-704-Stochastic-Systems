function jumpopt
% optimizes a jumping movement
	global sim
	
	% Some constants
	ndof = 9;
	nmus = 16;
	nstates = 2*ndof + 2*nmus;
	
	% Set initial condition, a sort of squat
	% (it would be better to start in a static equilibrium)
	sim.xinit = [0 0.68 -pi/4+0.1 pi/2 -pi/2 pi/4-0.1 pi/2 -pi/2 pi/4-0.1 ...
	            zeros(1,9)   ...
				1.4*ones(1,16) ...
				zeros(1,16)];
    sim.Neval = 0;

	% do the optimization
	p0 = 0.4 + zeros(8,1);			% initial guess for controller parameters
	tic
	p = fminsearch(@simulate, p0);
	disp('optimal parameters:');
	p
    disp('Hit ENTER to show the optimal solution');
    pause
    simulate(p);
end
%============================================================================
function [minusheight] = simulate(p)
% Simulates a jump
	global sim
	
	% Run the simulation
	tend = 1.0;
	sim.par = p;			% store controller parameters in global variable so controller can access it
	[t,x] = ode23(@odefun, [0 tend], sim.xinit);
    sim.Neval = sim.Neval + 1;
	if (toc > 3)            % do this every 3 seconds
		anim(interp1(t,x,(0:0.005:tend)')');
		hold on; plot(x(:,1),x(:,2),'k');drawnow;
		tic;
	end
	
	height = max(x(:,2));
	fprintf('%d: jump height: %6.4f m (p=%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f)\n', ...
        sim.Neval, height,p);
	minusheight = -height;
end
%=====================================================================================
function [xdot] = odefun(t,x);
	u = controller(t,x);
	xdot = gait2de(x,u);
end
%=====================================================================================
function [u] = controller(t,x);
	global sim
		% the 16 muscles are:
		% R.Iliopsoas			hip flexion	
		% R.Glutei				hip extension			
		% R.Hamstrings			hip extension, knee flexion
		% R.Rectus				hip flexion, knee extension	
		% R.Vasti				knee extension
		% R.Gastroc				knee flexion, ankle plantarflexion
		% R.Soleus				ankle plantarflexion
		% R.TibialisAnt			ankle dorsiflexion
		% L.Iliopsoas		
		% L.Glutei			
		% L.Hamstrings		
		% L.Rectus			
		% L.Vasti			
		% L.Gastroc			
		% L.Soleus			
		% L.TibialisAnt
	
	u = zeros(16,1);
		% parameter 1 determines when the Iliopsoas turns on, etc.
	for i=1:8
		if (t > sim.par(i))
			u(i) = 1;				% turn the muscle on
			u(i+8) = 1;				% and the one on the other side
		end
	end
end
%=====================================================================================
function anim(x);
% animate a simulation result x(t)

	R = [1:6 4];			% right stick points
	L = [2 7:10 8];			% left stick points
	nframes = size(x,2);
	u = zeros(16,1);
	for i=0:5:nframes-1
		clf
		plot([-1 1],[0 0],'k');
		hold on
		[~, ~, d] = gait2de(x(:,i+1),u);    % get the stick figure coordinates for state x
		d = reshape(d',2,10)';
		plot(d(R,1),d(R,2),'r',d(L,1),d(L,2),'b','LineWidth',2);
		axis('equal');
		axis([-1 1 -0.4 2]);
		axis('off');
		drawnow
	end

end