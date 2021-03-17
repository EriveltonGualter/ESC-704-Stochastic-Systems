function tlp_passive
	% simulates a 2-link pendulum
	
    global tlp
    
 	% set pendulum parameters for stance simulation, somewhat human-like
	bodymass = 80.0;						% total body mass
	tlp.L1 = 0.5;							% length of shank segment
	tlp.m1 = 2*0.0465 * bodymass;  			% mass of left shank + right shank using Winter's anthropometry Table
	tlp.d1 = 0.567*tlp.L1;					% center of mass
	tlp.I1 = tlp.m1*(0.302*tlp.L1)^2;		% moment of inertia from Winter's anthropometry Table
	tlp.L2 = 0.5;							% length of thigh segment
	tlp.m2 = bodymass-tlp.m1;               % rest of body is all in the thigh
	tlp.d2 = tlp.L2;    					% center of mass of thigh is at the hip (due to rest of body)
	tlp.I2 = 0.2*bodymass*(0.323*tlp.L2)^2;	% moment of inertia from Winter's anthropometry Table
	
	% simulation settings
    duration = 6.0;
    interval = 0.05;            % time interval for sampling of simulation
	x = [pi/3 pi/3 0 0]';		% initial state (q1,q2,q1d,q2d)
	tlp.controller = @controller_passive;

	% run a simulation of the system
	t = 0:interval:duration;				% times where simulation must be sampled
	nsamples = numel(t);
	tlp_draw(x,0);
	disp('Hit enter to run simulation');
	pause
    for i = 2:nsamples
        [t1,x1] = ode23(@tlp_dyn,[t(i-1) t(i)],x);    
        x = x1(end,:)';                           % only use the final state from this time interval
		tlp_draw(x,t(i));
	end
	
	% determine the stability of the system when in upright position
	disp('Hit ENTER to start stability analysis...');
	pause
    df_dx = zeros(4,4);         % initialize Jacobian
    x = [pi/2; 0; 0; 0];        % the state we are interested in
    f = tlp_dyn(0, x);              % the function f at this state
    hh = 1e-7;                  % a small finite difference
    for i=1:4
       xsave = x;
       x(i) = x(i) + hh;                	% add hh to state variable i
       df_dx(:,i) = (tlp_dyn(0,x) - f)/hh;  % column i of the Jacobian
       x = xsave;
    end
    [eigenvectors,eigenvalues] = eig(df_dx)         
    lambda = diag(eigenvalues) + 1e-10*j;	% add a small imaginary part so "plot" will use complex plane
    figure(1)
    semilogx(lambda,'o');
    title(['eigenvalues at x = [' num2str(x') ']']);
    xlabel('REAL (s^-1)');
    ylabel('IMAGINARY (s^-1)')

end
%=======================================================================
function u = controller_passive(t,x)
	% controller, generates torques u as a function of time t and state x
	u = [0;0];		% no torques
end
%=======================================================================
function u = controller_pd(t,x)
	% controller, generates torques u as a function of time t and state x
    q = x(1:2);			% angles
    qd = x(3:4);		% angular velocities
    qr = [pi/2 ; 0];	% desired posture
    Kp = 1000;          % proportional gain
    Kd = 100;           % derivative gain
	u = -Kp*(q-qr) - Kd*qd;
end
