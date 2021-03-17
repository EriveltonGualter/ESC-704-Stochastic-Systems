function [xdot] = tlp_dyn(t,x)
	% two link pendulum dynamics in state space form with x = (q1,q2,q1d,q2d)
	global tlp
	
	g = 9.81;
	m1 = tlp.m1;
	m2 = tlp.m2;
	d1 = tlp.d1;
	d2 = tlp.d2;
	L1 = tlp.L1;
	I1 = tlp.I1;
	I2 = tlp.I2;

	% extract q and qd from x
	q = x(1:2);
	qd = x(3:4);

	% get torques from controller
	u = tlp.controller(t,x);

	% add passive torques when q<qmin or q>qmax
    k = 10000;               % stiffness (Nm/rad)
    b = 100;                 % damping (Nms/rad)
    qmin = [pi/6 ; 0];
    qmax = [5*pi/6 ; 5*pi/6];
    u = u + (q<qmin).*(k*(qmin-q) - b*qd) + ... 		                
                (q>qmax).*(k*(qmax-q) - b*qd);

	% equations of motion
	M = [I1+I2+m1*d1^2+m2*(d2^2+L1^2+2*d2*L1*cos(q(2))) , ...
                I2+d2*m2*(d2+L1*cos(q(2))); ...
                I2+d2*m2*(d2+L1*cos(q(2))), ...
                I2+m2*d2^2];
	C = [-d2*L1*m2*sin(q(2))*qd(2) -d2*L1*m2*sin(q(2))*(qd(1)+qd(2)); ...
          L1*sin(q(2))*qd(1)        0];
	g = [d1*g*m1*cos(q(1))+g*m2*(L1*cos(q(1))+d2*cos(q(1)+q(2))); ...
		 d2*m2*g*cos(q(1)+q(2))];
	qdd = M\(u - C*qd - g);

	% state derivative
	xdot = [qd ; qdd];
end
