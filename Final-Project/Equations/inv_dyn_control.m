function u = inv_dyn_control(t,q,qd,kp,kd)

    hp      = 1.70;     % Person Height, m
    A       = 0.2465;   % Forearm length, m
    a       = 0.1060;   % Forearm CM location, m
    B       = 0.3196;   % Upper arm length, m
    b       = 0.1393;   % Upper arm CM location, m
    R1      = 0.2794;	% Handrim radius , m
    R2      = 0.3048;	% Rear wheel radius (R2), m
    Jr      = 0.1395;   % Moment of inertia of rear wheel, kg
    Jb      = 0.020;    % Forearm moment of inertia, kg
    Ja      = 0.021;    % Upper arm moment of inertia, kg
    mp      = 70;       % Person mass, kg
    ma      = 1.96;     % Upper arm mass, kg
    mb      = 1.54;     % Forearm mass (with hand), kg
    mr      = 3.00;     % Rear wheel mass, kg
    mcc     = 12.0;     % Complete wheelchair mass, kg
    g       = 9.81;     % Gravity acceleration, m/s2
    h       = 0.05;     % Shoulder to axle distance - horizontal, m
    v       = 0.73;     % Shoulder to axle distance - vertical, m
    Frol    = 15;       % Total rear wheels rolling resistance, N
    
    x       = q(1);
    bt      = q(2);
    ap      = q(3);
    xdot    = q(4);
    btdot   = q(5);
    apdot   = q(6);

    xd       = qd(1);
    btd      = qd(2);
    apd      = qd(3);
    xdotd    = qd(4);
    btdotd   = qd(5);
    apdotd   = qd(6);
    xddotd    = qd(7);
    btddotd   = qd(8);
    apddotd   = qd(9);
    
    gamma = x/R2;
    gammad = xd/R2;
    ni = (0*pi/180);
    gx = g*sin(ni);
    gy = g*cos(ni);

    
    epos = [x-xd; bt-btd; ap-apd];
    evel = [xdot-xdotd; btdot-btdotd; apdot-apdotd];

%     PID: kp*ePos + kd*eVel + ki*eAcc;
    ad = [xddotd; btddotd; apddotd] - kp*epos - kd*evel;
    
    % Desired
    [M,K,kg,G,H,Q] = dyn_propulsion(A,B,Frol,Ja,Jb,Jr,R1,R2,a,ap,apdot,b,bt,btdot,gamma,gx,gy,ma,mb,mcc,mr);
%     [M,K,kg,G,H,Q] = dyn_propulsion(A,B,Frol,Ja,Jb,Jr,R1,R2,a,apd,apdotd,b,btd,btdotd,gammad,gx,gy,ma,mb,mcc,mr);

    ud = M*ad + K - kg - Q;
    
%     u = horzcat(sum(2*G,2), [0;1;0], [0;0;1])\ud;
%     u(2:4) = u;
    u = ud;
    
%     ad = [qd(4:6).'; M\(kg + Q - K)];
    
%     gamma = x/R2;
%     ni = (0*pi/180);
%     gx = g*sin(ni);
%     gy = g*cos(ni);
    
%     [M2,K2,kg2,G2,H2,Q2] = dyn_propulsion(A,B,Frol,Ja,Jb,Jr,R1,R2,a,ap,apdot,b,bt,btdot,gamma,gx,gy,ma,mb,mcc,mr);

    
end
