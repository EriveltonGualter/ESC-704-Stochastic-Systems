function qdot = ode_dyn(t,q,u)

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

    Fx      = u(1);
    Fy      = u(2);
    tauo    = u(3);
    tauc    = u(4);
    
    gamma   = 0;
    gx      = 0;
    gy      = 0;

	[M,K,kg,G,H,Q] = dyn_propulsion(A,B,Frol,Ja,Jb,Jr,R1,R2,a,ap,apdot,b,bt,btdot,gamma,gx,gy,ma,mb,mcc,mr);
    
    clear qdot
    qdot(1,1) = q(4);
    qdot(2,1) = q(5);
    qdot(3,1) = q(6);
    qdot(4:6,1) = M\(kg + G*[Fx; Fy] + H*[tauo; tauc] + Q - K);
 
end