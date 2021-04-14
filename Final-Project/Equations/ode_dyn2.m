function qdot = ode_dyn2(t,q,u,p)

    hp      = 1.70;     % Person Height, m
    R1      = 0.2794;	% Handrim radius , m
    R2      = 0.3048;	% Rear wheel radius (R2), m
    mp      = 70;       % Person mass, kg
    g       = 9.81;     % Gravity acceleration, m/s2
    Frol    = 15;       % Total rear wheels rolling resistance, N
    
    A       = p(1); % Forearm length, m
    a       = p(2); % Forearm CM location, m
    B       = p(3); % Upper arm length, m
    b       = p(4); % Upper arm CM location, m
    Jr      = p(5); % Moment of inertia of rear wheel, kg
    Jb      = p(6); % Forearm moment of inertia, kg
    Ja      = p(7); % Upper arm moment of inertia, kg
    ma      = p(8); % Upper arm mass, kg
    mb      = p(9); % Forearm mass (with hand), kg
    mr      = p(10);% Rear wheel mass, kg
    mcc     = p(11);% Complete wheelchair mass, kg
    h       = p(12);% Shoulder to axle distance - horizontal, m
    v       = p(13);% Shoulder to axle distance - vertical, m

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
    
    qdot = [q(4); q(5); q(6); M\(kg + G*[Fx; Fy] + H*[tauo; tauc] + Q - K)];
end
