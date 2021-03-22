clear all
clc

syms t real

syms x bt ap xc a b A B R1 R2 xdot btdot apdot

% gamma = x/R2;

rc = [x+xc];
rb = [x+b*sin(bt); -b*cos(bt)];
ra = [x+B*sin(bt)+a*sin(bt+ap); -B*cos(bt)-a*cos(bt+ap)];
rr = [x+xc];
thb = bt;
tha = bt+ap;
thr = -x/R2;

F = [rc; rb; ra; rr; thb; tha; thr];
J = jacobian(F,[x,bt,ap]);
Jt = subs(J,{x,bt,ap},{str2sym('x(t)'),str2sym('bt(t)'),str2sym('ap(t)')});
Jdot = diff(Jt,t);

x = sym(str2sym('x(t)'));
bt = sym(str2sym('bt(t)'));
ap = sym(str2sym('ap(t)'));
Jdot = subs(Jdot,{diff(x,t),diff(bt,t),diff(ap,t)},{xdot,btdot,apdot});
Jdot = subs(Jdot,{x,bt,ap},{sym('x'),sym('bt'),sym('ap')});

q = [x; bt; ap];
qdot = diff(q,t);
qddot = diff(qdot,t);

rc = [x+xc];
rb = [x+b*sin(bt); -b*cos(bt)];
ra = [x+B*sin(bt)+a*sin(bt+ap); -B*cos(bt)-a*cos(bt+ap)];
rr = [x+xc];
thb = bt;
tha = bt+ap;
thr = -x/R2;

rcdot = diff(rc,t);
rbdot = diff(rb,t);
radot = diff(ra,t);
rrdot = diff(rr,t);
thbdot = diff(thb,t);
thadot = diff(tha,t);
thrdot = diff(thr,t);

rcddot = diff(rcdot,t);
rbddot = diff(rbdot,t);
raddot = diff(radot,t);
rrddot = diff(rrdot,t);
thbddot = diff(thbdot,t);
thaddot = diff(thadot,t);
thrddot = diff(thrdot,t);

syms mcc mb ma mr Jb Ja Jr gx gy gamma Fx Fy tauo tauc Frol
syms x bt ap 

Mb = diag([mcc mb mb ma ma mr Jb Ja Jr]);

fe1 = [-mcc*gx; -mb*gx; -mb*gy; -ma*gx; -ma*gy; -mr*gx; 0; 0; 0];
fe2 = [zeros(3,2); -eye(2); diag([1 0]); -(A-a)*cos(ap) -(A-a)*sin(ap); -R1*sin(gamma) -R1*cos(gamma)];
fe3 = [zeros(6,2); 1 -1; 0 1; 0 0];
fe4 = [1; zeros(8,1)];

M = simplify(J.'*Mb*J);
K = simplify(J.'*Mb*Jdot*[xdot; btdot; apdot]);
kg = simplify(J.'*fe1);
G = simplify(J.'*fe2);%*[Fx; Fy];
H = simplify(J.'*fe3);%*[tauo; tauc];
Q = simplify(J.'*fe4)*Frol;

Ke = simplify( kg + G + H + Q);

syms h v
c = [b*sin(bt)+a*sin(ap)-h+R1*cos(x/R2); b*cos(bt)+a*cos(ap)-v+R1*sin(x/R2)]; 


dyn = matlabFunction(M,K,kg,G,H,Q,'File','dyn_propulsion');

help dyn_propulsion

