function Tsave = compute_kinetic(I1,I2,Q2,Q1dot,Q2dot,l1,m1,m2,r1,r2)
%COMPUTE_KINETIC
%    TSAVE = COMPUTE_KINETIC(I1,I2,Q2,Q1DOT,Q2DOT,L1,M1,M2,R1,R2)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    24-Feb-2021 22:57:00

t2 = Q1dot+Q2dot;
t3 = Q1dot.^2;
t4 = r2.^2;
t5 = cos(Q2);
Tsave = I1.*t3.*(1.0./2.0)+I2.*t2.^2.*(1.0./2.0)+m2.*(t3.*t4+Q2dot.^2.*t4+l1.^2.*t3+Q1dot.*Q2dot.*t4.*2.0+l1.*r2.*t3.*t5.*2.0+Q1dot.*Q2dot.*l1.*r2.*t5.*2.0).*(1.0./2.0)+m1.*r1.^2.*t3.*(1.0./2.0);
