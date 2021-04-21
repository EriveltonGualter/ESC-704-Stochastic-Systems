function out = getAngle2(P1, P2, P3, P4, i)
    u(1) = P1.x(i) - P2.x(i);
    u(2) = P1.y(i) - P2.y(i);
    u(3) = P1.z(i) - P2.z(i);
    
    v(1) = P4.x(i) - P3.x(i);
    v(2) = P4.y(i) - P3.y(i); 
    v(3) = P4.z(i) - P3.z(i);
    
    out = atan2(u,v);
end