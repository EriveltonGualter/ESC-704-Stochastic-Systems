function [px,py,pz] = getSegment(A,B,t,time,i)
    if i==0
        max = interp1(t',A.x',time')';
        may = interp1(t',A.y',time')';
        maz = interp1(t',A.z',time')';

        mbx = interp1(t',B.x',time')';
        mby = interp1(t',B.y',time')';
        mbz = interp1(t',B.z',time')';
    else
        max = A.x(i);
        may = A.y(i);
        maz = A.z(i);

        mbx = B.x(i);
        mby = B.y(i);
        mbz = B.z(i);        
    end
    
    px = [mbx max];
    py = [mby may];
    pz = [mbz maz];
end