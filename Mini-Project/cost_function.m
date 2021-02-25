function [S] = cost_function(z, u, p)

    % Unpack
    R   = p.R;    
    dt  = p.dt;
    z   = z - p.xT.';
    
    % Quadractic Cost
    S = (z.'*R*z)*dt;    
end

