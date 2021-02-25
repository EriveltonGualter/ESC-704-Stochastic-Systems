function dz = sode_function(z, u, param)

% Return acceleration component from acrobot dynamics
acc = acrobotDynamics(z,u).';

dz = [z(2); acc(1); z(4); acc(2)];
end

