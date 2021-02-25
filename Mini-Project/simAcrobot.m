function simAcrobot(X_sys,p)
    clearvars -global

    % Parameters
    l1 = 1;     % Lenght, m
    l2 = 1;     % Lenght, m

    % Simulation parameters
    dt          = p.dt;              % Integration step size [s]
    iteration   = p.iteration;
    t           = 0:dt:dt*iteration; % array time

    % Unpack Angles
    q1 = X_sys(1,:);
    q2 = X_sys(3,:);

    % Forward Kinematics
    x1 = l1*cos(q1);
    y1 = l1*sin(q1);
    x2 = l1*cos(q1) + l2*cos(q1+q2);
    y2 = l1*sin(q1) + l2*sin(q1+q2);

    pos = [x1; y1; x2; y2];

    % Create figure;
    figure; box on;

    time = 0;
    tic;
    while time < t(end)

        % Compute the position of the system at the current real world time
        posDraw = interp1(t',pos',time')';
        px1 = posDraw(1);
        py1 = posDraw(2);
        px2 = posDraw(3);
        py2 = posDraw(4);

        if not(isnan(posDraw))
            draw(px1, py1, px2, py2, [0.2, 0.7, 0.2]);
        end

        % Update current time
        time = toc;
    end


