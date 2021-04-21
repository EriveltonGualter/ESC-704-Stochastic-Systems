function dx = ode_test(t,x,u,p)
    dx = [p*x(2); u];
end