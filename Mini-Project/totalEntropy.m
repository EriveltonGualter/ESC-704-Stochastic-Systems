function [entropy] = totalEntropy(Sk , del_uk, p)
    % Entropy code from : Williams, G., Drews, P., Goldfain, B., Rehg, J.
    % M., & Theodorou, E. A. (2016, May). Aggressive driving with model
    % predictive path integral control. In 2016 IEEE International
    % Conference on Robotics and Automation (ICRA) (pp. 1433-1440). IEEE.
        
    n = length(Sk);
    lambda = p.lambda;
    sum1 = 0;
    sum2 = 0;
    for i = 1:n
        sum1 = sum1 + exp(-(1/lambda)*Sk(i))*del_uk(i);
        sum2 = sum2 + exp(-(1/lambda)*Sk(i));
    end
    entropy = sum1/sum2;

end