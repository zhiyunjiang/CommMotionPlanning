function [x_star, opt] = newton(f, x0, max_iter, tolerance)

    x = x0;
    iter = 0;
    converged = 0;
    step=1e-12;
    dx = step*ones([length(x0), 1]);
    g = @(x) (f(x+dx) - f(x-dx))/(2*step);
    while (iter <= max_iter) && ~converged
       x_next = x - f(x)/g(x);
       if abs(x_next-x) < tolerance
           converged = 1;
           x_start = x;
       end
       x_next = x;
    end
    opt = f(x_start);
end