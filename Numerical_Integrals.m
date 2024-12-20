Problem1()
Problem2()
Problem3()
Problem4()
Problem1A()
function Problem1()
    % Main script part
    a = 0;  
    b = pi;  
    h = 0.01;  
    f = @(x) cos(x);  
    first_derivative = first_derv(a, b, f, h);
    second_derivative = second_derv(a, b, f, h);
    % Define the x values for plotting
    x_values = a:h:b;
    % Plot the functions
    figure;
    hold on;  
    plot(x_values, f(x_values), 'b', 'LineWidth',2);
    plot(x_values, first_derivative, 'r');
    plot(x_values, second_derivative, 'g');
    % Add title and labels
    title('Function and its Derivatives');
    xlabel('x');
    ylabel('y');
    legend('show');  
    hold off;  
end
function Problem2()
    f = @(x) 1 ./ (x.^2 + 0.1);
    a = 0;  
    b = 2;  
    n = 10;  %n has to be even  
    T = trap_int(a, b, f, n);
    fprintf('Trapezoidal Rule Approximation: %.4f\n', T);
    S = simp_int(a, b, f, n);
    fprintf('Simpson''s Rule Approximation: %.4f\n', S);
end
function Problem3()
    f = @(x) cos(x);
    a = 0;  
    b = 2;  
    n = 10;  %n has to be even  
    T = trap_int(a, b, f, n);
    fprintf('Trapezoidal Rule Approximation: %.4f\n', T);
    S = simp_int(a, b, f, n);
    fprintf('Simpson''s Rule Approximation: %.4f\n', S);
end
function Problem4()
    f = @(x) 1 ./(sqrt(x.^3+1));
    a = 0;  
    b = 2;  
    n = 10;  %n has to be even  
    T = trap_int(a, b, f, n);
    fprintf('Trapezoidal Rule Approximation: %.4f\n', T);
    S = simp_int(a, b, f, n);
    fprintf('Simpson''s Rule Approximation: %.4f\n', S);
end
% this is the problem for the integration equation
function Problem1A()
    a = 0;  
    b = 1;  
    h = 0.01;
    f = @(x,t) -x+sin(t);
    first_derivative = first_derv(a, b, f, h);
    second_derivative = second_derv(a, b, f, h);
    % Define the x values for plotting
    x_values = a:h:b;
    % Plot the functions
    figure;
    hold on;  
    plot(x_values, f(x_values), 'b', 'LineWidth',2);
    plot(x_values, first_derivative, 'r');
    plot(x_values, second_derivative, 'g');
    % Add title and labels
    title('Function Derivatives');
    xlabel('x');
    ylabel('t');
    legend('show');  
    hold off;   
end
% this is the function for the 
function d = first_derv(a,b,f,h)
    arguments
        a double = 0;
        b double = 1;
        f function_handle = @(x) sin(x);
        h = 0.01;
    end
    x = a:h:b;
    d = (f(x+h)-f(x-h))/(2*h);
end
function d = second_derv(a,b,f,h)
    arguments
        a double = 0;
        b double = 1;
        f function_handle= @(x) sin(x);
        h = 0.01;
    end
    x = a:h:b;
    d = (f(x+2*h)-2*f(x)+f(x-2*h))/(4*h^2);
end
% this the function that does the trap integration  
function approx = trap_int(a,b,f,n)
    h = (b - a) / n;  
    x = linspace(a, b, n+1); % a:n+1:b 
    y = f(x);  %defined in line1
    approx = (h/2) * (y(1) + 2*sum(y(2:end-1)) + y(end)); 
end
% this is the function that does simpson integration
function approx = simp_int(a,b,f,n)
    h = (b - a) / n; 
    x = linspace(a, b, n+1);  
    y = f(x);  
    approx = (h/3) * (y(1) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2)) + y(end));
end
%{
function eu = euler_forwardmeth(a,b,f,n)
    x = a:n:b;
    t = a:n:b;
    while i<b
       k1 = n*f(t(i),x(t(i)));
       x(i+1) = x(i)+k1;
    end
    eu=x(i+1)
end
function eu = Runge_kutta_meth2nd(a,b,f,n)
    x = a:n:b;
    t = a:n:b;
    while i<b
        k1 = n*f(t(i),x(t(i)));
        k2 = n*f(t(i)+(h/2),x(t(i))+(k1/2));
        k2 = n*f(t(i),x(t(i)));
        k3 = n*f(t(i)+(h/2),x(t(i))+(k1/2));
        x(t(i+1)) = x(t(i))+k2;
    end
    eu=x(t(i+1))
end
function eu = Runge_kutta_meth4th(a,b,f,n)
    x = a:n:b;
    t = a:n:b;
    while i<b
        k1 = n*f(t(i),x(t(i)));
        k2 = n*f(t(i)+(h/2),x(t(i))+(k1/2));
        k2 = n*f(t(i),x(t(i)));
        k3 = n*f(t(i)+(h/2),x(t(i))+(k1/2));
        x(t(i+1)) = x(t(i))+k2;
    end
    eu=x(i+1)
end         
%}