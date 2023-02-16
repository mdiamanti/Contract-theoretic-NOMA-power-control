function [p_uni,q_uni,r_uni,cell_utility_uni,user_utility_uni,t] = continuum_contract(N,g,I0,B)

t = g;                                                                  % user's type
t = (t-min(t))*(2-1)/(max(t)-min(t)) + 1;                               % normalized user's type attempts

f_uni = @(t1)unifpdf(t1,1,2);                                           % pdf of uniform distribution over the interval [1,2]
F_uni = @(t2)unifcdf(t2,1,2);                                           % cdf of uniform distribution over the interval [1,2]

% constants
c = 0.1;                                                                % BS's unit cost
h = 1;                                                                  % user's unit cost
a = 1/h;

r_uni = zeros(1,N);                                                     % user's contract theoretic reward
q_uni = zeros(1,N);                                                     % user's contract theoretic effort
p_uni = zeros(1,N);                                                     % user's uplink transmission power

% uniform users' type distribution
%{
for i = 1:N
    fun = @(x) -h*log(1+t(i)*x)*f_uni(t(i)) + c*x*f_uni(t(i)) + h*(1-F_uni(t(i)))*(x/(1+t(i)*x)); % objective function
    x0 = 1;                                                             % optimization's starting point  
    options = optimset('Display','iter');                               % optimization options
    [x,fval,exitflag,output] = fminsearch(fun,x0,options);              % call fminsearch
    r_uni(i) = x;
end
%}

r_uni = (t - 2*a*c + sqrt(t.^2 + 4*a*c.*t - 4*a*c*2)) ./ (2*a*c.*t);
q_uni = 0.2.*r_uni;                                                     
p_uni = 0.5e-3./q_uni;                                                

cell_utility_uni = q_uni - c*r_uni;                                      % cell's utility
user_utility_uni = log(1+t.*r_uni) - (1/h)*q_uni;                        % user's utility

end
