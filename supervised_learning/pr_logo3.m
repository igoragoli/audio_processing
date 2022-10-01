function [x,y] = pr_loqo3(c, H, A, b, l, u, use_octave)
%[X,Y] = PR_LOQO3(c, H, A, b, l, u)
%
% loqo solves the quadratic programming problem
%
% minimize   c' * x + 1/2 x' * H * x
% subject to A*x = b
%           l <= x <= u
%
% for a documentation see R. Vanderbei, LOQO: an Interior Point Code for Quadratic Programming
%
%%%%%%%%%%%%%%%%
% modifié par Stéphane Rossignol pour "octave"
% - ajout de use_octave à cause de la fonction "error" ci-dessous
% - ajout des "fprintf(1,);

%%%fprintf(1,"pr_loqo3 - on commence\n");

if use_octave==0
  error('')
end;

% the fudge factors
%%%fprintf(1,"pr_loqo3 - the fudge factor\n");
margin = 0.05;
%bound  = 100;
bound  = 0.5*max(u-l); %for the starting point
sigfig_max = 8;   % It has been verified that at deg=10^2, we need this number
                  % of sig figures for the 2 class gaussian cluster to correctly
                  % determine which are margin and non-margin SVs
counter_max = 50;

% gather some starting data
%%%fprintf(1,"pr_loqo3 - gather some starting data\n");
[m, n] = size(A); 
% this is done in order not to mess up H but only to tamper with H_x
H_x    = H;
H_diag = diag(H);

b_plus_1 = 1;
c_plus_1 = norm(c) + 1;
one_x = -ones(n,1);
one_y = -ones(m,1);

% starting point
for i = 1:n 
  H_x(i,i) = H_diag(i) + 1;
end;
H_y = eye(m);
c_x = c;
c_y = 0;

% and solve the system [-H_x A'; A H_y] [x, y] = [c_x; c_y]
%%%fprintf(1,"pr_loqo3 - chol (1/3)\n");
R = chol(H_x);
%size(R)
%size(A)
%size(c_x)
H_Ac = R \ ([A; c_x'] / R)';
H_A = H_Ac(:,1:m);
H_c = H_Ac(:,(m+1):(m+1));
A_H_A = A * H_A;
A_H_c = A * H_c;
H_y_tmp = (A_H_A + H_y);
y = H_y_tmp \ (c_y + A_H_c);
x = H_A * y - H_c;
%DEBUG : comment out line below
%x %xxx

g = max(abs(x - l), bound);
z = max(abs(x), bound);
t = max(abs(u - x), bound);
s = max(abs(x), bound);

mu = (z' * g + s' * t)/(2 * n);

% set some default values
sigfig = 0;
counter = 0;
alfa = 1;

while ((sigfig < sigfig_max) * (counter < counter_max)),

% update the iteration counter
counter = counter + 1;

%central path (predictor)
H_dot_x = H * x;

rho = - A * x + b;%%%
%DEBUG : Next three lines commented out
%A*x
%b
%rho %xxx
nu = l - x + g;
tau = u - x - t;
sigma = c - A' * y - z + s + H_dot_x;

gamma_z = - z;
gamma_s = - s;

% instrumentation
x_dot_H_dot_x = x' * H_dot_x;

primal_infeasibility = norm([tau; nu]) / b_plus_1;
dual_infeasibility = norm([sigma]) / c_plus_1;

primal_obj = c' * x + 0.5 * x_dot_H_dot_x;
dual_obj = - 0.5 * x_dot_H_dot_x + l' * z - u' * s + b'*y; %%%

old_sigfig = sigfig;
sigfig = max(-log10(abs(primal_obj - dual_obj)/(abs(primal_obj) + 1)), 0);

report = sprintf('counter %i p_i %e d_ii %e sigfig %f, alpha %f, p_o %f d_o %f mu %e', counter, primal_infeasibility, dual_infeasibility, sigfig, alfa, primal_obj, dual_obj, mu);
%disp(report);

% some more intermediate variables (the hat section)
hat_nu = nu + g .* gamma_z ./ z;
hat_tau = tau - t .* gamma_s ./ s;

% the diagonal terms
d = z ./ g + s ./ t;

% initialization before the big cholesky
for i = 1:n H_x(i,i) = H_diag(i) + d(i); end;
H_y = 0;
c_x = sigma - z .* hat_nu ./ g - s .* hat_tau ./ t;
c_y = rho;
%c_y %debug : line removed

% and solve the system [-H_x A'; A H_y] [delta_x, delta_y] = [c_x; c_y]
%%%fprintf(1,"pr_loqo3 - chol (2/3) %d\n",counter);
R = chol(H_x);
H_Ac = R \ ([A; c_x'] / R)';
H_A = H_Ac(:,1:m);
H_c = H_Ac(:,(m+1):(m+1));
A_H_A = A * H_A;
A_H_c = A * H_c;
H_y_tmp = (A_H_A + H_y);
%DEBUG : 5 lines commented out
%size(H_y_tmp)
%size(c_y)
%size(A_H_c)
delta_y = H_y_tmp \ (c_y + A_H_c);
%size(H_A)
%size(delta_y)
delta_x = H_A * delta_y - H_c;

%backsubstitution
delta_s = s .* (delta_x - hat_tau) ./ t;
delta_z = z .* (hat_nu - delta_x) ./ g;

delta_g = g .* (gamma_z - delta_z) ./ z;
delta_t = t .* (gamma_s - delta_s) ./ s;

%central path (corrector)
gamma_z = mu ./ g - z - delta_z .* delta_g ./ g;
gamma_s = mu ./ t - s - delta_s .* delta_t ./ t;

% some more intermediate variables (the hat section)
hat_nu = nu + g .* gamma_z ./ z;
hat_tau = tau - t .* gamma_s ./ s;

% the diagonal terms
%d = z ./ g + s ./ t;

% initialization before the big cholesky
%for i = 1:n H_x(i,i) = H_diag(i) + d(i);
%H_y = 0;
c_x = sigma - z .* hat_nu ./ g - s .* hat_tau ./ t;
c_y = rho;

% and solve the system [-H_x A'; A H_y] [delta_x, delta_y] = [c_x; c_y]
% R = chol(H_x);
%%%fprintf(1,"pr_loqo3 - chol (3/3) %d\n",counter);
H_Ac = R \ ([A; c_x'] / R)';
H_A = H_Ac(:,1:m);
H_c = H_Ac(:,(m+1):(m+1));
A_H_A = A * H_A;
A_H_c = A * H_c;
H_y_tmp = (A_H_A + H_y);
delta_y = H_y_tmp \ (c_y + A_H_c);
delta_x = H_A * delta_y - H_c;

%backsubstitution

delta_s = s .* (delta_x - hat_tau) ./ t;
delta_z = z .* (hat_nu - delta_x) ./ g;

delta_g = g .* (gamma_z - delta_z) ./ z;
delta_t = t .* (gamma_s - delta_s) ./ s;

%compute the updates
alfa = - (1-margin) / min([delta_g ./ g; delta_t ./ t;
                    delta_z ./ z; delta_s ./ s; -1]);

mu = (z' * g + s' * t)/(2 * n);
mu = mu * ((alfa - 1) / (alfa + 10))^2;

x = x + delta_x * alfa;
g = g + delta_g * alfa;
t = t + delta_t * alfa;
y = y + delta_y * alfa;
z = z + delta_z * alfa;
s = s + delta_s * alfa;

end

%disp(report);

%%% les commentaires qui suivent proviennent du code original et pas
%%% du Lagis ou de CentraleSupélec
%Warning: Empty matrix multiplication with unequal inner dimensions.
%> In /home/neuro/bs/matlab/1pr/pr_loqo3.m at line 114
%  In /home/neuro/bs/matlab/1pr/usps.m at line 56
%??? Error using ==> -
%Matrix dimensions must agree.
%
%Error in ==> /home/neuro/bs/matlab/1pr/pr_loqo3.m
%On line 114  ==> delta_x = H_A * delta_y - H_c;
%
%Error in ==> /home/neuro/bs/matlab/1pr/usps.m
%On line 56  ==>   alpha = pr_loqo3(c, H, A, b, bl, bu);

