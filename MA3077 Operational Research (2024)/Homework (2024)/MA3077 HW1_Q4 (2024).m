%% clean up and set parameters
clear; close all;
n = 50; % number of modes
nz = 5;  % number of nonzero weights
m = 23; % number of measurements

%% generate random signal
rng(10,'twister');          % seed random number generator
r = sort(randperm(n,nz)');   % select random modes
w = 2*rand(1,nz)-1;          % select random weights
t_ = linspace(0,pi);
y_ = w*cos(r*t_);           % create signal
plot(t_,y_,'b','LineWidth',3); hold on, axis tight
%title('Original signal','Interpreter','latex','FontSize',18)
%exportgraphics(gcf,'signal.png')


%% extract random measurements
idx = sort(randi([1 length(t_)], 1, m));
t = t_(idx)';
y = y_(idx)';
%figure
plot(t, y, 'r.','MarkerSize',20)
%title('Measurements','Interpreter','latex','FontSize',18)
%exportgraphics(gcf,'samples.png')

%% set up and solve linear programming problems
% control variable is x = [w; t]
% objective function is f(w,t) = t
f = [zeros(1, n) 1];

% inequality constraints: - w - t*ones(n,1) <= 0, w - t*ones(n,1) <= 0
% we rewrite them as [-I -ones(n,1); I -ones(n,1)]*x <= 0, where
% I is the nxn identity matrix
Aie =[-eye(n) -ones(n,1); eye(n) -ones(n,1)];
bie = zeros(2*n, 1);

% equality constraints: Ax + 0t = y
% we rewrite it as [A 0]*x = y
A = cos(t*(1:n));
Aeq = [A zeros(m,1)];
beq = y;
lb = [-inf*ones(n,1);0];
ub = [inf*ones(n,1);inf];

% solve the problem
x = linprog(f, Aie, bie, Aeq, beq,lb,ub);
y_recovered = x(1:n)'*cos((1:n)'*t_);    % recovered signal
plot(t_, y_recovered,'m','linewidth',1.5)

% compare the original and the retrieved weights
fprintf("The number of nonzero weights of the solution is %d \n",sum(abs(x(1:50)) > 1e-14))

%exportgraphics(gcf,'l_inf_min.png')


