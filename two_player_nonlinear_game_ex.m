%% Nonlinear two-player game searching for Nash equilibrium
% min_x max_y (1/2)<Qx,x> + <Kx,y>
clf; clc

% Some parameters

seed = 37;
m = 500; % Dimension of y

rng(seed,'twister');

k = 1000;
n = 2500;

% Generate problem
A = randn(k, n); % independent, from std normal dist
K = rand(m, n)*2 - 1; % independently, uniformly from [-1,1]
Q = A'*A; % This is n by n

%% Algorithm specifications

N = 20000; % Number of iterations

% Specify algorithm's starting points 
update2_x = rand(n,1);
update2_x = update2_x/sum(update2_x);
update2_y = rand(m,1);
update2_y = update2_y/sum(update2_y);

update2 = [update2_x; update2_y];
anchor = update2;

z_temp = zeros(n+m,1);

% Gamma, R (Lipschitz constant)

gamma = zeros(1,N);
c_matrix = zeros(1,N);
delta = (1/100)*(exp(1) - 1); % Mess around with this at 1/10 and -gamma
% rho = 0, don't need to include it here

G = [Q, K'; -K, zeros(m,m)]; % saddle gradient operator
R = norm(G); % Lipschitz constant of linear operator is its spectral norm
alpha = 1/(2*R);
c_matrix(1) = exp((pi^2)/6); % computed explicitly, and then update2d in algorithm

%% Algorithm
tic;
for j=2:N
    % algorithm steps
    T = update2(:,j-1) - (alpha)*G*update2(:,j-1);
    T(T<0) = 0;
    T(1:n) = T(1:n)/sum(T(1:n));
    T(n+1:end) = T(n+1:end)/sum(T(n+1:end));
    z_temp = update2(:,j-1) + (1/j)*(anchor(:,j-1) - update2(:,j-1)) - (1 - 1/j)*alpha*(update2(:,j-1) - T);
    % 'Normalizing' back to n- and m- simplices
    % z_temp(z_temp<0) = 0; % Return to nonnegative orthant
    % z_temp(1:n) = z_temp(1:n)/sum(z_temp(1:n)); % Return x to n-simplex
    % z_temp(n+1:end) = z_temp(n+1:end)/sum(z_temp(n+1:end)); % Return y to m-simplex
    
    z_temp = z_temp - (alpha)*G*z_temp;
    z_temp(z_temp<0) = 0;
    z_temp(1:n) = z_temp(1:n)/sum(z_temp(1:n));
    z_temp(n+1:end) = z_temp(n+1:end)/sum(z_temp(n+1:end));
    update2(:,j) = update2(:,j-1) + (1/j)*(anchor(:,j-1)) - alpha*(z_temp - G*z_temp);
    % 'Normalizing' back to n- and m- simplices
    % update2(update2<0) = 0; % Return to nonneg orthant, can be done more efficiently?
    % update2(1:n,j) = update2(1:n,j)/sum(update2(1:n,j)); % Return x to n-simplex
    % update2(n+1:end,j) = update2(n+1:end,j)/sum(update2(n+1:end,j)); % Return y to m-simplex

    gamma(j) = -((j+1)*delta)/(c_matrix(j-1)); % Tunable, primary anchor parameter
    % Do an if-else for how we want delta to be update2d
    % if j<400
    % misc terms; may comment out for fixed anchor iteration
        c_matrix(j) = (c_matrix(j-1)/(1 + delta)); % Parameter for gamma
        delta = (1/100)*(exp(1/(j*j)) - 1); % Needed for both gamma and c, which are in turn needed for anchor
    % else
    %     c_matrix(j) = (c_matrix(j-1)/(1 + delta));
    %     delta = (1/100)*(exp(1/(j*j)) - 1);      
    % end
    % Anchor steps
    update2(:,j) = update2(:,j) - (1/(2*R))*G*update2(:,j);
    update2(update2<0) = 0;
    update2(1:n,j) = update2(1:n,j)/sum(update2(1:n,j));
    update2(n+1:end,j) = update2(n+1:end,j)/sum(update2(n+1:end,j));
    anchor(:,j) = anchor(:,j-1) - gamma(j)*(update2(:,j) - G*update2(:,j));
    
    % For plotting; make sure to change each color when you change
    % gamma/update2 name
    % plot(log(j),2*log(norm(G*update2(:,j))),'k*')
    % hold on;
end
toc
%%
for j=1:N
    plot(log(j),2*log(norm(G*update(:,j))),'b',log(j),log(norm(G*update1(:,j))),'k',log(j),log(norm(G*update2(:,j))),'r')
end
title('FEG anchors on 2 player nonlinear game')