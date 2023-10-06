% clf;clc

% Testing simple convex-concave function example
% FEG - Moving anchor5
% L(x,y) = \epsilon(||x||^2)/2 + <x,y> - \epsilon(||y||^2)/2
% x,y real numbers.
% Instances to check:
% postive, negative, and zero rho;
% positive and negative gamma term;
% epsilon = 0.01 and epsilon = 0.

%% FEG with moving anchor5, gamma values change, rho values change (may be 0)

N = 2000; % iterations
update5 = zeros(2,N); % Matrix where x and y values are stored, is plotted.
anchor5 = zeros(2,N); % similar as above but for anchor5 terms

eps = 0.01; % Set equal to zero at some point?

% Initial point:
update5(1,1) = 1;
update5(2,1) = 0;
anchor5(1,1) = update5(1,1);
anchor5(2,1) = update5(2,1);

% gamma, R, and rho
gamma = zeros(1,N);
c_matrix = zeros(1,N);
delta = (exp(1) - 1);
R = sqrt(2 + 2*eps);
rho = 0;
alpha = (1/R); % Luckily is constant throughout this time around
c_matrix(1) = 10; % Check notes, need to see what this can be

for j=2:N
    % update5 primary terms
    update5(1,j) = update5(1,j-1) + (1/j)*(anchor5(1,j-1) - update5(1,j-1)) - (1/R)*(eps*(update5(1,j-1) + (1/j)*(anchor5(1,j-1) - update5(1,j-1)) - (1 - (1/j))*((1/R) + 2*rho)*(eps*update5(1,j-1) + update5(2,j-1))) + (update5(2,j-1) + (1/j)*(anchor5(2,j-1) - update5(2,j-1)) - (1 - (1/j))*((1/R) + 2*rho)*(eps*update5(2,j-1) - update5(1,j-1)))) - 2*rho*(1 - (1/j))*(eps*update5(1,j-1) + update5(2,j-1));
    update5(2,j) = update5(2,j-1) + (1/j)*(anchor5(2,j-1) - update5(2,j-1)) - (1/R)*(eps*(update5(2,j-1) + (1/j)*(anchor5(2,j-1) - update5(2,j-1)) - (1 - (1/j))*((1/R) + 2*rho)*(eps*update5(2,j-1) - update5(1,j-1))) - (update5(1,j-1) + (1/j)*(anchor5(1,j-1) - update5(1,j-1)) - (1 - (1/j))*((1/R) + 2*rho)*(eps*update5(1,j-1) + update5(2,j-1)))) - 2*rho*(1 - (1/j))*(eps*update5(2,j-1) - update5(1,j-1));

    % update5 some misc terms and constants we need
    % Terms for anchor5 update5
    gamma(j) = 0;%((j+1)*delta)/(c_matrix(j-1)); % primary parameter for updating anchor5
    c_matrix(j) = (c_matrix(j-1)/(1 + delta)); % Parameter for gamma
    delta = exp(1/(j*j)) - 1; % Needed for both gamma and c, which are in turn needed for anchor5 update5

    % update5 the anchor5
    anchor5(1,j) = anchor5(1,j-1) + gamma(j)*(eps*update5(1,j) + update5(2,j));
    anchor5(2,j) = anchor5(2,j-1) + gamma(j)*(eps*update5(2,j) - update5(1,j));
end

% % This code is for plotting the algorithm's error
% plot(log(1:N),2*log(sqrt(update5(1,:).^2 + update5(2,:).^2)),'r');
% title('Error, zero \gamma, \epsilon=0.01, \rho<0') % update5 with each running

% % This code is for plotting the algorithm iterations
% plot(update5(1,:),update5(2,:),'LineWidth',2);
% hold on;
% for j=1:N
%     plot(update5(1,j),update5(2,j),'r*');
% end
% plot(anchor5(1,:),anchor5(2,:),'c');
% for j=1:N
%     plot(anchor5(1,j),anchor5(2,j),'g*');
% end
% plot(0,0,'k*');
% title('FEG moving anchor5, positive \gamma, \epsilon=0.01, \rho=0')
% set(gca,'linewidth',2)