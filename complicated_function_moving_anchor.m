% clf;clc

% Testing more involved nonconvex-nonconcave function example
% FEG - Moving anchor
% L(x,y) = (-1/6)x^2 + (2sqrt(2)/3)xy + (1/6)y^2
% x,y real numbers.
% Instances to check:
% positive, negative, and zero gamma term.

%% FEG complicated fxn all anchoring options

N = 2000; % Iterations
update3 = zeros(2,N); % Matrix of x and y update3s; is plotted
anchor = zeros(2,N); % Matrix of anchor points
z_temp = zeros(2,1); % Matrix of half-steps

% Initial point

update3(1,1) = 1;
update3(2,1) = 0;
anchor(1,1) = update3(1,1);
anchor(2,1) = update3(2,1);

% Gamma, R, and \rho

gamma = zeros(1,N);
c_matrix = zeros(1,N);
delta = (exp(1) - 1);
R = 1; % We call this L in the onenote notes, is Lipschitz constant
rho = -1/3;
alpha = 1;
c_matrix(1) = exp((pi^2)/6); % computed explicitly

% Saddle Gradient operator

A = [(-1/3) (2*sqrt(2))/3; -(2*sqrt(2))/3 (-1/3)];

for j=2:N
    % update3 z_temp, the intermediate update3
    z_temp = update3(:,j-1) + (1/j)*(anchor(:,j-1) - update3(:,j-1)) - ((j-1)/(j))*(1/3)*A*update3(:,j-1);
    update3(:,j) = update3(:,j-1) + (1/j)*(anchor(:,j-1)) - A*(z_temp - ((j-1)/j)*(2/3)*update3(:,j-1));
    
    % update3 misc
    gamma(j) = ((j+1)*delta)/(c_matrix(j-1)); % primary parameter for updating anchor
    c_matrix(j) = (c_matrix(j-1)/(1 + delta)); % Parameter for gamma
    delta = exp(1/(j*j)) - 1; % Needed for both gamma and c, which are in turn needed for anchor update3
    
    % update3 anchor
    anchor(:,j) = anchor(:,j-1) - gamma(j)*A*update3(:,j);
    
end

% % This code is for plotting the algorithm iterations
% plot(update3(1,:),update3(2,:),'LineWidth',2);
% hold on;
% for j=1:N
%     plot(update3(1,j),update3(2,j),'r*');
% end
% plot(anchor(1,:),anchor(2,:),'c');
% for j=1:N
%     plot(anchor(1,j),anchor(2,j),'g*');
% end
% plot(0,0,'k*');
% title('FEG m. anchor, positive \gamma, 2nd fxn')
% set(gca,'linewidth',2)