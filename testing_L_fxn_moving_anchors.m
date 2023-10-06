% clf;clc

% Testing simple convex-concave function example
% L(x,y) = \epsilon(||x||^2)/2 + <x,y> - \epsilon(||y||^2)/2
% x,y real numbers.
% Instances to check:
% Stationary anchor (Dr. Ryu et al);
% Mobile anchor with positive \gamma term;
% Mobile anchor with negative \gamma term.
% See associated README file for further info

%% Moving anchor, positive/negative gamma

N = 5000; % number of iterations, choose 2k for visualization, 5k for grad norm squared figure
update3 = zeros(2,N); % Matrix where we store x and y values as we iterate
                     % through the algorithm; will also be used for
                     % plotting
                     
anchor = zeros(2,N); % Matrix for anchor update3s.

eps = 0.01; % Can play with this if desired
% Initial point; note that the anchor and update3 start in the same location
update3(1,1) = 1;
update3(2,1) = 0;
anchor(1,1) = update3(1,1);
anchor(2,1) = update3(2,1);

% Alpha, beta, gamma, and R (R depends on epsilon)
alpha = zeros(1,N);
gamma = zeros(1,N);
c_matrix = zeros(1,N);
delta = (exp(1)-1);
R = sqrt(2 + 2*eps);
alpha(1) = 3/8; % Need to specify, is related to eps.
c_matrix(1) = 20; % Check notes for reason for this specific parameter
% Note that it seems the first entry of gamma is immaterial; we don't need
% it.

for j=2:N
    % update3 primary terms
    update3(1,j) = update3(1,j-1) + (1/(j+1))*(anchor(1,j-1) - update3(1,j-1)) - alpha(j-1)*(eps*(update3(1,j-1) + (1/(j+1))*(anchor(1,j-1) - update3(1,j-1)) - alpha(j-1)*(eps*update3(1,j-1) + update3(2,j-1))) + (update3(2,j-1) + (1/(j+1))*(anchor(2,j-1) - update3(2,j-1)) - alpha(j-1)*(eps*update3(2,j-1) - update3(1,j-1)))); % WHEW finally got the parenthesis right
    update3(2,j) = update3(2,j-1) + (1/(j+1))*(anchor(2,j-1) - update3(2,j-1)) - alpha(j-1)*(eps*(update3(2,j-1) + (1/(j+1))*(anchor(2,j-1) - update3(2,j-1)) - alpha(j-1)*(eps*update3(2,j-1) - update3(1,j-1))) - (update3(1,j-1) + (1/(j+1))*(anchor(1,j-1) - update3(1,j-1)) - alpha(j-1)*(eps*update3(1,j-1) + update3(2,j-1))));
    
    % update3 some misc terms and constants
    % For regular update3
    alpha(j) = alpha(j-1)*(1 - (1/(j)*(j+2))*((alpha(j-1)*alpha(j-1)*R*R)/(1 - alpha(j-1)*alpha(j-1)*R*R)));
    % For anchor
    
    gamma(j) = -((j+1)*delta)/(c_matrix(j-1)); % primary parameter for updating anchor
    c_matrix(j) = (c_matrix(j-1)/(1 + delta)); % Parameter for gamma
    delta = exp(1/(j*j)) - 1; % Needed for both gamma and c, which are in turn needed for anchor update3
    
    
    % update3 the anchor
    anchor(1,j) = anchor(1,j-1) + gamma(j)*(eps*update3(1,j) + update3(2,j));
    anchor(2,j) = anchor(2,j-1) + gamma(j)*(eps*update3(2,j) - update3(1,j));
end

% % This code is for plotting the gradient norm squared
% plot(log(1:N), 2*log((eps*update3(1,:) + update3(2,:)).^2+(eps*update3(2,:) - update3(1,:)).^2))
% title('Gradient norm squared, \gamma positive, moving anchor: EAG setting, \epsilon=0.01')

% % This code is for plotting the algorithm's error
% plot(log(1:N),2*log(sqrt(update3(1,:).^2 + update3(2,:).^2)),'r');
% title('Error, moving anchor, negative \gamma, conv-conc., \epsilon=0.01') % update3 with each running

% % This code is for plotting the algorithm iterations with anchor
% iterations
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
% title('Negative \gamma, \epsilon=0')
% set(gca,'linewidth',2)