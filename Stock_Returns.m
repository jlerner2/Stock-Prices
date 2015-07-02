function [trace_min ,nuc_min] = Stock_Returns(  )
% Author: Jeremy Lerner, Stony Brook University

% Calculate the daily returns, the difference in price from the current day
% and the previous day. Where sp(i).price(j) is the price of stock i 
% on day j
tic
load('test_data.mat')
% returns = zeros(numel(sp),numel(sp(1).price));
total_stocks =  100;%numel(sp); % "= numel(sp)" for all stocks
start_stocks = 100;
stop_stocks = start_stocks + total_stocks - 1 ;
start_time = 2000;
total_time =  1000;%numel(sp(1).price); % "= numel(sp(1).price)" for all days
end_time = start_time + total_time - 1;
% total_time = numel(sp(1).price); % "= numel(sp(1).price)" for all days
returns = zeros(total_stocks,total_time);
for i = start_stocks:stop_stocks
    for j = start_time:end_time
        returns(i-start_stocks+1,j) = (sp(i).price(j)  - sp(i).price(j-1))/sp(i).price(j-1); 
    end    
end
returns = log(1+returns);
Sigma = returns*returns';
% Sigma = rand(100, 10)*rand(10,100) + 0.1*diag(rand(100,1));
% Sigma = Sigma*Sigma';
[m,n] = size(Sigma);

% epsilon = 5e-3;
% W = inv(Sigma+epsilon*eye(n,n));
% W_old = ones(n,n);
% D_old = ones(n,1);
% D = zeros(n,1);
% trace_min_old = 0;
% trace_min = 1;

%  norm(W_old - W)
%  norm(D_old - D)
%  abs(trace_min - trace_min_old)
% 
% while ( norm(W_old - W) > 1e-2 && norm(D_old - D) > 1e-2 && abs(trace_min - trace_min_old) > 1e-2 )
%     tic
%     W_old = W;
%     D_old = D;
%     trace_min_old = trace_min;
%     
%     % Use CVX to optimize the trace, with variables Sigma_h (n x n matrix)
%     % and D (n x 1 diagonal matrix)
%     cvx_begin
%         variables D(n) Sigma_h(n,n)
%         minimize trace(W*Sigma_h)
%         subject to
%             %This is equivalent to D>=0, as D is represented as a vector
%             diag(D) == semidefinite(n);
%             Sigma_h == semidefinite(n);
%             Sigma == Sigma_h + diag(D);
%     cvx_end
%     trace_min =  trace(W*Sigma_h);
%     W = inv(Sigma - diag(D) + epsilon*eye(n,n));
%     toc
%     Sigma_hat = Sigma_h;
% end
% 
% 
% 
toc
% The first and only weight is (Sigma)^-1
W = inv(Sigma);

cvx_begin
    variables D(n) Sigma_h(n,n)
    minimize trace(W*(Sigma_h))
    subject to
        diag(D) == semidefinite(n);
        Sigma_h == semidefinite(n);
        Sigma == Sigma_h + diag(D);
cvx_end
trace_min = trace(W*Sigma_h);


% 
% Solve the convex optimization problem in CVX
cvx_begin
    variables Sigma_hat(m,n) Y(m,m) Z(n,n) D(n)
    minimize trace(Y)/2 + trace(Z)/2
    subject to
        [Y W*Sigma_hat; (W*Sigma_hat)' Z] == semidefinite(n+m);   
        Sigma == Sigma_hat + diag(D);
        Sigma_hat == semidefinite(n);
        diag(D) == semidefinite(n);
cvx_end
nuc_min = trace(Y)/2 + trace(Z)/2;








% epsilon = 5e-3;
% W = inv(Sigma+epsilon*eye(n,n));
% W_old = ones(n,n);
% D_old = ones(n,1);
% D = zeros(n,1);
% trace_min_old = 0;
% trace_min = 1;
% 
% 
% while ( norm(W_old - W) > 1e-2 && norm(D_old - D) > 1e-2 && abs(trace_min - trace_min_old) > 1e-2 )
%     W_old = W;
%     D_old = D;
%     trace_min_old = trace_min;
%     
%     % Use CVX to optimize the trace, with variables Sigma_h (n x n matrix)
%     % and D (n x 1 diagonal matrix)
%     cvx_begin
%         variables Sigma_hat(m,n) Y(m,m) Z(n,n) D(n)
%         minimize trace(Y)/2 + trace(Z)/2
%         subject to
%             [Y W*Sigma_hat; (W*Sigma_hat)' Z] == semidefinite(n+m);   
%             Sigma == Sigma_hat + diag(D);
%             Sigma_hat == semidefinite(n);
%             diag(D) == semidefinite(n);
%     cvx_end
%     trace_min = trace(Y)/2 + trace(Z)/2;
%     W = inv(Sigma - diag(D) + epsilon*eye(n,n));
% end
toc
% % 
% % % stocks 100-299: 168.2668
% % toc