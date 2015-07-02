function [trace_min, trace_min_iter, nuc_min, nuc_min_iter] = Stock_Returns_Iter_Full_Noise_Const(stock_pick, time  )
% Author: Jeremy Lerner, Stony Brook University
% Note: put the stock pick on line 51, pick from the below stocks. Energy
% stocks seem the most correlated

% Calculate the daily returns, the difference in price from the current day
% and the previous day. Where sp(i).price(j) is the price of stock i 
% on day j
tic
load('test_data.mat')
% returns = zeros(numel(sp),numel(sp(1).price));

% % OLD 37 Technology stocks, questions: 214/Lockheed Martin, 241/motorola,
% % 283/Parker-Hannifin, 335/ ATT, 345/Thermo Fisher Scientific, 368/Waters, 
% tech_stocks = [2 , 7, 8, 11, 21, 23,  24, 63, 91, 92, 96, 103, 124, 172, 176,...
%     180, 183, 211, 214, 217, 218, 227, 236, 240, 241, 243, 258, 266, 283,... 
%     299 , 320, 333, 335, 345, 352, 368,  369];

% 43 Technology stocks, questions: 214/Lockheed Martin, 241/motorola,
% 283/Parker-Hannifin, 335/ ATT, 345/Thermo Fisher Scientific, 368/Waters, 
tech_stocks = [2 , 7, 8, 11, 21, 23,  24, 63, 91, 92, 96, 103, 124, 172, 176,...
    180, 183, 190, 201, 211, 214, 217, 218, 227, 236, 240, 241, 243, 258, 266, 283,... 
    299 , 305,  320, 333, 335,339, 345, 352, 366, 368,  369, 382, ];

% 35 Healthcare stocks
healthcare_stocks = [3,4 , 6,17,  22,25, 43, 47, 48,57, 59, 65, 72, 73, 76, 105, ...
    117, 130, 148, 164, 210, 212, 228, 230, 237, 246, 276, 280, 294,...
    326, 332, 341, 355, 360, 382] ;

% 12 Insurance stocks
insurance_stocks = [15,16,18, 20,  29, 68, 77, 168, 179, 282, 347, 381] ;

% 9 Real estate stocks
real_estate_stocks = [36,61, 106, 128, 165, 200, 286, 322, 364];

% 24 Banking/finance stocks
banking_stocks = [39, 42,45,  50,  53, 62, 86, 141, 142, 143,  163, 196, 199, 239, 242, 259, 270, 302, 312, 325, 350, 358, 371, 388] ;

% % OLD! 39 Energy stocks
% energy_stocks = [ 12, 13, 14, ...
%     120, 131, 134, 140, 150, 160, 167, 238, 244, 247, 248, 249, 252, 253, 256, 264, 268, 273, 277, 291 , 293, 298, ...
%     301 , 308, 311, 317, 321, 330, 337, 338, 349, 362, 370, 375, 380, 383] ;

% 60 Energy stocks
energy_stocks = [ 12, 13, 14, 30, 31, 52, 66, 74, 84, 85, 87, 88, 99, 100, 111, 115, 116, ...
    120, 122, 127, 129, 131, 134, 140, 150, 160, 167, 172, 238, 244, 247, 248, 249, 250, ...
    252, 253, 256, 260, 264, 268, 273, 277, 290, 291 , 293, 298, ...
    301 , 308, 311, 317, 321, 330, 337, 338, 349, 362, 370, 375, 380, 383] ;

%********************

% 12 Insurance stocks
insurance_stocks = [5, 15, 16, 18, 20, 29, 68, 76, 77, 168, 207, 215, 282, 344,  356, 381];
    
% 9 Metals and mining
metal_stocks = [ 1, 19, 35, 79, 137, 251, 261, 275, 379];

% 2 Agriculture stocks
agriculture_stocks = [9, 102];

% 5 Big Retailers
big_retailer_stocks = [ 26, 89, 109, 138, 376];

% 11 Clothing retailer
clothing_retailer_stocks = [ 28, 147, 157, 192, 197, 206, 254, 296, 304, 307, 361];

% 30 Oil stocks
oil_stocks = [ 30, 31, 52, 66, 74, 87, 88, 99, 111, 127, 131, 140, 160, 167, 172, 238, 244, 247, ...
    248, 249, 252, 256, 268, 298, 308, 317, 330, 349, 362, 383];

% 10 Real estate investment trusts
real_estate_investment_stocks = [ 36, 61, 128, 164, 165, 200, 295, 322, 364, 365];

% 2 Cosmetics/beauty
cosmetic_stocks = [ 37, 123];

% 10 Super regional banks
super_regional_banks_stocks = [ 42, 81, 86, 143, 163, 199, 288, 325, 358, 371];

% 8 Auto parts and cars
auto_stocks = [ 27, 40, 60, 169, 191, 203, 267, 272];

% 8 Defense/aerospace stocks
defense_stocks = [ 41, 152, 170, 214, 255, 309, 353, 359]; 

% 4 Liquor stocks
beer_stocks = [49, 51, 328, 336];

% 9 Food stocks
food_stocks = [64, 90, 104, 154, 175, 198, 232, 316, 348];

% 11 Entertainment
entertainment_stocks = [ 69, 71, 82, 97, 108, 118, 151, 162, 182, 225, 351];

% 7 Transportation/Logistics stocks
transport_stocks = [ 75, 93, 135, 139, 257, 300, 357 ];

% Large Cap stocks (S&P 500)
% large_cap = [1, 
%     
% % Mid Cap stocks (S&P Midcap 400)
mid_cap = [24, 28, 33, 37, 56, 141, 190, 192, 193, 248, 301, 318, 340, 379 ];
% 
% % Small Cap stocks (Russell 200)
% small_cap = [ 

%***************************************************************%
% total_stocks =  8;%numel(sp); % "= numel(sp)" for all stocks
% start_stocks = 201;
% stop_stocks = start_stocks + total_stocks - 1 ;
% start_time =  2;
% total_time =  numel(sp(1).price); %1000 "= numel(sp(1).price)" for all days
% end_time = start_time + total_time - 2;
% returns = zeros(total_stocks,total_time-1);
% for i = start_stocks:stop_stocks
%     for j = start_time:end_time
%         returns(i-start_stocks+1,j) = (sp(i).price(j) - sp(i).price(j-1))/sp(i).price(j-1); 
%     end    
% end
%***************************************************************%

% %***************************************************************%
% % Put the stock pick here
% stock_pick = 201:300;
% start_time =  2;%1926;%2;
% total_time =  110;%numel(sp(1).price); %1000 "= numel(sp(1).price)" for all days
% end_time = start_time + total_time - 2;
% returns = zeros(numel(stock_pick),total_time-1);
% k = 1;
% for i = stock_pick
%     for j = start_time:end_time
%         returns(k,j) = (sp(i).price(j) - sp(i).price(j-1))/sp(i).price(j-1); 
%     end    
%     k = k+1;
% end
% %***************************************************************%

%***************************************************************%
% Put the stock pick here
% stock_pick = 201:230;
% time = 3020:3523;%numel(sp(1).date) ; %2:10:3926;
returns = zeros(numel(stock_pick),numel(time));
k = 1;

for i = stock_pick
    inner_var = 1;
    for j = time
        returns(k,inner_var) = (sp(i).price(j) - sp(i).price(j-1))/sp(i).price(j-1); 
        inner_var = inner_var + 1;
    end    
    k = k+1;
end
%***************************************************************%


returns = log(1+returns);
Sigma = returns*returns';

[m,n] = size(Sigma);
const = norm(Sigma,2)*(2/100);

epsilon = 5e-3;
W = inv(Sigma+epsilon*eye(n,n));
W_old = ones(n,n);
D_old = ones(n,n);
D = zeros(n,n);
trace_min_old = 0;
trace_min_iter = 1;
toc

% % Trace iterative
% while ( norm(W_old - W) > 1e-2 && norm(D_old - D) > 1e-2 && abs(trace_min_iter - trace_min_old) > 1e-2 )
%     W_old = W;
%     D_old = D;
%     trace_min_old = trace_min_iter;
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
%     trace_min_iter =  trace(W*Sigma_h);
%     W = inv(Sigma - diag(D) + epsilon*eye(n,n));
%     toc
% end

% Note: this (commented out) version of the iterative trace minimization
% only optimizes over diag(D), instead of over diag(D) AND Sigma_hat, as
% above. The results are very similar (nearly identical)

% epsilon = 5e-3;
% W = inv(Sigma+epsilon*eye(n,n));
% W_old = ones(n,n);
% D_old = ones(n,1);
% D = zeros(n,1);
% trace_min_old = 0;
% trace_min_iter = 1;
% toc
% 
% Trace iterative v2
while ( norm(W_old - W) > 1e-2 && norm(D_old - D) > 1e-2 && abs(trace_min_iter - trace_min_old) > 1e-2 )
    W_old = W;
    D_old = D;
    trace_min_old = trace_min_iter;
    
    % Use CVX to optimize the trace, with variables Sigma_h (n x n matrix)
    % and D (n x 1 diagonal matrix)
    cvx_begin
        variable D(n,m)
        minimize trace( W*(Sigma - D) )
        subject to
            %This is equivalent to D>=0, as D is represented as a vector
            D == semidefinite(n);
            Sigma - D == semidefinite(n);
            norm(D,2) <= const;
    cvx_end
    trace_min_iter =  trace(W*(Sigma - D));
    W = inv(Sigma - D + epsilon*eye(n,n));
    toc
end


W = inv(Sigma+epsilon*eye(n,n));
W_old = ones(n,n);
D_old = ones(n,n);
D = zeros(n,n);
nuc_min_old = 0;
nuc_min_iter = 1;


% %Nuclear norm iterative
% while ( norm(W_old - W) > 1e-2 && norm(D_old - D) > 1e-2 && abs(nuc_min_iter - nuc_min_old) > 1e-2 )
%     W_old = W;
%     D_old = D;
%     nuc_min_old = nuc_min_iter;
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
%     nuc_min_iter = trace(Y)/2 + trace(Z)/2;
%     W = inv(Sigma - diag(D) + epsilon*eye(n,n));
% end

%Nuclear norm iterative
while ( norm(W_old - W) > 1e-2 && norm(D_old - D) > 1e-2 && abs(nuc_min_iter - nuc_min_old) > 1e-2 )
    W_old = W;
    D_old = D;
    nuc_min_old = nuc_min_iter;
    % Use CVX to optimize the trace, with variables Sigma_h (n x n matrix)
    % and D (n x 1 diagonal matrix)
    cvx_begin
        variables Y(m,m) Z(n,n) D(n,n)
        minimize trace(Y)/2 + trace(Z)/2
        subject to
            [Y W*(Sigma - D); (W*(Sigma - D))' Z] == semidefinite(n+m);   
            Sigma - D ==  semidefinite(n);
            D == semidefinite(n);
            norm(D,2) <= const;
    cvx_end
    nuc_min_iter = trace(Y)/2 + trace(Z)/2;
    W = inv(Sigma - D + epsilon*eye(n,n));
end


% % The first and only weight is (Sigma)^-1
% W = inv(Sigma);
% %Trace 
% cvx_begin
%     variables D(n) Sigma_h2(n,n)
%     minimize trace(W*(Sigma_h2))
%     subject to
%         diag(D) == semidefinite(n);
%         Sigma_h2 == semidefinite(n);
%         Sigma == Sigma_h2 + diag(D);
% cvx_end
% trace_min = trace(W*Sigma_h2);

% The first and only weight is (Sigma)^-1
W = inv(Sigma);
%Trace 
cvx_begin
    variables D(n,n)
    minimize trace(W*(Sigma - D))
    subject to
        D == semidefinite(n);
        Sigma - D == semidefinite(n);
        norm(D,2) <= const;
cvx_end
trace_min = trace(W*(Sigma-D));


% 
% % Solve the convex optimization problem in CVX
% % Nuclear norm
% cvx_begin
%     variables Sigma_hat(m,n) Y(m,m) Z(n,n) D(n)
%     minimize trace(Y)/2 + trace(Z)/2
%     subject to
%         [Y W*Sigma_hat; (W*Sigma_hat)' Z] == semidefinite(n+m);   
%         Sigma == Sigma_hat + diag(D);
%         Sigma_hat == semidefinite(n);
%         diag(D) == semidefinite(n);
% cvx_end
% nuc_min = trace(Y)/2 + trace(Z)/2;


% Nuclear norm
cvx_begin
    variables Y(m,m) Z(n,n) D(n,n)
    minimize trace(Y)/2 + trace(Z)/2 
    subject to
        [Y W*(Sigma - D); (W*(Sigma - D))' Z] == semidefinite(n+m);   
        Sigma - D ==  semidefinite(n);
        D == semidefinite(n);
        norm(D,2) <= const;
cvx_end
nuc_min = trace(Y)/2 + trace(Z)/2;


toc
