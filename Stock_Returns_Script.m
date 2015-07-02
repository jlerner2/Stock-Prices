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






stock_pick = oil_stocks;
time =  3524:3926;


trace_min = zeros(4,1);
trace_min_iter = zeros(4,1);
nuc_min = zeros(4,1);
nuc_min_iter = zeros(4,1);

i = 1;

[trace_min(i), trace_min_iter(i), nuc_min(i), nuc_min_iter(i)] = Stock_Returns_Iter_Full_Noise_Const(stock_pick, time  );

num = numel(stock_pick);
i = 2;
[trace_min(i), trace_min_iter(i), nuc_min(i), nuc_min_iter(i)] = Stock_Returns_Iter_Full_Noise_Const((1:num), time  );

i = 3;
[trace_min(i), trace_min_iter(i), nuc_min(i), nuc_min_iter(i)] = Stock_Returns_Iter_Full_Noise_Const((1:num)+100, time  );

i = 4;
[trace_min(i), trace_min_iter(i), nuc_min(i), nuc_min_iter(i)] = Stock_Returns_Iter_Full_Noise_Const((1:num)+200, time  );

A = [ trace_min , trace_min_iter , nuc_min , nuc_min_iter]';

matrix2lyx(A, 'save');

