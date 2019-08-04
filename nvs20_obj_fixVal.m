function xy = nvs20_obj_fixVal(x)

% This is the objective function 
% It contains 11 continuous variables and 5 integer variables 
% X(1), X(2), X(3), X(4), X(5) are integer varibales

global MX;
X = x.*MX;

xy = 1;

