function [xy] = nvs20_val(x)

% This code calculates the objective function value with the obtained
% solution variables where x is the vector

global MX;
X = x.*MX; 

xy = 46 + 6*X(4) + X(4)*X(15) + X(4)*X(15).^2 + X(4)*X(11) + X(4)*X(11).^2 + X(4)*X(7) + X(4)*X(7).^2 ...
+ X(4)*X(1) + X(4)*X(1).^2 + 7*X(4).^2 + X(4).^2*X(15) + X(4).^2*X(15).^2 + X(4).^2*X(11) + X(4).^2*X(11).^2 ...
+ X(4).^2*X(7) + X(4).^2*X(7).^2 + X(4).^2*X(1) + X(4).^2*X(1).^2 + 2*X(4).^3 + X(4).^4 + 5*X(13) + X(13)*X(11) ...
+ X(13)*X(11).^2 + X(13)*X(14) + X(13)*X(14).^2 + X(13)*X(7) + X(13)*X(7).^2 + 6*X(13).^2 + X(13).^2*X(11) + ...
X(13).^2*X(11).^2 + X(13).^2*X(14) + X(13).^2*X(14).^2 + X(13).^2*X(7) + X(13).^2*X(7).^2 + 2*X(13).^3 + X(13).^4 + 5 *X(12) ...
+ X(12)*X(5) + X(12)*X(5)^2 + X(12)*X(14) + X(12)*X(14)^2 + X(12)*X(9) + X(12)*X(9)^2 + 6*X(12)^2 + X(12)^2*X(5) + ...
X(12).^2*X(5).^2 + X(12).^2*X(14) + X(12).^2*X(14).^2 + X(12).^2*X(9) + X(12).^2*X(9).^2 + 2*X(12).^3 + X(12).^4 + ...
5*X(6) + X(6)*X(5) + X(6)*X(5).^2 + X(6)*X(15) + X(6)*X(15).^2 + X(6)*X(8) + X(6)*X(8).^2 + 6*X(6).^2 + X(6).^2 *X(5) + ...
X(6).^2*X(5).^2 + X(6).^2*X(15) + X(6).^2*X(15).^2 + X(6).^2*X(8) + X(6).^2*X(8).^2 + 2*X(6).^3 + X(6).^4 + 6*X(5) + X(5)*X(10) ...
+ X(5)*X(10).^2 + X(5)*X(16) + X(5)*X(16).^2 + 7*X(5).^2 + X(5).^2*X(10) + X(5).^2*X(10).^2 +X(5).^2*X(16) + X(5).^2*X(16).^2 + ...
2*X(5).^3 + X(5).^4 + 5*X(15) + X(15)*X(8) + X(15)*X(8).^2 + 6*X(15).^2 + X(15).^2*X(8) + X(15).^2*X(8).^2+ 2*X(15).^3 + X(15).^4 ...
+ 5*X(11) + X(11)*X(7) + X(11)*X(7).^2 + 6*X(11).^2 + X(11).^2*X(7) + X(11).^2*X(7).^2 + 2*X(11).^3 + X(11).^4 + 6*X(14) + X(14)*X(10) ...
+ X(14)*X(10).^2 + X(14)*X(3) + X(14)*X(3).^2 + 7*X(14).^2 + X(14).^2*X(10) + X(14).^2*X(10).^2 + X(14).^2*X(3) + X(14).^2*X(3).^2 ...
+ 2*X(14).^3 + X(14).^4 + 5*X(9) + X(9)*X(3) + X(9)*X(3).^2 + X(9)*X(16) + X(9)*X(16).^2 + 6*X(9).^2 + X(9).^2*X(3) + X(9).^2*X(3).^2 ...
+ X(9).^2*X(16) + X(9).^2*X(16).^2 + 2*X(9).^3 + X(9).^4 + 7*X(10) + X(10)*X(3) + X(10)*X(3).^2 + X(10)*X(2) + X(10)*X(2).^2 + X(10)*X(8) ...
+ X(10)*X(8).^2 + 8*X(10).^2 + X(10).^2*X(3) + X(10).^2*X(3).^2 + X(10).^2*X(2) + X(10).^2*X(2).^2 + X(10).^2*X(8) + X(10).^2*X(8).^2 + ...
2*X(10)^3 + X(10)^4 + 7*X(3) + X(3)*X(2) + X(3)*X(2)^2 + X(3)*X(7) + X(3)*X(7)^2 + 8*X(3)^2 + X(3)^2*X(2) + X(3)^2*X(2)^2 + ...
X(3).^2*X(7) + X(3).^2*X(7).^2 + 2*X(3).^3 + X(3).^4 + 5*X(2) + X(2)*X(7) + X(2)*X(7).^2 + 6*X(2).^2 + X(2).^2*X(7) + X(2).^2*X(7).^2 + 2* ...
X(2).^3 + X(2).^4 + 5*X(16) + X(16)*X(1) + X(16)*X(1).^2 + 6*X(16).^2 + X(16).^2*X(1) + X(16).^2*X(1).^2 + 2*X(16).^3 + X(16).^4 + 6*X(8) ...
+ X(8)*X(1) + X(8)*X(1).^2 + 7*X(8).^2 + X(8).^2*X(1) + X(8).^2*X(1).^2 + 2*X(8).^3 + X(8).^4 + 8*X(7) + X(7)*X(1) + X(7)*X(1).^2 + ...
9*X(7).^2 + X(7).^2*X(1) + X(7).^2*X(1).^2 + 2*X(7).^3 + X(7).^4 + 6*X(1) + 7*X(1).^2 + 2*X(1).^3 + X(1).^4;

