function [g, h] = nvs20_CNL(x)

% this is the constraint function
% It contains 11 continuous variables and 5 integer variables 
% X(1), X(2), X(3), X(4), X(5) are integer varibales

global MX;
X = x.*MX;

% equations of the model :
e1= 0.22*X(1) + 0.2*X(2) + 0.19*X(3) + 0.25*X(4) + 0.15*X(5) + 0.11*X(6) + 0.12*X(7) + 0.13*X(8) + X(9);
e2= - 1.46*X(1) - 1.3*X(3) + 1.82*X(4) - 1.15*X(5) + 0.8*X(7) + X(10); 
e3= 1.29*X(1) - 0.89*X(2) - 1.16*X(5) - 0.96*X(6) - 0.49*X(8) + X(11); 
e4= - 1.1*X(1) - 1.06*X(2) + 0.95*X(3) - 0.54*X(4) - 1.78*X(6) - 0.41*X(7) + X(12); 
e5= - 1.43*X(4) + 1.51*X(5) + 0.59*X(6) - 0.33*X(7) - 0.43*X(8) + X(13); 
e6= - 1.72*X(2) - 0.33*X(3) + 1.62*X(5) + 1.24*X(6) + 0.21*X(7) - 0.26*X(8) + X(14); 
e7= 1.12*X(1) + 0.31*X(4) + 1.12*X(7) - 0.36*X(9) + X(15);
e8= 0.45*X(2) + 0.26*X(3) - 1.1*X(4) + 0.58*X(5) - 1.03*X(7) + 0.1*X(8) + X(16); 

% inequality constraints :
global LG UG MG;

g(1) = (LG(1) - e1)/MG(1);
g(2) = (e1 - UG(1))/MG(1);
g(3) = (LG(2) - e2)/MG(2);
g(4) = (e2 - UG(2))/MG(2);
g(5) = (LG(3) - e3)/MG(3);
g(6) = (e3 - UG(3))/MG(3);
g(7) = (LG(4) - e4)/MG(4);
g(8) = (e4 - UG(4))/MG(4);
g(9) = (LG(5) - e5)/MG(5);
g(10) = (e5 - UG(5))/MG(5);
g(11) = (LG(6) - e6)/MG(6);
g(12) = (e6 - UG(6))/MG(6);
g(13) = (LG(7) - e7)/MG(7);
g(14) = (e7 - UG(7))/MG(7);
g(15) = (LG(8) - e8)/MG(8);
g(16) = (e8 - UG(8))/MG(8);

h = [];

% call number counter for nb_CNL_eval :
global nb_CNL_eval raz_nb_CNL_eval ;
persistent cpt;
if isempty(cpt) || raz_nb_CNL_eval == 1
    cpt = 1;
    raz_nb_CNL_eval = 0;
else
    cpt = cpt + 1;
end
nb_CNL_eval = cpt;

% call number counter for nb_CNL_eval1 :
global nb_CNL_eval1 raz_nb_CNL_eval1;
persistent cpt1;
if isempty(cpt1) || raz_nb_CNL_eval1 == 1
    cpt1 = 1;
    raz_nb_CNL_eval1 = 0;
else
    cpt1 = cpt1 + 1;
end
nb_CNL_eval1 = cpt1;