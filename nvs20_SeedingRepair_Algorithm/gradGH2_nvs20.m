function [G, H, DG, DH] = gradGH2_nvs20(Xs)
% constraints and their gradients
% all calculated with analytic formula

global kr Ma bfer deltap rhocu; % constantes du modèle

% calculate de G

% calculate e1,...,e8 :
e1= 0.22*Xs(1) + 0.2*Xs(2) + 0.19*Xs(3) + 0.25*Xs(4) + 0.15*Xs(5) + 0.11*Xs(6) + 0.12*Xs(7) + 0.13*Xs(8) + Xs(9);
e2= - 1.46*Xs(1) - 1.3*Xs(3) + 1.82*Xs(4) - 1.15*Xs(5) + 0.8*Xs(7) + Xs(10); 
e3= 1.29*Xs(1) - 0.89*Xs(2) - 1.16*Xs(5) - 0.96*Xs(6) - 0.49*Xs(8) + Xs(11); 
e4= - 1.1*Xs(1) - 1.06*Xs(2) + 0.95*Xs(3) - 0.54*Xs(4) - 1.78*Xs(6) - 0.41*Xs(7) + Xs(12); 
e5= - 1.43*Xs(4) + 1.51*Xs(5) + 0.59*Xs(6) - 0.33*Xs(7) - 0.43*Xs(8) + Xs(13); 
e6= - 1.72*Xs(2) - 0.33*Xs(3) + 1.62*Xs(5) + 1.24*Xs(6) + 0.21*Xs(7) - 0.26*Xs(8) + Xs(14); 
e7= 1.12*Xs(1) + 0.31*Xs(4) + 1.12*Xs(7) - 0.36*Xs(9) + Xs(15);
e8= 0.45*Xs(2) + 0.26*Xs(3) - 1.1*Xs(4) + 0.58*Xs(5) - 1.03*Xs(7) + 0.1*Xs(8) + Xs(16); 

G(1) = e1;
G(2) = e2;
G(3) = e3;
G(4) = e4;
G(5) = e5;
G(6) = e6;
G(7) = e7;
G(8) = e8;

H=[];

% calcul du gradient de G

%    1/ partial derivative of e1
DG(1,1) = 0.22; %de1/dx1
DG(1,2) = 0.2; %de1/dx2
DG(1,3) = 0.19; %de1/dx3
DG(1,4) = 0.25; %de1/dx4
DG(1,5) = 0.15; %de1/dx5
DG(1,6) = 0.11; %de1/dx6
DG(1,7) = 0.12; %de1/dx7
DG(1,8) = 0.13; %de1/dx8
DG(1,9) = 1; %de1/dx9
DG(1,10) = 0; %de1/dx10
DG(1,11) = 0; %de1/dx11
DG(1,12) = 0; %de1/dx12
DG(1,13) = 0; %de1/dx13
DG(1,14) = 0; %de1/dx14
DG(1,15) = 0; %de1/dx15
DG(1,16) = 0; %de1/dx16

%    2/ partial derivative of e2
DG(2,1) = -1.46; %de2/dx1
DG(2,2) = 0; %de2/dx2
DG(2,3) = -1.3; %de2/dx3
DG(2,4) = 1.82; %de2/dx4
DG(2,5) = -1.15; %de2/dx5
DG(2,6) = 0; %de2/dx6
DG(2,7) = 0; %de2/dx7
DG(2,8) = 0; %de2/dx8
DG(2,9) = 0; %de2/dx9
DG(2,10) = 1; %de2/dx10
DG(2,11) = 0; %de2/dx11
DG(2,12) = 0; %de2/dx12
DG(2,13) = 0; %de2/dx13
DG(2,14) = 0; %de2/dx14
DG(2,15) = 0; %de2/dx15
DG(2,16) = 0; %de2/dx16

%    3/ partial derivative of e3
DG(3,1) = 1.29; %de3/dx1
DG(3,2) = -0.89; %de3/dx2
DG(3,3) = 0; %de3/dx3
DG(3,4) = 0; %de3/dx4
DG(3,5) = -1.16; %de3/dx5
DG(3,6) = -0.96; %de3/dx6
DG(3,7) = 0; %de3/dx7
DG(3,8) = -0.49; %de3/dx8
DG(3,9) = 0; %de3/dx9
DG(3,10) = 0; %de3/dx10
DG(3,11) = 1; %de3/dx11
DG(3,12) = 0; %de3/dx12
DG(3,13) = 0; %de3/dx13
DG(3,14) = 0; %de3/dx14
DG(3,15) = 0; %de3/dx15
DG(3,16) = 0; %de3/dx16

%    4/ partial derivative of e4
DG(4,1) = -1.1; %de4/dx1
DG(4,2) = -1.06; %de4/dx2
DG(4,3) = 0.95; %de4/dx3
DG(4,4) = -0.54; %de4/dx4
DG(4,5) = 0; %de4/dx5
DG(4,6) = -1.78; %de4/dx6
DG(4,7) = -0.41; %de4/dx7
DG(4,8) = 0.13; %de4/dx8
DG(4,9) = 1; %de4/dx9
DG(4,10) = 0; %de4/dx10
DG(4,11) = 0; %de4/dx11
DG(4,12) = 1; %de4/dx12
DG(4,13) = 0; %de4/dx13
DG(4,14) = 0; %de4/dx14
DG(4,15) = 0; %de4/dx15
DG(4,16) = 0; %de4/dx16

%    5/ partial derivative of e5
DG(5,1) = 0; %de5/dx1
DG(5,2) = 0; %de5/dx2
DG(5,3) = 0; %de5/dx3
DG(5,4) = -1.43; %de5/dx4
DG(5,5) = 1.51; %de5/dx5
DG(5,6) = 0.59; %de5/dx6
DG(5,7) = -0.33; %de5/dx7
DG(5,8) = -0.43; %de5/dx8
DG(5,9) = 0; %de5/dx9
DG(5,10) = 0; %de5/dx10
DG(5,11) = 0; %de5/dx11
DG(5,12) = 0; %de5/dx12
DG(5,13) = 1; %de5/dx13
DG(5,14) = 0; %de5/dx14
DG(5,15) = 0; %de5/dx15
DG(5,16) = 0; %de5/dx16

%    6/ partial derivative of e6
DG(6,1) = 0; %de6/dx1
DG(6,2) = -1.72; %de6/dx2
DG(6,3) = -0.33; %de6/dx3
DG(6,4) = 0; %de6/dx4
DG(6,5) = 1.62; %de6/dx5
DG(6,6) = 1.24; %de6/dx6
DG(6,7) = 0.21; %de6/dx7
DG(6,8) = -0.26; %de6/dx8
DG(6,9) = 0; %de6/dx9
DG(6,10) = 0; %de6/dx10
DG(6,11) = 0; %de6/dx11
DG(6,12) = 0; %de6/dx12
DG(6,13) = 0; %de6/dx13
DG(6,14) = 1; %de6/dx14
DG(6,15) = 0; %de6/dx15
DG(6,16) = 0; %de6/dx16

%    7/ partial derivative of e7
DG(7,1) = 1.12; %de1/dx1
DG(7,2) = 0; %de1/dx2
DG(7,3) = 0; %de1/dx3
DG(7,4) = 0.31; %de1/dx4
DG(7,5) = 0; %de1/dx5
DG(7,6) = 0; %de1/dx6
DG(7,7) = 1.12; %de1/dx7
DG(7,8) = 0; %de1/dx8
DG(7,9) = -0.36; %de1/dx9
DG(7,10) = 0; %de1/dx10
DG(7,11) = 0; %de1/dx11
DG(7,12) = 0; %de1/dx12
DG(7,13) = 0; %de1/dx13
DG(7,14) = 0; %de1/dx14
DG(7,15) = 1; %de1/dx15
DG(7,16) = 0; %de1/dx16

%    8/ partial derivative of e8
DG(8,1) = 0; %de8/dx1
DG(8,2) = 0.45; %de8/dx2
DG(8,3) = 0.26; %de8/dx3
DG(8,4) = -1.1; %de8/dx4
DG(8,5) = 0.58; %de8/dx5
DG(8,6) = 0; %de8/dx6
DG(8,7) = -1.03; %de8/dx7
DG(8,8) = 0.1; %de8/dx8
DG(8,9) = 0; %de8/dx9
DG(8,10) = 0; %de8/dx10
DG(8,11) = 0; %de8/dx11
DG(8,12) = 0; %de8/dx12
DG(8,13) = 0; %de8/dx13
DG(8,14) = 0; %de8/dx14
DG(8,15) = 0; %de8/dx15
DG(8,16) = 1; %de8/dx16

DH = [];
end

