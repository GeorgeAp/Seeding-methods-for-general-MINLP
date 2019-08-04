function [X, nb_GHgrad_eval] = Xvalide_nvs20_c()
% generation of feasible individuals for MAPSE
% method [Chootinan2006]
% see the file:
% The method of generation of individuals at the frontier of the feasible space_NS.docx
% attention: 16 variables only (excluded), not scaled
% gradients calculated by function call
%
% Difference from Xvalide_nvs20: the number of calls to the
% calculating function G, H and their gradients is counted
% (variable nb_GHgrad_eval)

% variables Xs : 
% X(1), X(2), X(3), X(4), X(5) are integer varibales
global LXs UXs;

global LG3 UG3 C;

I = max(length(LXs),length(UXs)); % I = 16
J = max(length(LG3),length(UG3)); % J = 8
% K = length(C); % K = 2
size_contraint = I + J;
tol = 1.0e-6;

global kr Ma bfer deltap rhocu; % constantes du modèle

% génération d'un individu initial Xs entre les bornes UXs et LXs
Xs = zeros([1 I]);
for i=1:I
    Xs(i) = rand*(UXs(i)-LXs(i))+LXs(i);
end

nb_GHgrad_eval = 0;
max_iter = 100;
for iter = 1:max_iter
    [G, H, DG, DH] = gradGH2_nvs20(Xs);
    nb_GHgrad_eval = nb_GHgrad_eval + 1;
    
    e1      = G(1);
    e2      = G(2);
    e3      = G(3);
    e4      = G(4);
    e5      = G(5);
    e6      = G(6);  
    e7      = G(7);
    e8      = G(8);    
    
    deltaV = [min(0,UXs(1)-Xs(1))+max(0,LXs(1)-Xs(1));...
              min(0,UXs(2)-Xs(2))+max(0,LXs(2)-Xs(2));...
              min(0,UXs(3)-Xs(3))+max(0,LXs(3)-Xs(3));...
              min(0,UXs(4)-Xs(4))+max(0,LXs(4)-Xs(4));...
              min(0,UXs(5)-Xs(5))+max(0,LXs(5)-Xs(5));...
              min(0,UXs(6)-Xs(6))+max(0,LXs(6)-Xs(6));...
              min(0,UXs(7)-Xs(7))+max(0,LXs(7)-Xs(7));...
              min(0,UXs(8)-Xs(8))+max(0,LXs(8)-Xs(8));...
              min(0,UXs(9)-Xs(9))+max(0,LXs(9)-Xs(9));...
              min(0,UXs(10)-Xs(10))+max(0,LXs(10)-Xs(10));...
              min(0,UXs(11)-Xs(11))+max(0,LXs(11)-Xs(11));...
              min(0,UXs(12)-Xs(12))+max(0,LXs(12)-Xs(12));...
              min(0,UXs(13)-Xs(13))+max(0,LXs(13)-Xs(13));...
              min(0,UXs(14)-Xs(14))+max(0,LXs(14)-Xs(14));...
              min(0,UXs(15)-Xs(15))+max(0,LXs(15)-Xs(15));...
              min(0,UXs(16)-Xs(16))+max(0,LXs(16)-Xs(16));...
              min(0,UG3(1)-e1)+max(0,LG3(1)-e1);...
              min(0,UG3(2)-e2)+max(0,LG3(2)-e2);...
              min(0,UG3(3)-e3)+max(0,LG3(3)-e3);...
              min(0,UG3(4)-e4)+max(0,LG3(4)-e4);...
              min(0,UG3(5)-e5)+max(0,LG3(5)-e5);...
              min(0,UG3(6)-e6)+max(0,LG3(6)-e6);...
              min(0,UG3(7)-e7)+max(0,LG3(7)-e7);...
              min(0,UG3(8)-e8)+max(0,LG3(8)-e8);]; % deltaV non compressé

    % test de violation et compression
    continuer = 0;
    % flag_cpr : flag de compression :
    % flag_cpr = 1 si ligne à conserver, = 0 si ligne à éliminer
    flag_cpr = zeros(size_contraint,1); 
    deltaVcpr = [];
    for icomp = 1:size_contraint
        if abs(deltaV(icomp))>tol
            continuer = 1;
            flag_cpr(icomp) = 1;
            deltaVcpr = [deltaVcpr;deltaV(icomp)]; % deltaV compressé
        end
    end
    if continuer == 0         % solution found
        X = Xs;
        return
    end
    
    gradV = [eye(I);DG;DH];
    
    % other calculations
    gradVcpr = [];
    for icomp = 1:size_contraint % compression of gradV
        if flag_cpr(icomp) == 1
            gradVcpr = [gradVcpr;gradV(icomp,:)]; % compression gradV 
        end
    end
    invgradVcpr = pinv(gradVcpr);
    deltaX = invgradVcpr*deltaVcpr;
    Xs = Xs+deltaX';
    if ~isreal(Xs)             % Erreur : Xs : complexe
        X = [];
        return
    end
    
end

% Nb.maxi d'itérations atteint
X = [];
return
