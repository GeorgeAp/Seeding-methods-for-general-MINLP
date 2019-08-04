% script : " ga_seeding_chootinan_mapseVu_s3cx "
% 04/12/2013
% Lance Nexec ex�cutions pour minimiser Vu (MAPSE) en 3 phases :
%    1- seeding par m�thode "chootinan"
%    2- ga
%    3- sqp avec " fmincon "
% A chaque phase :
%    - comptage de toute fonction "de base" appel�e
%    - chronom�trage
% Ecriture : dans un fichier xlsx

%%%%%%%%%% SECTION A DOCUMENTER AVANT EXECUTION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; 

PC = 'PC-GSCOP-288';
problem = 'nvs20';

v = ver;
matlab_version = regexprep(v(1,1).Release, '(','');
matlab_version = {regexprep(matlab_version, ')','')};

% change the integer (last argument) :
RandStream.setGlobalStream(RandStream('mcg16807', 'Seed', 5))
 
Npop = 20;  % number of individuals in GA
tolfun_ga = 1e-6;
tolcon_ga = 1e-6;

tolfun_naive_ga = 1e-3;
tolcon_naive_ga = 1e-3;

gen_ext_ga_sqp = 100;
stall_gen_limit_ga = 50;

aaaammjj = datestr(now,'yyyymmdd');
fich_res = ['res_seedingChootinan_', problem, '_s3cx_pop', num2str(Npop), ...
    '_gen', num2str(gen_ext_ga_sqp),'_stall', num2str(stall_gen_limit_ga), ...
    '_', aaaammjj, '.xlsx'];

list = ls(fich_res);
if ~isempty(list)
    disp('result file: already exists: can not continue');
    return
end

I = 16; % number of variables

% terminals LX, UX, and reduced terminals lx, ux, on input var. :
global MX LX UX LXs UXs;
LX = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; 
UX = [200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200];
MX = (LX + UX)*0.5; % average of the terminals
lx = LX./MX; % (lb_r)
ux = UX./MX; % (ub_r)
LXs=LX;
UXs=LX;

% bounds of nonlinear inequality constraints
% and constant C non-linear equality constraints:
global MG LG UG LG3 UG3;
LG = [2.5 1.1 -3.1 -3.5 1.3 2.1 2.3 -1.5];
UG = [100 100 100 100 100 100 100 100];
MG = (LG + UG)*0.5; % average of the bounds
LG3 = LG;
UG3 = UG;

tol_eq = 1.0e-4; % tol sur les contraintes d'�galit�
handle_CNL = @nvs20_CNL; % constraint function
handle_obj = @nvs20_obj; % objective function
tous_xi = zeros(Npop,I); % Npop initial individuals (1 per line)

nb_valides_ga = 0;
resu_br = [];
resu_val = [];
global nb_CNL_eval raz_nb_CNL_eval;

% boucle for sur les Nexec ex�cutions
Nexec = 2;
for iexec = 1:Nexec
    tic; % lance chrono seeding
    
    % GENERATION DE LA POPULATION INITIALE
    nb_grad_eval_sd = 0;
    tous_xi = [];
    for ipop = 1:Npop
        Xi = [];
        while isempty(Xi)
            [Xi, ngr] = Xvalide_nvs20_c();
            nb_grad_eval_sd = nb_grad_eval_sd + ngr;
        end
        xi = Xi./MX;
        tous_xi = [tous_xi;xi];
    end
    % fin de la g�n�ration de la pop.initiale
    
    tElapsed_sd = toc; % dur�e du seeding
    nb_obj_eval_sd = 0; % pas d'appel � la fct.calcul de vu
    nb_CNL_eval_sd = 0; % pas d'appel � la fct.calcul des CNL
    resu1 = [iexec, tElapsed_sd, nb_grad_eval_sd, nb_obj_eval_sd, nb_CNL_eval_sd];
    
    % �criture de la structure "options_ga"
    options_ga = [];
    options_ga = gaoptimset(@ga);
    options_ga.InitialPopulation = tous_xi;
    options_ga.Display = 'off';
    options_ga.MigrationFraction = 0.0; % pas de migration
    options_ga.PopulationSize = Npop;
    options_ga.EliteCount = 2;
    options_ga.CrossoverFraction = 0.8;
    options_ga.Generations = 400;
    options_ga.StallGenLimit = 200;
    options_ga.TolFun = 1.0e-6;
    options_ga.TolCon = 1.0e-6;
    options_ga.InitialPenalty = 10;
    options_ga.PenaltyFactor = 100;
    options_ga.FitnessScalingFcn = @fitscalingrank; % fonction de scaling
    options_ga.SelectionFcn = @selectionstochunif; % fct.de s�lection
    options_ga.CrossoverFcn = @crossoverscattered; % fct.de croisement
    options_ga.MutationFcn = @mutationadaptfeasible; % fonction de mutation
    % fin de l'�criture de la structure "options"
    initoptions_ga = options_ga; % m�morisation de la structure options

    % lancement de ga :
    raz_nb_CNL_eval = 1;
    tic; % lance chrono ga
    [x_ga, fval_ga, exitflag_ga, output_ga] = ga(handle_obj, I, ...
        [], [], [], [], lx, ux, handle_CNL, options_ga);
    tElapsed_ga = toc; % dur�e de ga

    % teste si ga a op�r� des changements dans la structure options_ga
    ok_ga = isequal(initoptions_ga,options_ga); % = 1 ssi structures identiques
    if ~ok_ga
       disp('***la structure options_ga a �t� chang�e par l''ex�cution de ga***')
    end

    gamesslength = length(regexprep(output_ga.message,'\n',''));
    resu2 = [exitflag_ga, output_ga.generations, output_ga.funccount, ...
        nb_CNL_eval, gamesslength, output_ga.maxconstraint, ...
        tElapsed_ga, ok_ga];

    % d�termination compl�te de la solution :
%     [vu_ga, va_ga, pj_ga, sigmaem_ga, ech_ga, kf_ga, be_ga, p_ga, c_ga] = mapse_3val_sd2(x_ga);
    vu_ga = nvs20_CNL(x_ga);
    X_ga = x_ga.*MX; % conversion var.r�duites --> var.non r�duites
    resu3 = [X_ga, vu_ga];
    
    % algo de descente locale (post-traitement de la sol.ga)
    options_sqp = [];
    options_sqp = optimoptions('fmincon');
    options_sqp.Algorithm = 'sqp';
    options_sqp.Display = 'off';
    options_sqp.MaxFunEvals = 100000;
    options_sqp.MaxIter = 100000;
    options_sqp.TolCon = 1.0e-6;
    options_sqp.TolFun = 1.0e-6;
    options_sqp.TolX =   1.0e-6;
    %options_sqp.TypicalX = x_ga; % peut causer une erreur
    options_sqp.TypicalX = ones(1,7);
    initoptions_sqp = options_sqp; % m�mo.de la structure options_sqp

    % lancement de sqp :
    raz_nb_CNL_eval = 1;
    tic; % lance chrono sqp
    [x_sqp, fval_sqp, exitflag_sqp, output_sqp] = fmincon(handle_obj, ...
        x_ga, [], [], [], [], lx, ux, handle_CNL, options_sqp);
    tElapsed_sqp = toc; % dur�e de sqp

    % teste si fmincon a op�r� des changements dans la struct.options_sqp
    ok_sqp = isequal(initoptions_sqp,options_sqp); % = 1 ssi struct.identiques
    if ~ok_sqp
       warning('la structure options_sqp a �t� chang�e par l''ex�cution de fmincon')
    end

    sqpmesslength = length(regexprep(output_sqp.message,'\n',''));
    resu4 = [ exitflag_sqp, ...
        output_sqp.funcCount, nb_CNL_eval, sqpmesslength, ...
        output_sqp.constrviolation, tElapsed_sqp, ok_sqp];
    
    % d�termination compl�te de la solution :
%     [vu_sqp, va_sqp, pj_sqp, sigmaem_sqp, ech_sqp, kf_sqp, be_sqp, p_sqp, c_sqp] = mapse_3val_sd2(x_sqp);
    vu_sqp = nvs20_CNL(x_sqp);
    X_sqp = x_sqp.*MX; % conversion var.r�duites --> var.non r�duites
    resu5 = [ X_sqp, vu_sqp];

    
    % �valuation des �carts entre sol. ga et sqp :
    %     (formules v�rifi�es 20131102)
    ec_geno = norm(X_sqp - X_ga); % �cart absolu en X : g�notype
    ec_rel_geno = ec_geno/norm(X_sqp); % �cart relatif en X : g�notype
    ec_geno_reduit = norm(x_sqp - x_ga); % �cart absolu en x (var.r�duites) : g�notype
    ec_rel_geno_reduit = ec_geno_reduit/norm(x_sqp); % �cart relatif en x (var.r�duites) : g�notype
    ec_p = p_sqp - p_ga;% �cart absolu en p (montre que sqp ne fait qu'arrondir p)
    ec_pheno = vu_sqp - vu_ga; % �cart absolu en vu : ph�notype
    ec_rel_pheno = ec_pheno/vu_sqp; % �cart relatif en vu : ph�notype
    
    duree_totale = tElapsed_sd + tElapsed_ga + tElapsed_sqp;
    resu6 = [ec_geno_reduit, ec_rel_geno_reduit, ec_p, ec_pheno, ec_rel_pheno, duree_totale];
    resu = [resu1, resu2, resu3, resu4, resu5, resu6];
    resu_br = [resu_br; resu];
    if exitflag_ga > 0
        resu_val = [resu_val; resu];
        nb_valides_ga = nb_valides_ga + 1;
    end
    
    fprintf(1, ...
        'iexec = %d ; generations ga = %d ; exitflag_ga = %d ; exitflag_sqp = %d ; \n', ...
        iexec, output_ga.generations, exitflag_ga, exitflag_sqp);
    fprintf(1, ...
        'tElapsed_ga = %g ; tElapsed_sqp = %g ; ec_geno_reduit = %d ; ec_pheno = %d . \n\n', ...
        tElapsed_ga, tElapsed_sqp, ec_geno_reduit, ec_pheno);
end

warning off MATLAB:xlswrite:AddSheet

% �criture du fichier r�sultat au format .xlsx : feuille 'param'
noms_params = {'param�tres ga :'; 'TolFun'; 'TolCon'; 'population'; ...
    'g�n�rations'; 'stall gen.';''; 'param�tres sqp :'; 'TolCon'; ...
    'TolFun'; 'TolX'; 'MaxFunEvals'; 'MaxIter'; ''; 'Matlab'};
xlswrite(fich_res, noms_params, 'param', 'A1');
val_param_ga = [options_ga3step.TolFun; options_ga3step.TolCon; Npop; ...
    options_ga3step.Generations; options_ga3step.StallGenLimit];
xlswrite(fich_res, val_param_ga, 'param', 'B2');
val_param_sqp = [options_sqp.TolCon; options_sqp.TolFun; ...
    options_sqp.TolX; options_sqp.MaxFunEvals; options_sqp.MaxIter];
xlswrite(fich_res, val_param_sqp, 'param', 'B9');
xlswrite(fich_res, matlab_version, 'param', 'b15');
xlswrite(fich_res, {aaaammjj}, 'param', 'b17');

% �criture du fichier r�sultat au format .xlsx : feuille 'bruts'
entete = {'iexec','tElapsed_sd','nb_grad_eval_sd','nb_obj_eval_sd', ...
    'nb_CNL_eval_sd','exitflag_ga','output_ga.generations', ...
    'output_ga.funccount','nb_CNL_eval_ga','gamesslength', ...
    'output_ga.maxconstraint','tElapsed_ga','ok_ga','X_ga(1)', ...
    'X_ga(2)','X_ga(3)','X_ga(4)','X_ga(5)', ...
    'X_ga(6)','X_ga(7)','X_ga(8)','X_ga(9)','X_ga(10)','X_ga(11)','X_ga(12)',...
    'X_ga(13)','X_ga(14)','X_ga(15)','X_ga(16)','sol_ga', ...
    'exitflag_sqp','output_sqp.funccount','nb_CNL_eval_sqp', ...
    'sqpmesslength','output_sqp.constrviolation','tElapsed_sqp', ...
    'ok_sqp','X_sqp(1)','X_sqp(2)','X_sqp(3)', ...
    'X_sqp(4)','X_sqp(5)','X_sqp(6)','X_sqp(7)','X_sqp(8)','X_sqp(9)','X_sqp(10)',...
    'X_sqp(11)','X_sqp(12)','X_sqp(13)','X_sqp(14)','X_sqp(15)','X_sqp(16)', ...
    'sol_sqp','ec_geno_reduit','ec_rel_geno_reduit', ...
    'ec_d','ec_pheno','ec_rel_pheno','duree totale'};

ncol_res = length(entete);
xlswrite(fich_res, entete, 'bruts', 'b1');
xlswrite(fich_res, resu_br, 'bruts', 'b2');

txt_colA = {'moyenne';'minimum';'maximum';'';'% invalides';'duree totale'};
xlrange = ['a', num2str(Nexec + 3)];
xlswrite(fich_res, txt_colA, 'bruts', xlrange);

moyminmax = [mean(resu_br); min(resu_br); max(resu_br)];
xlrange = ['b', num2str(Nexec + 3)];
xlswrite(fich_res, moyminmax, 'bruts', xlrange);

pcent_invalides_ga = (Nexec - nb_valides_ga)*100/Nexec;
xlrange = ['b', num2str(Nexec + 7)];
xlswrite(fich_res, pcent_invalides_ga, 'bruts', xlrange);

xlrange = ['c', num2str(Nexec + 7)];
xlswrite(fich_res, '%', 'bruts', xlrange);

duree_totale = sum(resu_br(:,ncol_res))/3600; % en heures
xlrange = ['b', num2str(Nexec + 8)];
xlswrite(fich_res, duree_totale, 'bruts', xlrange);

hPC = ['h sur ',PC];
xlrange = ['c', num2str(Nexec + 8)];
xlswrite(fich_res, {hPC}, 'bruts', xlrange);

xlrange = ['a', num2str(Nexec + 10)];
xlswrite(fich_res, {aaaammjj}, 'bruts', xlrange);

% �criture du fichier r�sultat au format .xlsx : feuille 'valides ga tries vu_sqp'
xlswrite(fich_res, entete, 'valides ga tries vu_sqp', 'b1');

col_vu_sqp = 54;
resu_val_sorted = sortrows(resu_val,col_vu_sqp);
xlswrite(fich_res, resu_val_sorted, 'valides ga tries vu_sqp', 'b2');

xlrange = ['a', num2str(nb_valides_ga + 3)];
xlswrite(fich_res, txt_colA, 'valides ga tries vu_sqp', xlrange);

moyminmax = [mean(resu_val); min(resu_val); max(resu_val)];
xlrange = ['b', num2str(nb_valides_ga + 3)];
xlswrite(fich_res, moyminmax, 'valides ga tries vu_sqp', xlrange);

pcent_invalides_ga = 0;
xlrange = ['b', num2str(nb_valides_ga + 7)];
xlswrite(fich_res, pcent_invalides_ga, 'valides ga tries vu_sqp', xlrange);

xlrange = ['c', num2str(nb_valides_ga + 7)];
xlswrite(fich_res, '%', 'valides ga tries vu_sqp', xlrange);

duree_totale = sum(resu_val(:,ncol_res))/3600; % en heures
xlrange = ['b', num2str(nb_valides_ga + 8)];
xlswrite(fich_res, duree_totale, 'valides ga tries vu_sqp', xlrange);

xlrange = ['c', num2str(nb_valides_ga + 8)];
xlswrite(fich_res, {hPC}, 'valides ga tries vu_sqp', xlrange);

xlrange = ['a', num2str(nb_valides_ga + 10)];
xlswrite(fich_res, {aaaammjj}, 'valides ga tries vu_sqp', xlrange);

disp('calculation finished');