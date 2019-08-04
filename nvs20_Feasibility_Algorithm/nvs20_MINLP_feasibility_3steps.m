% Algorithm for nonvonvex mixed integer nonlinear programming problem
% the algorithm is divided into 3 parts
% A. Generate initial feasible solution using naive algorithm
% B. Solve the problem using GA
% C. Refine solution using SQP
% 
% 
% This script calls the naive_algo_imple function which solves the Step A
% The Steps B and C are then solved afterwards
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
fich_res = ['res_feasibility_', problem, '_s3cx_pop', num2str(Npop), ...
    '_gen', num2str(gen_ext_ga_sqp),'_stall', num2str(stall_gen_limit_ga), ...
    '_', aaaammjj, '.xlsx'];

list = ls(fich_res);
if ~isempty(list)
    disp('result file: already exists: can not continue');
    return
end

I = 16; % number of variables

% terminals LX, UX, and reduced terminals lx, ux, on input var. :
global MX lx ux LX UX;
LX = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; 
UX = [200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200];
MX = (LX + UX)*0.5; % average of the terminals

lx = LX./MX;
ux = UX./MX;


% bounds of nonlinear inequality constraints:
global LG UG MG;

LG = [2.5 1.1 -3.1 -3.5 1.3 2.1 2.3 -1.5];
UG = [100 100 100 100 100 100 100 100];
MG = (LG + UG)*0.5; % average of the bounds

handle_CNL = @nvs20_CNL; % constraint function
% handle_CNL_fixInt6 = @mapse_pvar_CNL_fixInt; % constraint function for fixed integer
handle_obj = @nvs20_obj; % objective function
% handle_obj_fixInt6 = @mapse_pvar_vu_fixInt; % objective function with fixed p
tous_xi = zeros(Npop,I); % Npop initial individuals (1 per line)

global resu_br1 nb_valides_ga resu_val_sorted;
nb_valides_ga = 0;
resu_br = [];
resu_br1 = [];
resu_val = []; 

global nb_CNL_eval raz_nb_CNL_eval;

%STEP A
% Initial feasible solution
% for loop on the Nexec executions% generation of the pop.ini.random (or seeding) :
seed_iter = 100;
Nexec = 100;
for iexec = 1:Nexec
    
    % input for function 'naive_algo_imple' is number of iterations
    [tous_yi_3step, pcent_invalides_ga] = nvs20_feasibility_sd(tolfun_naive_ga, tolcon_naive_ga, Npop, seed_iter, I, gen_ext_ga_sqp, handle_CNL, handle_obj); 
    new_tous_xi = tous_yi_3step;
    
    nb_grad_eval_sd = 0; % no call to fct.calculation of gradients of g and h
    resu_seeding = resu_val_sorted;
%     tElapsed_sd = sum(resu_seeding(1,52));
    tElapsed_sd = sum(resu_seeding(1,29));% for feasibility
%     nb_obj_eval_sd = (sum(resu_seeding(1,8)) + sum(resu_seeding(1,31)));
%     nb_CNL_eval_sd = (sum(resu_seeding(1,7)) + sum(resu_seeding(1,30)));
    nb_obj_eval_sd = sum(resu_seeding(1,8));% for feasibility
    nb_CNL_eval_sd = sum(resu_seeding(1,7));% for feasibility
    resu1 = [iexec, tElapsed_sd, nb_grad_eval_sd, nb_obj_eval_sd, nb_CNL_eval_sd, pcent_invalides_ga]; 
%     end
%     p = X(2);
    
% STEP B
% Augmented Lagrangian Genetic Algorithm (ALGA)

    % writing the structure of "options_ga"
    options_ga3step = [];
    options_ga3step = gaoptimset(@ga);
    options_ga3step.InitialPopulation = tous_yi_3step;
    options_ga3step.Display = 'off';
    options_ga3step.MigrationFraction = 0.0; % no migration
    options_ga3step.PopulationSize = Npop;
    options_ga3step.EliteCount = 2;
    options_ga3step.CrossoverFraction = 0.8;
    options_ga3step.Generations = gen_ext_ga_sqp;
    options_ga3step.StallGenLimit = stall_gen_limit_ga;
    options_ga3step.TolFun = tolfun_ga;
    options_ga3step.TolCon = tolcon_ga;
    options_ga3step.InitialPenalty = 10;
    options_ga3step.PenaltyFactor = 100;    
    options_ga3step.FitnessScalingFcn = @fitscalingrank; % function for scaling
    options_ga3step.SelectionFcn = @selectionstochunif; % function for selection
    options_ga3step.CrossoverFcn = @crossoverscattered; % function for crossing
    options_ga3step.MutationFcn = @mutationadaptfeasible; % croisement mutation
    initoptions_ga3step = options_ga3step; % memorisation of the structure of options.ga
    
    % launch of GA :
    raz_nb_CNL_eval = 1;
    tic; % lauch of timer for GA
    [x_ga, fval_ga, exitflag_ga, output_ga] = ga(handle_obj, I, ...
        [], [], [], [], lx, ux, handle_CNL, options_ga3step);
    tElapsed_ga = toc; % duration of GA
    x_ga (:,1:5) = int32(x_ga(:,1:5));
    
    % tests if ga has made changes in the structure options_ga
    ok_ga = isequal(initoptions_ga3step,options_ga3step); % = 1 if identical structures
    if ~ok_ga
       disp('***the options_ga structure has been changed by running ga ****')
    end

    gamesslength = length(regexprep(output_ga.message,'\n',''));
    resu2 = [exitflag_ga, output_ga.generations, output_ga.funccount, ...
        nb_CNL_eval, gamesslength, output_ga.maxconstraint, ...
        tElapsed_ga, ok_ga];

    % determination of complete solution :
    [xy_sol] = handle_obj(x_ga);
%      vu_ga = xy_sol; va_ga = 1; pj_ga = 1; sigmaem_ga = 1; ech_ga = 1; kf_ga = 1; be_ga = 1; d_ga = 1; c_ga = 1;
    vu_ga = xy_sol;
%     X_ga = x_ga.*MX; % conversion var.reduced --> var.non reduced
    X_ga = x_ga.*MX;
    resu3 = [X_ga, vu_ga];
    
% STEP C
% Refinement
% Sequential Quadratic Programming

    % algorithm of local descent (post-treatment of the sol.ga)
    options_sqp = [];
    options_sqp = optimoptions('fmincon');
    options_sqp.Algorithm = 'sqp';
    options_sqp.Display = 'off';
    options_sqp.MaxFunEvals = 100000;
    options_sqp.MaxIter = 100000;
    options_sqp.TolCon = tolcon_ga;
    options_sqp.TolFun = tolfun_ga;
    options_sqp.TolX =   1.0e-6;
    options_sqp.TypicalX = ones(1,I);
    initoptions_sqp = options_sqp; % memorisation of the structure options_sqp
    

    % Round X(2) to the nearest integer 
    % calculate the variable tol_text which is used to test tolerance
    
    % p = X(2);
%     X(2) = int32(x_ga(2)); 

    % launch of sqp :
 
    raz_nb_CNL_eval = 1;
    tic; % launch chrono sqp
    [x_sqp, fval_sqp, exitflag_sqp, output_sqp] = fmincon(handle_obj, ...
        x_ga, [], [], [], [], lx, ux, handle_CNL, options_sqp);
    tElapsed_sqp = toc; % duration of sqp
    x_sqp(:,1:5) = x_ga(:,1:5); % fixed integers

    % tests if fmincon made changes in the struct.options_sqp
    ok_sqp = isequal(initoptions_sqp,options_sqp); % = 1 if struct.identiques
    if ~ok_sqp
       warning('the options_sqp structure has been changed by running fmincon')
    end

    sqpmesslength = length(regexprep(output_sqp.message,'\n',''));
    resu4 = [ exitflag_sqp, ...
        output_sqp.funcCount, nb_CNL_eval, sqpmesslength, ...
        output_sqp.constrviolation, tElapsed_sqp, ok_sqp];
    
    % complete determination of the solution :
%     [vu_sqp, va_sqp, pj_sqp, sigmaem_sqp, ech_sqp, kf_sqp, be_sqp, d_sqp, c_sqp] = mapse_pvar_3val(x_sqp);
    [xy_sol] = handle_obj(x_sqp);
%     vu_sqp = xy_sol; va_sqp = 1; pj_sqp = 1; sigmaem_sqp = 1; ech_sqp = 1; kf_sqp = 1; be_sqp = 1; d_sqp = 1; c_sqp = 1;
    vu_sqp = xy_sol; 

%     X_sqp = x_sqp.*MX; % conversion var.reduced --> var.not reduced
    X_sqp = x_sqp.*MX;
    resu5 = [ X_sqp, vu_sqp];

    
    % evaluation of differences between solution of ga and sqp :
    %     (formules vérifiées 20180404)
    ec_geno = norm(X_sqp - X_ga); % absolute difference in X: genotype
    ec_rel_geno = ec_geno/norm(X_sqp); % relative difference in X: genotype
    ec_geno_reduit = norm(x_sqp - x_ga); % absolute deviation in x (var.reduced): genotype
    ec_rel_geno_reduit = ec_geno_reduit/norm(x_sqp); % relative difference in x (var.reduced): genotype
    ec_d = 0;% Absolute deviation in p (shows that sqp only rounds p)
    ec_pheno = vu_sqp - vu_ga; % absolute difference in view: phenotype
    ec_rel_pheno = ec_pheno/vu_sqp; % relative difference in view: phenotype
    
    duree_totale = tElapsed_sd + tElapsed_ga + tElapsed_sqp;
    resu6 = [ec_geno_reduit, ec_rel_geno_reduit, ec_d, ec_pheno, ec_rel_pheno, duree_totale];
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
    nx = (x_sqp + ux)*0.5;

end

warning off MATLAB:xlswrite:AddSheet

% écriture du fichier résultat au format .xlsx : feuille 'param'
noms_params = {'paramètres ga :'; 'TolFun'; 'TolCon'; 'population'; ...
    'générations'; 'stall gen.';''; 'paramètres sqp :'; 'TolCon'; ...
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

% écriture du fichier résultat au format .xlsx : feuille 'bruts'
entete = {'iexec','tElapsed_sd','nb_grad_eval_sd','nb_obj_eval_sd', ...
    'nb_CNL_eval_sd','percent_valid_sd','exitflag_ga','output_ga.generations', ...
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

if nb_valides_ga > 0
    % écriture du fichier résultat au format .xlsx : feuille 'valides ga tries vu_sqp'
    xlswrite(fich_res, entete, 'valides ga tries vu_sqp', 'b1');

    col_vu_sqp = 55;
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
end

disp('calculation finished');