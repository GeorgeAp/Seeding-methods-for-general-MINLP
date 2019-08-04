function [resu_best, pcent_invalides_ga, Nexec] = nvs20_feasibility_sd(tolfun_naive_ga, tolcon_naive_ga, Npop, seed_iter, I, gen_ext_ga_sqp, handle_CNL, handle_obj_fixVal, tous_xi)

% The feasibility algorithm solves the general MINLP as follows:
% 0. Initialised parameter values
% 1. Generate x0 in [Lx, Ux], y0 in [Ly, Uy] n N at random
% 2. Run GA with fixed objective value 
% 3. Stop when number number of individuals > population number 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 0
tic; % lance chrono seeding
global MX X lx ux LX UX;

global resu_val_sorted;

handle_obj2 = @nvs20_obj_fixVal;
nb_valides_ga = 0;
resu_br1 = [];
resu_val1 = [];

global nb_CNL_eval1 raz_nb_CNL_eval1 nb_obj_eval raz_nb_obj_eval;


% Step 1
% randomly generate an initial population 

% generation of the pop.ini.random (or seeding) :
    for c = 1:Npop
        for i = 1:I
            tous_xi(c,i) = (LX(i) + (UX(i) - LX(i)) * rand)/MX(i);
        end
    end

nb_grad_eval_sd = 0; % no call to fct.calculation of gradients of g and h
nb_obj_eval_sd = 0; % no call to the fct.calcul of vu
nb_CNL_eval_sd = 0; % no call to the fct.calcul of CNL
% resu1 = [j, tElapsed_sd, nb_grad_eval_sd, nb_obj_eval_sd, nb_CNL_eval_sd];
    
% writing the "options_ga" structure
options_ga = [];
options_ga = gaoptimset(@ga);
options_ga.Display = 'off';
options_ga.MigrationFraction = 0.0; 
options_ga.PopulationSize = Npop;
options_ga.EliteCount = 2;
options_ga.CrossoverFraction = 0.8;
options_ga.Generations = gen_ext_ga_sqp;
options_ga.StallGenLimit = 50;
options_ga.TolFun = tolfun_naive_ga;
options_ga.TolCon = tolcon_naive_ga;
options_ga.InitialPenalty = 10;
options_ga.PenaltyFactor = 100;
options_ga.FitnessScalingFcn = @fitscalingrank; % function of scaling
options_ga.SelectionFcn = @selectionstochunif; % fct.de selection
options_ga.CrossoverFcn = @crossoverscattered; % fct.de crossing
options_ga.MutationFcn = @mutationadaptfeasible; % function of the mutation
initoptions_ga = options_ga; % memorization of the options_ga structure

% Step 2
% solve the integer relaxed problem with GA to obtain y0
j = 1;
while j <= seed_iter 
    
    options_ga.InitialPopulation = tous_xi;
    raz_nb_CNL_eval1 = 1;
    raz_nb_obj_eval = 1;
    
    % launch GA with relaxed integrality
    [y_ga, fval_ga, exitflag_ga1, output_ga1] = ga(handle_obj2, I, ...
        [], [], [], [], lx, ux, handle_CNL, options_ga);
    y_ga (:,1:5) = int32(y_ga(:,1:5));
    % tests if g has made changes in the structure options_ga
    ok_ga = isequal(initoptions_ga,options_ga); % = 1 if identical structures
    if ~ok_ga
       disp('***the options_ga structure has been changed by running ga ****')
       
    end
    
    gamesslength1 = length(regexprep(output_ga1.message,'\n',''));
    resu2sd = [exitflag_ga1, output_ga1.generations, output_ga1.funccount, ...
        nb_CNL_eval1, nb_obj_eval, gamesslength1, output_ga1.maxconstraint];
   
    resu3sd = [y_ga, fval_ga];
    
    resu1sd = [j, nb_grad_eval_sd, nb_obj_eval_sd, nb_CNL_eval_sd]; 
    
    tElapsed_sd = toc; % duration of seeding
    resusd = [resu1sd, resu2sd, resu3sd, tElapsed_sd];%, resu4sd, resu5sd, tElapsed_sd];
    resu_br1 = [resu_br1; resusd];
    exitflag_ga = exitflag_ga1;
    if exitflag_ga > 0
        resu_val1 = [resu_val1; resusd];
        nb_valides_ga = nb_valides_ga + 1;
    end
    
    % step 3
    tous_xi = y_ga;
    j = j+1;

end % ends while statement

% col_vu_sqp = 28;
col_vu_sqp = 29; %feasibility algo number of columns
pcent_invalides_ga = (nb_valides_ga/seed_iter)*100;
resu_val_sorted = sortrows(resu_val1,col_vu_sqp);
resu_best = [resu_val_sorted(1,11), resu_val_sorted(1,12), resu_val_sorted(1,13), resu_val_sorted(1,14), resu_val_sorted(1,15), resu_val_sorted(1,16), resu_val_sorted(1,17), resu_val_sorted(1,18), resu_val_sorted(1,19), resu_val_sorted(1,20), resu_val_sorted(1,21), resu_val_sorted(1,22), resu_val_sorted(1,23), resu_val_sorted(1,24), resu_val_sorted(1,25), resu_val_sorted(1,26)];
% fprintf('The best initial feasible solution from the naive algorithm is %d \n', resu_best);

disp('calculation finished');

