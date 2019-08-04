function [resu_best, pcent_invalides_ga, Nexec] = nvs20_naive_sd(tolfun_naive_ga, tolcon_naive_ga, Npop, seed_iter, I, gen_ext_ga_sqp, handle_CNL, handle_obj, tous_xi)

% The naive algorithm solves the nonconvex MINLP as follows:
% 0. Initialised parameter values
% 1. Initialise x0 in [Lx, Ux], y0 in [Ly, Uy] n N
% 2. Solve the integrality relaxed problem with GA using x0,y0 to yield x,y
% 3. Round y to the nearest integer to obtain y_tilda
% 4. Test if y in [Ly, Uy] n N with tolerance epsilon STOP else go to 5, epsilon = abs(int(y)-y)
% 5. Solve main problem with the integer variable fixed to obtain x_tilda using GA
% 6. Set x0 = x_tilda, y0 = y_tilda and go to step 1
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
options_ga.InitialPopulation = tous_xi;
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
    
    % tests if g has made changes in the structure options_ga
    ok_ga = isequal(initoptions_ga,options_ga); % = 1 if identical structures
    if ~ok_ga
       disp('***the options_ga structure has been changed by running ga ****')
       
    end
    
    gamesslength1 = length(regexprep(output_ga1.message,'\n',''));
    resu2sd = [exitflag_ga1, output_ga1.generations, output_ga1.funccount, ...
        nb_CNL_eval1,nb_obj_eval, gamesslength1, output_ga1.maxconstraint];
   
    resu3sd = y_ga;
      
    
% Step 3
% Round X(2) to the nearest integer 
% calculate the variable tol_text which is used to test tolerance    
    
    for i = 1: 5
        
        X(:,i) = int32(y_ga(:,i));
%     X(2) = int32(y_ga(2));
        tol_test(i) = abs(X(:,i) - y_ga(:,i));
%     x_sqp (:,1:5) = int32(x_sqp(:,1:5)); % convert X(2): p to integer
    end
    
% Step 4
% Error test
    
    if tol_test < tolfun_naive_ga
        
% Step 5
        fprintf('Solution found');
        
        resu1sd = [j, nb_grad_eval_sd, nb_obj_eval_sd, nb_CNL_eval_sd]; 
    
        tElapsed_sd = toc; % duration of seeding
        resu4sd = [0, 0, 0, 0, 0, 0, 0];
        resu5sd = x_ga;
        resuga2 = nvs20_obj(x_ga);
        resusd = [resu1sd, resu2sd, resu3sd, resu4sd, resu5sd, tElapsed_sd, resuga2];
        resu_br1 = [resu_br1; resusd];
        if exitflag_ga1 > 0
            resu_val1 = [resu_val1; resusd];
            nb_valides_ga = nb_valides_ga + 1;
        end
        resu_val_sorted = resu_val1;
        pcent_invalides_ga = (nb_valides_ga/seed_iter)*100; 
        resu_best = resu3sd;
        Nexec = j;
        return  

    else 
        % fix the integer variable X(2) and solve with GA to obtain y0
        tous_xi = y_ga;
        options_ga.InitialPopulation = tous_xi; % assign the new initial population to y_ga
         
       % launch of GA   

        raz_nb_CNL_eval1 = 1;
        raz_nb_obj_eval = 1;
        
        [x_ga, fval_ga2, exitflag_ga2, output_ga2] = ga(handle_obj2, I, ...
            [], [], [], [], lx, ux, handle_CNL, options_ga);
        
        % tests if g has made changes in the structure options_ga
        ok_ga = isequal(initoptions_ga,options_ga); % = 1 if identical structures
        if ~ok_ga
           disp('***the options_ga structure has been changed by running ga ****')
           
        end
    
        gamesslength2 = length(regexprep(output_ga2.message,'\n',''));
        resu4sd = [exitflag_ga2, output_ga2.generations, output_ga2.funccount, ...
            nb_CNL_eval1,nb_obj_eval, gamesslength2, output_ga2.maxconstraint];
   
        resu5sd = [x_ga, fval_ga2]; 
        %resuga2 = nvs20_obj(x_ga);
       
    end % end of if statement fot tolerance test
    
    resu1sd = [j, nb_grad_eval_sd, nb_obj_eval_sd, nb_CNL_eval_sd]; 
    
    tElapsed_sd = toc; % duration of seeding
    resusd = [resu1sd, resu2sd, resu3sd, resu4sd, resu5sd, tElapsed_sd];
    resu_br1 = [resu_br1; resusd];
    exitflag_ga = exitflag_ga1 + exitflag_ga2;
    if exitflag_ga > 0
        resu_val1 = [resu_val1; resusd];
        nb_valides_ga = nb_valides_ga + 1;
    end
    
% Step 6

tous_xi = x_ga; % sets the initial population to the new obtained solution for x_ga

j = j+1;

end % ends while statement

col_vu_sqp = 51;
pcent_invalides_ga = (nb_valides_ga/seed_iter)*100;

resu_val_sorted = sortrows(resu_val1,col_vu_sqp);
resu_best = [resu_val_sorted(1,35), resu_val_sorted(1,36), resu_val_sorted(1,37), resu_val_sorted(1,38), resu_val_sorted(1,39), resu_val_sorted(1,40), resu_val_sorted(1,41), resu_val_sorted(1,42), resu_val_sorted(1,43), resu_val_sorted(1,44), resu_val_sorted(1,45), resu_val_sorted(1,46), resu_val_sorted(1,47), resu_val_sorted(1,48), resu_val_sorted(1,49), resu_val_sorted(1,50)];
% fprintf('The best initial feasible solution from the naive algorithm is %d \n', resu_best);

disp('calculation finished');