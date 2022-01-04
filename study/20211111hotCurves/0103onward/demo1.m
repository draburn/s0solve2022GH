clear;
 ## Example for default optimization (Levenberg/Marquardt with
 ## BFGS), one non-linear equality constraint. Constrained optimum is
 ## at p = [0; 1].
 objective_function = @ (p) p(1)^2 + p(2)^2;
 
 switch 4
 case 1
 pin = [-2; 5];
 constraint_function = @ (p) p(1)^2 + 1 - p(2);
 [p, objf, cvg, outp] = nonlin_min (objective_function, pin, optimset ("equc", {constraint_function}))
 
 case 2
 pin = [-2; 5];
 [p, objf, cvg, outp] = nonlin_min (objective_function, pin)
 
 case 3
 constraint_function = @ (p) 4.0*p(1)^2 + p(2)^2 - 1.04;
 pin = [0.1;1.0]
 [p, objf, cvg, outp] = nonlin_min (objective_function, pin, optimset ("equc", {constraint_function}))
 
 case 4
 %constraint_function = @ (p) 4.0*p(1)^2 + p(2)^2 - 1.04;
 %constraint_function = @ (p) ( (p(1)-1)^2 + (p(2)-1)^2 - 1 );
 constraint_function = @ (p) -( (p(1)-1)^2 + (p(2)-1)^2 - 1.1 );
 pin = [0.0;1.0]
 [p, objf, cvg, outp] = nonlin_min (objective_function, pin, optimset ("inequc", {constraint_function}))
 end
 