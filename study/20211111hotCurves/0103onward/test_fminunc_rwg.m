clear;
%
%fcn = @rwg
%fcn = @(x)rwg(x);
c = 10.0
fcn = @(x)rwgc(x,c);
vecX0 = [ 0.0; 0.0 ]
fcn0 = fcn(vecX0)
opts = optimset( 'GradObj', 'on' );
vecXf = fminunc( fcn, vecX0, opts )
