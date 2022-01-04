clear;
options = optimset('GradObj', 'on', 'MaxIter', 400 );
x0 = 1.0
fcn = @(x)( (x-2.0).^2 );
fcn(x0)
funfun(x0)
vecXF = fminunc( fcn, x0, options )
