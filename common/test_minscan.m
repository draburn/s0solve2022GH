clear;
commondefs;
thisFile = "test_minscan";
tic();
%
x0 = 1.0/sqrt(2.0);
funchOmega = @(x)( (x-x0).^2 );
prm.verbLev = VERBLEV__COPIOUS;
xBest = minscan( 0.0, 1.0, funchOmega, prm );
funchOmega(xBest)
toc();
