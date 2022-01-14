clear;
commondefs;
thisFile = "muSlipDemo";
setprngstates(0);
numFigs = 0;
tic();
%
penal = @(x,mu)( (x.^2)*mu/2.0 );
omega = @(x)(-( x+((x.^3)/3.0) ));
xi = @(x,mu)( omega(x) + penal(x,mu) );
xlo = 0.0;
xhi = 3.0;
x = linspace( xlo, xhi, 301 );
%
numFigs++; figure(numFigs);
plot( ...
  x, xi(x,4.0), 'o-', 'linewidth', 3, ...
  x, xi(x,3.0), 'o-', 'linewidth', 3, ...
  x, xi(x,2.0), 'o-', 'linewidth', 3, ...
  x, xi(x,1.0), 'o-', 'linewidth', 3 );
grid on;
xlabel( "x" );
ylabel( "xi" );
title( "muSlipDemo: xi vs x" );
axis( [ xlo, xhi, -xhi, xhi ] );
%
numFigs++; figure(numFigs);
plot( ...
  cent(x), diff(xi(x,4.0))./diff(x), 'o-', 'linewidth', 3, ...
  cent(x), diff(xi(x,3.0))./diff(x), 'o-', 'linewidth', 3, ...
  cent(x), diff(xi(x,2.0))./diff(x), 'o-', 'linewidth', 3, ...
  cent(x), diff(xi(x,1.0))./diff(x), 'o-', 'linewidth', 3 );
grid on;
xlabel( "x" );
ylabel( "dxi/dx" );
title( "muSlipDemo: dxidx vs x" );
axis( [ xlo, xhi, -xhi, xhi ] );
%
toc();
