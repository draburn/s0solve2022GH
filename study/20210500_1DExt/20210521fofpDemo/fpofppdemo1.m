clear;
%setprngstates(79763216);
%setprngstates(57865504);
setprngstates(47337552);
%setprngstates();
numFigs = 0;

bigX = randn*exp(2*randn);
bigG0 = abs(randn);
cPL = abs(2+randn(3,1));
cGL = abs(randn(3,1));
cPR = abs(2+randn(3,1));
cGR = abs(randn(3,1));
funch_gL = @(x)( ...
   cGL(1)*abs(x-bigX).^cPL(1) ...
 + cGL(2)*abs(x-bigX).^cPL(2) ...
 + cGL(3)*abs(x-bigX).^cPL(3) );
funch_gR = @(x)( ...
   cGR(1)*abs(x-bigX).^cPR(1) ...
 + cGR(2)*abs(x-bigX).^cPR(2) ...
 + cGR(3)*abs(x-bigX).^cPR(3) );
funch_g = @(x)( bigG0 + (x<bigX).*funch_gL(x) + (x>bigX).*funch_gR(x) );
%funch_g = @(x)( exp(-(x-bigX).^2) )
%funch_g = @(x)( exp(abs(x-bigX)) - 1.0 )

x = linspace(bigX-0.01,bigX+0.01,1001);
%x = linspace(bigX-1.0,bigX+1.0,1001);

g = funch_g(x);
numFigs++; figure(numFigs);
plot( ...
  x, g, 'o-', ...
  bigX, bigG0, 'x', 'linewidth', 2, 'markersize', 20 );
grid on;

cg = cent(g);
dg = diff(g)./diff(x);
cx = cent(x);
ch1 = cg./dg;
%ch1 = cap( ch1, -5.0, 5.0 );
numFigs++; figure(numFigs);
plot( ...
  cx, ch1, 'o-', ...
  bigX, 0.0, 'x', 'linewidth', 2, 'markersize', 20 );
grid on;

ddg = diff(dg)./diff(cx);
ccx = cent(cx);
cdg = cent(dg);
cch2 = cdg./abs(ddg);
%cch2 = cap( cch2, -5.0, 5.0 );
numFigs++; figure(numFigs);
plot( ...
  ccx, cch2, 'o-', ...
  bigX, 0.0, 'x', 'linewidth', 2, 'markersize', 20 );
grid on;
