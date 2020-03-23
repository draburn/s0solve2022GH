clear;
commondefs;
tic();
numFigs = 0;
%
funchZ = @(x)( x.^0.25 );
funchF = @(x)( [ funchZ(x).*cos(2*pi*funchZ(x)); funchZ(x).*sin(2*pi*funchZ(x)) ] );
%
x0 = 0.0;
x1 = 1.0;
numValsRequested = 30;
%
prm.funchFSupportsMultiArg = true;
rvecXDACLin = daclinspace( x0, x1, numValsRequested, funchF, prm );
rvecXLin = linspace( x0, x1, numValsRequested );
%
matFLin = funchF(rvecXLin);
matFDACLin = funchF(rvecXDACLin);
%
numFigs++; figure( numFigs );
plot( ...
  rvecXLin, matFLin(1,:), 'o-', ...
  rvecXDACLin, matFDACLin(1,:), 'x-' );
grid on;
%
numFigs++; figure( numFigs );
plot( ...
  rvecXLin, matFLin(2,:), 'o-', ...
  rvecXDACLin, matFDACLin(2,:), 'x-' );
grid on;
%
numFigs++; figure( numFigs );
plot( ...
  matFLin(1,:), matFLin(2,:), 'o-', ...
  matFDACLin(1,:), matFDACLin(2,:), 'o-' );
grid on;
