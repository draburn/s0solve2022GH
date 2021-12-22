clear;
tic();
numFigs = 0;
%
x0 = 0.0;
x1 = 2.0;
y0 = -1.0;
y1 = 1.0;
z0 = 0.0;
z1 = 1.0;
numXVals = 51;
numYVals = 53;
numZVals = 9;
xVals = linspace( x0, x1, numXVals );
yVals = linspace( y0, y1, numYVals );
zVals = linspace( z0, z1, numZVals );
[ aryX, aryY, aryZ ] = ndgrid( xVals, yVals, zVals );
%
aryF = (aryX.^2).*sin(2*pi*(aryY-0.5*aryZ));
fMax = max(max(max(aryF)));
fMin = min(min(min(aryF)));
%
sizeC = 256;
aryV = 1.0+((sizeC-1.0)*(aryF-fMin)/(fMax-fMin));
%
for z=1:numZVals
	numFigs++; figure(numFigs);
	colormap(jet(sizeC));
	%image( xVals, yVals, aryV(:,:,z)' );
	contourf( aryX(:,:,z), aryY(:,:,z), aryV(:,:,z), 31 );
	axis equal;
	title(sprintf( "z = %g", zVals(z) ));
end
