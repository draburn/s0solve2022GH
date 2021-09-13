clear;
numFigs = 0;
%
use12Label = false;
funch_fx = @(x,y)( 0.5 - x .* (x-1.0) .* (x+1.0) );
funch_fy = @(x,y)( y );
funch_f = @(x,y)[ funch_fx(x,y); funch_fy(x,y) ];
funch_omega = @(x,y)( sum(funch_f(x,y).^2, 1) );
multiArgLevel_fx = 2;
multiArgLevel_fy = 2;
multiArgLevel_omega = 1;
ax = [ -1.4, 1.4, -1.4, 1.4 ];
%ax = [ -1.0, 1.3, -0.6, 0.6 ];
sizeX = 201;
sizeY = 203;
[ gridX, gridY, gridFX ] = gridfunch( funch_fx, multiArgLevel_fx, ax, sizeX, sizeY );
[ gridX, gridY, gridFY ] = gridfunch( funch_fy, multiArgLevel_fy, ax, sizeX, sizeY );
[ gridX, gridY, gridOmega ] = gridfunch( funch_omega, multiArgLevel_omega, ax, sizeX, sizeY );
%
numColors = 21;
numContours = 20;
assert( numColors <= numContours+1 );
%
numFigs++; figure(numFigs);
gridF = gridFX;
contourf( gridX, gridY, gridF, numContours );
cMap = 0.6 + (0.4*jet(numColors));
z0 = ( 0.0 - min(min(gridF)) ) / ( max(max(gridF)) - min(min(gridF)) );
i0 = 1 + round( (numColors-1)*z0 );
if ( 1 <= i0 && i0 <= numColors )
	%cMap(i0,:) *= 0.25;
	%cMap(i0,:) = 0.75 + 0.25*cMap(i0,:);
	cMap(i0,:) = 0.50 - 0.50*cMap(i0,:);
end
colormap(cMap);
if (use12Label)
	xlabel( "x_1" );
	ylabel( "x_2" );
	title( "F_1 vs x_1, x_2" );
else
	xlabel( "x" );
	ylabel( "y" );
	title( "F_x vs x, y" );
end
grid on;
%
numFigs++; figure(numFigs);
gridF = gridFY;
contourf( gridX, gridY, gridF, numContours );
cMap = 0.6 + (0.4*jet(numColors));
z0 = ( 0.0 - min(min(gridF)) ) / ( max(max(gridF)) - min(min(gridF)) );
i0 = 1 + round( (numColors-1)*z0 );
if ( 1 <= i0 && i0 <= numColors )
	%cMap(i0,:) *= 0.25;
	%cMap(i0,:) = 0.75 + 0.25*cMap(i0,:);
	cMap(i0,:) = 0.50 - 0.50*cMap(i0,:);
end
colormap(cMap);
if (use12Label)
	xlabel( "x_1" );
	ylabel( "x_2" );
	title( "F_2 vs x_1, x_2" );
else
	xlabel( "x" );
	ylabel( "y" );
	title( "F_y vs x, y" );
end
grid on;
%
numFigs++; figure(numFigs);
gridF = sqrt(gridOmega);
contourf( gridX, gridY, gridF, numContours );
cMap = 0.6 + (0.4*jet(numColors));
z0 = ( 0.0 - min(min(gridF)) ) / ( max(max(gridF)) - min(min(gridF)) );
i0 = 1 + round( (numColors-1)*z0 );
if ( 0==i0 || -1==i0 ) % HA~ACK
	i0 = 1;
end
if ( 1 <= i0 && i0 <= numColors )
	%cMap(i0,:) *= 0.25;
	%cMap(i0,:) = 0.75 + 0.25*cMap(i0,:);
	cMap(i0,:) = 0.50 - 0.50*cMap(i0,:);
end
colormap(cMap);
if (use12Label)
	xlabel( "x_1" );
	ylabel( "x_2" );
	title( "||F|| vs x_1, x_2" );
else
	xlabel( "x" );
	ylabel( "y" );
	title( "||F|| vs x, y" );
end
grid on;
