clear
%ax = [ -0.5, 1.5, -1.0, 1.0 ];
ax = [ -0.3, 1.3, -0.8, 0.5 ];
if (0)
	funch_map_theta = @(x,y)( 0.2*pi*( (x-1.0).^2  + y.^2 ) );
	%ax = [ -0.183, -0.181, 0.233, 0.236 ];
elseif (1)
	funch_map_theta = @(x,y)( 0.4*pi*( (x-1.0).^2  + y.^2 ) );
	%ax = [ 0.223, 0.226, -0.248, -0.245 ];
elseif (1)
	funch_map_theta = @(x,y)( 0.5*pi*( (x-1.0).^2  + y.^2 ) );
	%ax = [ 0.129, 0.131, -0.320, -0.317 ];
	%ax += [ -0.5, 0.5, -0.5, 0.5 ];
elseif (0)
	funch_map_theta = @(x,y)( pi*( (x-1.0).^2  + y.^2 ) );
else
	funch_map_theta = @(x,y)( 0*x );
end
funch_map_x = @(x,y)( ...
   (cos(funch_map_theta(x,y)).*x) ...
 - (sin(funch_map_theta(x,y)).*y) );
funch_map_y = @(x,y)( ...
   (sin(funch_map_theta(x,y)).*x) ...
 + (cos(funch_map_theta(x,y)).*y) );
%
%
funch_fx = @(x,y)( 0.5 + (x.*((x-1.0).^2)) );
funch_fy = @(x,y)( y );
funch_r0 = @(x,y)( funch_fx(x,y).^2 + funch_fy(x,y).^2 );
funch_r1 = @(x,y)( funch_r0(funch_map_x(x,y),funch_map_y(x,y)) );
%funch_z = @(x,y)( sqrt(funch_r1(x,y)) );
%funch_z = @(x,y)( (funch_r1(x,y)).^0.1 );
%funch_z = @(x,y)( funch_fy(x,y).^2 );
funch_z0 = @(x,y)( asinh(10000.0*funch_r1(x,y))/10000.0 );
%
funch_ry = @(x,y)( funch_fy(funch_map_x(x,y),funch_map_y(x,y)).^2 );
funch_z1 = @(x,y)( asinh(10000.0*funch_ry(x,y))/10000.0 );
%
%
numFigs = 0;
numContours = 30;
numColors = numContours+1;
sizeX = 51;
sizeY = 53;
multiArgLevel = 1;
%
numFigs++; figure(numFigs);
[ gridX, gridY, gridZ ] = gridfunch( funch_z0, multiArgLevel, ax, sizeX, sizeY );
contourf( gridX, gridY, gridZ.^4, numContours );
cMap = 0.6 + (0.4*jet(numColors));
z0 = ( 0.0 - min(min(gridZ)) ) / ( max(max(gridZ)) - min(min(gridZ)) );
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
%axis equal;
grid on;
return;
%
numContours = 30;
numColors = numContours+1;
numFigs++; figure(numFigs);
[ gridX, gridY, gridZ ] = gridfunch( funch_z1, multiArgLevel, ax, sizeX, sizeY );
contourf( gridX, gridY, gridZ.^2, numContours );
%imagesc( gridZ' );
cMap = 0.6 + (0.4*jet(numColors));
z0 = ( 0.0 - min(min(gridZ)) ) / ( max(max(gridZ)) - min(min(gridZ)) );
i0 = 1 + round( (numColors-1)*z0 );
if ( 0==i0 || -1==i0 ) % HA~ACK
	i0 = 1;
end
if ( 1 <= i0 && i0 <= numColors )
	%cMap(i0,:) *= 0.25;
	%cMap(i0,:) = 0.75 + 0.25*cMap(i0,:);
	cMap(i0,:) = 0.50 - 0.50*cMap(i0,:);
	cMap(i0,:) = 0.50 - 0.50*cMap(i0+1,:);
end
colormap(cMap);
%axis equal;
grid on;

return;
