clear
ax = [ -0.5, 1.5, -1.0, 1.0 ];
if (0)
	funch_map_theta = @(x,y)( 0.2*pi*( (x-1.0).^2  + y.^2 ) );
	%ax = [ -0.183, -0.181, 0.233, 0.236 ];
elseif (0)
	funch_map_theta = @(x,y)( 0.4*pi*( (x-1.0).^2  + y.^2 ) );
	ax = [ 0.223, 0.226, -0.248, -0.245 ];
elseif (1)
	funch_map_theta = @(x,y)( 0.5*pi*( (x-1.0).^2  + y.^2 ) );
	ax = [ 0.129, 0.131, -0.320, -0.317 ];
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
cMap = jet(256);
cMap(1,:) *= 0.2;
cMap(1,:) += 0.8;
cMap(end,:) *= 0.2;
%
numFigs++; figure(numFigs);
%[ gridX, gridY, gridZ ] = gridfunch( funch_z0, 1, ax, 51, 51 );
%contourf( gridX, gridY, gridZ,21 );
[ gridX, gridY, gridZ ] = gridfunch( funch_z0, 1, ax, 201, 201 );
contourf( gridX, gridY, gridZ, 51 );
colormap(cMap);
%axis equal;
grid on;
%
numFigs++; figure(numFigs);
%[ gridX, gridY, gridZ ] = gridfunch( funch_z1, 1, ax, 51, 51 );
%contourf( gridX, gridY, gridZ, 21 );
[ gridX, gridY, gridZ ] = gridfunch( funch_z1, 1, ax, 201, 201 );
contourf( gridX, gridY, gridZ, 51 );
colormap(cMap);
%axis equal;
grid on;

return;
