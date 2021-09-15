clear;
thisFile = "prc_soCalledOrig_zoomOnOther"
tic();
numFigs = 0;
%
if (0)
	funch_map_theta = @(x,y)( 0.4*pi*( (x-1.0).^2  + y.^2 ) );
	vecRAltExt = [0.663846; -0.747869];
	vecRRoot = [-0.0345211640940307; 0.2951445062707295];
	ax = [ -0.3, 1.3, -0.8, 0.5 ]; % Base zoom.
	%ax = [ 0.663, 0.665, -0.749, -0.747 ]; % Zoom on alt ext.
else
	funch_map_theta = @(x,y)( 0.5*pi*( (x-1.0).^2  + y.^2 ) );
	%vecRAltExt = [0.784404521330891; -0.620300195394128];
	vecRAltExt = [ 0.7898325, -0.6133225 ];
	vecRRoot = [0.0213073050733866; 0.2963916144271111];
	ax = [ -0.3, 1.3, -0.8, 0.5 ]; % Base zoom.
	%ax = [ 0.789, 0.791, -0.614, -0.613 ]; % Zoom on alt ext.
	%ax = [ 0.78982, 0.78984, -0.61333, -0.61331 ];
	%ax = [ 0.789832, 0.789833, -0.613323, -0.613322 ]; % Over-zoom on alt ext.
end
funch_map_x = @(x,y)( ...
   (cos(funch_map_theta(x,y)).*x) ...
 - (sin(funch_map_theta(x,y)).*y) );
funch_map_y = @(x,y)( ...
   (sin(funch_map_theta(x,y)).*x) ...
 + (cos(funch_map_theta(x,y)).*y) );
funch_gx = @(x,y)( 0.5 + (x.*((x-1.0).^2)) );
funch_gy = @(x,y)( y );
funch_fx = @(x,y)( funch_gx(funch_map_x(x,y),funch_map_y(x,y)) );
funch_fy = @(x,y)( funch_gy(funch_map_x(x,y),funch_map_y(x,y)) );
%
funch_f = @(x,y)[ funch_fx(x,y); funch_fy(x,y) ];
funch_omega = @(x,y)( sqrt(sum(funch_f(x,y).^2, 1)) );
%
%%%funch_soCalledR0 = @(x,y)( funch_gx(x,y).^2 + funch_gy(x,y).^2 );
%%%funch_soCalledR1 = @(x,y)( funch_soCalledR0(funch_map_x(x,y),funch_map_y(x,y)) );
%%%funch_soCalledZ0 = @(x,y)( asinh(10000.0*funch_soCalledR1(x,y))/10000.0 );
%
%
%
epsX = sqrt(eps);
epsY = sqrt(eps);
vecR0 = [ 1.0; 0.0 ];
vecR = vecR0;
for n=1:10
	vecF = funch_f(vecR(1),vecR(2));
	omega = norm(vecF);
	msg( thisFile, __LINE__, sprintf( "%3d,  %g", n, omega ) );
	matJ = [ ...
	  ( funch_f(vecR(1)+epsX,vecR(2)) - funch_f(vecR(1)-epsX,vecR(2)) ) / (2.0*epsX), ...
	  ( funch_f(vecR(1),vecR(2)+epsY) - funch_f(vecR(1),vecR(2)-epsY) ) / (2.0*epsY) ];
	%%%vecDelta = - matJ \ vecF;
	vecDelta = - (matJ'*matJ) \ (matJ'*vecF);
	for lambda_trial = [ 1.0, 0.99, 0.9, 0.8, 0.5, 0.2, 0.1, ...
	  1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, ...
	  1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16 ]
		vecR_trial = vecR + lambda_trial*vecDelta;
		vecF_trial = funch_f(vecR_trial(1),vecR_trial(2));
		omega_trial = norm(vecF_trial);
		msg( thisFile, __LINE__, sprintf( "   %g,  %g", lambda_trial, omega_trial ) );
		if ( omega_trial < omega )
			echo__lambda_trial = lambda_trial
			break;
		end
	end
	if ( omega_trial >= omega )
		break;
	end
	vecR = vecR_trial;
end
msg( thisFile, __LINE__, sprintf( "%3d,  %g", n, omega ) );
vecRExt = vecR
vecFExt = funch_f(vecRExt(1),vecRExt(2))
fExtNorm = sqrt( vecFExt(1)^2 + vecFExt(2)^2 );
if ( fExtNorm > 0.0 )
	vecFExtHat = vecFExt/fExtNorm;
else
	vecFExtHat = vecFExt;
end
matJExt = [ ...
 (  funch_f(vecRExt(1)+epsX,vecRExt(2)) ...
  - funch_f(vecRExt(1)-epsX,vecRExt(2)) ) / (2.0*epsX), ...
 (  funch_f(vecRExt(1),vecRExt(2)+epsY)
  - funch_f(vecRExt(1),vecRExt(2)-epsY) ) / (2.0*epsY) ]
%
%
%
use12Label = true;
multiArgLevel_fx = 2;
multiArgLevel_fy = 2;
multiArgLevel_omega = 1;
%%%
%ax = [ vecR(1)-0.01, vecR(1)+0.01, vecR(2)-0.01, vecR(2)+0.01 ];
%ax = [ -0.5, 1.5, -1.0, 1.0 ];
%ax = [ 0.4, 0.9, -0.9, -0.6 ];
%ax = [ 0.663846, 0.663847, -0.747869, -0.747868 ];
%%%
sizeX = 201;
sizeY = 203;
[ gridX, gridY, gridFX ] = gridfunch( funch_fx, multiArgLevel_fx, ax, sizeX, sizeY );
[ gridX, gridY, gridFY ] = gridfunch( funch_fy, multiArgLevel_fy, ax, sizeX, sizeY );
[ gridX, gridY, gridOmega ] = gridfunch( funch_omega, multiArgLevel_omega, ax, sizeX, sizeY );
valsFX = reshape(gridFX,1,[]);
valsFY = reshape(gridFY,1,[]);
valsOmega = reshape(gridOmega,1,[]);
%
valsJExtTF = matJExt' * [ valsFX; valsFY ];
valsOmegaA = sqrt( valsJExtTF(1,:).^2 + valsJExtTF(2,:).^2 );
gridOmegaA = reshape( valsOmegaA, sizeX, sizeY );
%
valsFExtHatTF = vecFExtHat' * [ valsFX; valsFY ];
gridFExtHatTF = reshape( valsFExtHatTF, sizeX, sizeY );
%
valsOut = ( eye(2,2) - vecFExtHat * vecFExtHat' ) * [ valsFX; valsFY ];
valsOmegaB = sqrt( valsOut(1,:).^2 + valsOut(2,:).^2 );
gridOmegaB = reshape( valsOmegaB, sizeX, sizeY );
%
%
%
numColors = 51;
numContours = 50;
assert( numColors <= numContours+1 );
%
numFigs++; figure(numFigs);
gridViz = gridOmega;
contourf( gridX, gridY, gridViz, numContours );
cMap = 0.6 + (0.4*jet(numColors));
z0 = ( 0.0 - min(min(gridViz)) ) / ( max(max(gridViz)) - min(min(gridViz)) );
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
hold on;
plot( vecR(1), vecR(2), "x", "color", [0.8,0.0,0.0], "linewidth", 3, "markersize", 20 );
plot( vecRAltExt(1), vecRAltExt(2), "+", "color", [0.0,1.0,0.0], "linewidth", 3, "markersize", 20 );
plot( vecRRoot(1), vecRRoot(2), "*", "color", [0.0,0.0,1.0], "linewidth", 3, "markersize", 20 );
hold off;
grid on;
%
%
%
numFigs++; figure(numFigs);
gridViz = gridFExtHatTF;
contourf( gridX, gridY, gridViz, numContours );
cMap = 0.6 + (0.4*jet(numColors));
z0 = ( 0.0 - min(min(gridViz)) ) / ( max(max(gridViz)) - min(min(gridViz)) );
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
	title( "F_e*F / ||F_e|| vs x_1, x_2" );
else
	xlabel( "x" );
	ylabel( "y" );
	title( "F_e*F / ||F_e|| vs x, y" );
end
hold on;
plot( vecR(1), vecR(2), "x", "color", [0.8,0.0,0.0], "linewidth", 3, "markersize", 20 );
plot( vecRAltExt(1), vecRAltExt(2), "+", "color", [0.0,1.0,0.0], "linewidth", 3, "markersize", 20 );
plot( vecRRoot(1), vecRRoot(2), "*", "color", [0.0,0.0,1.0], "linewidth", 3, "markersize", 20 );
hold off;
grid on;
%
%
numFigs++; figure(numFigs);
gridViz = gridOmegaB;
contourf( gridX, gridY, gridViz, numContours );
cMap = 0.6 + (0.4*jet(numColors));
z0 = ( 0.0 - min(min(gridViz)) ) / ( max(max(gridViz)) - min(min(gridViz)) );
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
	%title( "|| ( I - F_{ext}^{hat} F_{ext}^{hat T} ) * F || vs x_1, x_2" );
	title( "|| F - F_e (F_e*F) / (F_e*F_e) ||  vs x_1, x_2" );
else
	xlabel( "x" );
	ylabel( "y" );
	%title( "|| ( I - F_{ext}^{hat} F_{ext}^{hat T} ) * F ||  vs x, y" );
	title( "|| F - F_e (F_e*F) / (F_e*F_e) ||  vs x, y" );
end
hold on;
plot( vecR(1), vecR(2), "x", "color", [0.8,0.0,0.0], "linewidth", 3, "markersize", 20 );
plot( vecRAltExt(1), vecRAltExt(2), "+", "color", [0.0,1.0,0.0], "linewidth", 3, "markersize", 20 );
plot( vecRRoot(1), vecRRoot(2), "*", "color", [0.0,0.0,1.0], "linewidth", 3, "markersize", 20 );
hold off;
grid on;
%
%
%
toc();
thisFile = [ "RETURN FROM " thisFile ];
return;
