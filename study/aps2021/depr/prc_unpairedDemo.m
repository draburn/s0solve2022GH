clear;
thisFile = "prc_unpairedDemo"
tic();
numFigs = 0;
%
%funch_fx = @(x,y)( x - 4.0./(1.0+(x-4.0).^2+(y-4.0).^2).^2 );
%funch_fy = @(x,y)( y + 4.0./(1.0+(x-4.0).^2+(y+4.0).^2).^2 );
funch_fx = @(x,y)( x - 4.0./(1.0+(x-4.0).^2+(y-2.8).^2).^2 );
funch_fy = @(x,y)( y - 4.0./(1.0+(x-2.8).^2+(y-4.0).^2).^2 );
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
%%%vecR0 = [ 3.8; 3.6 ];
%vecR0 = [ 2.5; 3.8 ];
vecR0 = [ 0.0; 0.0 ];
vecR = vecR0;
for n=1:10
	vecF = funch_f(vecR(1),vecR(2));
	omega = norm(vecF);
	msg( thisFile, __LINE__, sprintf( "%3d,  %g", n, omega ) );
	matJ = [ ...
	  ( funch_f(vecR(1)+epsX,vecR(2)) - funch_f(vecR(1)-epsX,vecR(2)) ) / (2.0*epsX), ...
	  ( funch_f(vecR(1),vecR(2)+epsY) - funch_f(vecR(1),vecR(2)-epsY) ) / (2.0*epsY) ];
	vecDelta = - matJ \ vecF;
	for lambda_trial = [ 1.0, 0.99, 0.9, 0.8, 0.5, 0.2, 0.1, ...
	  1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10  ]
		vecR_trial = vecR + lambda_trial*vecDelta;
		vecF_trial = funch_f(vecR_trial(1),vecR_trial(2));
		omega_trial = norm(vecF_trial);
		msg( thisFile, __LINE__, sprintf( "   %g,  %g", lambda_trial, omega_trial ) );
		if ( omega_trial < omega )
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
%ax = [ -5.0, 7.0, -6.0, 6.0 ];
%%%ax = [ 2.0, 4.5, 2.0, 4.5 ];
ax = [ 0.00, 0.01, 0.00, 0.01 ];
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
numColors = 31;
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
axis equal;
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
axis equal;
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
axis equal;
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
hold off;
grid on;
%
%
%
toc();
thisFile = [ "RETURN FROM " thisFile ];
return;
