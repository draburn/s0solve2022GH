if (1)
clear;
funch_fx = @(x)( x .* (x-1.0) .* (x+1.0) + 1.0 );
%%%funch_fy = @(y)( y );
funch_fy = @(y)( y );
funch_f = @(x,y)[ funch_fx(x); funch_fy(y) ];
eps_j = sqrt(eps);
funch_dfxdx = @(x)( (funch_fx(x+eps_j)-funch_fx(x-eps_j))/(2.0*eps_j) );
funch_dfydy = @(y)( (funch_fy(y+eps_j)-funch_fy(y-eps_j))/(2.0*eps_j) );
funch_j = @(x,y)( [ funch_dfxdx(x), 0.0; 0.0, funch_dfydy(y) ] );
funch_omega = @(x,y)( sum(funch_f(x,y).^2, 1) );
%
ax = [ -1.5, 1.5, -1.5, 1.5 ];
%ax = [ 0.57734, 0.57736, -0.00001, 0.00001 ];
[ gridX, gridY, gridZ ] = gridfunch( funch_omega, 1, ax, 201, 201 );
end
%
numFigs = 0;
cMap = jet(256);
cMap(1,:) *= 0.3;
cMap(1,:) += 0.7;
cMap(end,:) *= 0.3;
%
numFigs++; figure(numFigs);
contourf( gridX, gridY, asinh(100.0*gridZ)/100.0, 51 );
%%%imagesc( gridY', gridX', asinh(100.0*gridZ')/100.0 );
colormap(cMap);
axis equal;
grid on;
%
vecXBad = [ 0.577350269161; 0.0 ];
vecXGood = [ -1.324717957244746; 0.0 ];
%
%vecFBad = funch_f(vecXBad(1),vecXBad(2));
%funch_oldPRC = @(x,y)( vecFBad.'*funch_f(x,y) );
vecNu = [0;1];
funch_oldPRC = @(x,y)( vecNu.'*funch_f(x,y) );
numFigs++; figure(numFigs);
[ gridX, gridY, gridZ ] = gridfunch( funch_oldPRC, 1, ax, 201, 201 );
contourf( gridX, gridY, asinh(100.0*gridZ)/100.0, 51 );
colormap(cMap);
%
return;
%
eps_g = sqrt(eps);
funch_gx = @(x,y)(funch_omega(x+eps_g,y)-funch_omega(x-eps_g,y))/(2.0*eps_g);
funch_gy = @(x,y)(funch_omega(x,y+eps_g)-funch_omega(x,y-eps_g))/(2.0*eps_g);
funch_g = @(x,y)[ funch_gx(x,y); funch_gy(x,y) ];
%funch_g = @(x,y)[ ...
%  (funch_omega(x+eps_g,y)-funch_omega(x-eps_g,y))/(2.0*eps_g);
%  (funch_omega(x,y+eps_g)-funch_omega(x,y-eps_g))/(2.0*eps_g) ];
%
%vecX0 = [ 0.577; 0.2 ];
%%%vecX0 = [ 0.6; 0.5 ];
vecX0 = [ 0.5; 0.5 ];
%vecX0 = [ -0.2; 0.1 ];
epsD = 5e-3; epsD2 = 1e-3;
%epsD = 1e-4; epsD2 = 5e-4;
matX(:,1) = vecX0;
for n=1:100
	vecX = matX(:,n);
	%
	for m=1:10
	vecTX = -funch_g(vecX(1),vecX(2));
	matJ = funch_j(vecX(1),vecX(2));
	matS = diag(diag(matJ'*matJ));
	vecT0 = matS\vecTX;
	%! Use Marquardt-scaled gradient!
	%%%vecT0 = -funch_g(vecX(1),vecX(2));
	vecD = vecX-vecXBad;
	vecDHat = vecD / norm(vecD);
	vecT1 = vecT0 - (vecDHat'*vecT0)*vecDHat;
	%
	if (norm(vecT1)<0.1*epsD)
		epsD *= 1.1;
		if (epsD>1e-2)
			epsD = 1e-2;
		end
		break;
	end
	%
	vecDelta = epsD2 * (vecT1 / norm(vecT1) );
	vecX = vecX + vecDelta;
	%
	end
	%
	vecT0 = vecX-vecXBad;
	vecDelta = epsD * (vecT0 / norm(vecT0));
	vecX = vecX + vecDelta;
	%
	matX(:,n+1) = vecX;
end
%
%ax = [ -1.33, -1.32, -0.01, 0.01 ];
[ gridX, gridY, gridGX ] = gridfunch( funch_gx, 1, ax, 201, 201 );
numFigs++; figure(numFigs);
contourf( gridX, gridY, asinh(100.0*gridGX)/100.0, 51 );
%
[ gridX, gridY, gridGY ] = gridfunch( funch_gy, 1, ax, 201, 201 );
numFigs++; figure(numFigs);
contourf( gridX, gridY, asinh(100.0*gridGY)/100.0, 51 );
