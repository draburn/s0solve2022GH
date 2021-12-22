clear;
thisFile = "func1111_init2";
commondefs;
tic();
setprngstates(6352560); % "Good behavior" but not single ext.
%setprngstates();
numFigs = 0;
%
%msg( thisFile, __LINE__, "This simple case illustrates some issues with existing handling of 'bad points'." );
%
sizeX = 2;
sizeF = 2;
isBadMin = true;
forceSingleExt = false;
forceGoodBehavior = true; % Likely single ext, but not always.
%
x1Lo = -5.0; x1Hi = +5.0; x2Lo = -5.0; x2Hi = +5.0;
%x1Lo = -2.4; x1Hi = -1.8; x2Lo = -0.2; x2Hi = 0.4;
%
%%% THIS STUFF ISN'T USED!!!
funcPrm.ary3K = randn(sizeX,sizeX,sizeF);
matJPre = randn(sizeF,sizeX);
vecFE = randn(sizeF,1);
funcPrm.vecXE = randn(sizeX,1);
%
vecFEHat = vecFE/norm(vecFE);
if (forceSingleExt)
	funcPrm.ary3K(:,:,:) = 0.0;
	for n=1:sizeX
		funcPrm.ary3K(n,n,:) = abs(randn)*vecFEHat;
	end
elseif (forceGoodBehavior)
	funcPrm.ary3K(:,:,:) = 0.0;
	for n=1:sizeX
	for m=1:sizeX
		vecT1 = randn(sizeF,1);
		vecT2 = vecT1 - vecFEHat*(vecFEHat'*vecT1);
		if (n==m)
			funcPrm.ary3K(n,m,:) = vecT2 + abs(randn)*vecFEHat;
		else
			funcPrm.ary3K(n,m,:) = vecT2;
		end
	end
	end
end
%
if (isBadMin)
	funcPrm.vecFE = vecFE;
	funcPrm.matJ = matJPre - vecFE*(vecFE'*matJPre)/(vecFE'*vecFE);
else
	funcPrm.vecFE = zeros(sizeF,1);
	funcPrm.matJ = matJPre;
end
%%%funcPrm.matJ *= 0.0;
%
funcPrm.sizeX = sizeX;
funcPrm.sizeF = sizeF;
for n=1:sizeF
	matTemp = funcPrm.ary3K(:,:,n);
	funcPrm.ary3K(:,:,n) = ( matTemp' + matTemp ) / 2.0;
end
%
%
%
matI = eye(2,2);
numSVals = 100;
sVals = linspace(0.0,1.0,numSVals);
%
vecXLev0 = [ 4.0; 4.0 ];
[ vecFDemo, matJDemo, vecGDemo, matHDemo ] = func1111_evalDeriv( vecXLev0, funcPrm );
%
vecGLev0 = vecGDemo;
matHLev0 = matHDemo;
vecXLevJTJ0 = vecXLev0;
vecFLevJTJ0 = vecFDemo;
matJLevJTJ0 = matJDemo;
vecGLevJTJ0 = matJLevJTJ0'*vecFLevJTJ0; % Should be same as vecGLev0 if at same pt.
matHLevJTJ0 = matJLevJTJ0'*matJLevJTJ0;
for n=1:numSVals
	s = sVals(n);
	matDeltaLev(:,n) = -( s*matHLev0 + (1.0-s)*matI ) \ ( s*vecGLev0 );
	matXLev(:,n) = vecXLev0 + matDeltaLev(:,n);
	matDeltaLevJTJ(:,n) = -( s*matHLevJTJ0 + (1.0-s)*matI ) \ ( s*vecGLevJTJ0 );
	matXLevJTJ(:,n) = vecXLevJTJ0 + matDeltaLevJTJ(:,n);
end
%
vecXLevB0 = matXLev(:,end);
[ vecFDemo, matJDemo, vecGDemo, matHDemo ] = func1111_evalDeriv( vecXLevB0, funcPrm );
vecGLevB0 = vecGDemo;
matHLevB0 = matHDemo;
vecXLevJTJB0 = matXLevJTJ(:,end);
[ vecFDemo, matJDemo, vecGDemo, matHDemo ] = func1111_evalDeriv( vecXLevJTJB0, funcPrm );
vecFLevJTJB0 = vecFDemo;
matJLevJTJB0 = matJDemo;
vecGLevJTJB0 = matJLevJTJB0'*vecFLevJTJB0; % Should be same as vecGLev0 if at same pt.
matHLevJTJB0 = matJLevJTJB0'*matJLevJTJB0;
for n=1:numSVals
	s = sVals(n);
	matDeltaLevB(:,n) = -( s*matHLevB0 + (1.0-s)*matI ) \ ( s*vecGLevB0 );
	matXLevB(:,n) = vecXLevB0 + matDeltaLevB(:,n);
	matDeltaLevJTJB(:,n) = -( s*matHLevJTJB0 + (1.0-s)*matI ) \ ( s*vecGLevJTJB0 );
	matXLevJTJB(:,n) = vecXLevJTJB0 + matDeltaLevJTJB(:,n);
end
%
vecXLevC0 = matXLevB(:,end);
[ vecFDemo, matJDemo, vecGDemo, matHDemo ] = func1111_evalDeriv( vecXLevC0, funcPrm );
vecGLevC0 = vecGDemo;
matHLevC0 = matHDemo;
vecXLevJTJC0 = matXLevJTJB(:,end);
[ vecFDemo, matJDemo, vecGDemo, matHDemo ] = func1111_evalDeriv( vecXLevJTJC0, funcPrm );
vecFLevJTJC0 = vecFDemo;
matJLevJTJC0 = matJDemo;
vecGLevJTJC0 = matJLevJTJC0'*vecFLevJTJC0; % Should be same as vecGLev0 if at same pt.
matHLevJTJC0 = matJLevJTJC0'*matJLevJTJC0;
for n=1:numSVals
	s = sVals(n);
	matDeltaLevC(:,n) = -( s*matHLevC0 + (1.0-s)*matI ) \ ( s*vecGLevC0 );
	matXLevC(:,n) = vecXLevC0 + matDeltaLevC(:,n);
	matDeltaLevJTJC(:,n) = -( s*matHLevJTJC0 + (1.0-s)*matI ) \ ( s*vecGLevJTJC0 );
	matXLevJTJC(:,n) = vecXLevJTJC0 + matDeltaLevJTJC(:,n);
end
%
numGradSteps = 1000;
vecXHOTGrad0 = vecXLev0;
vecXHOTGrad = vecXHOTGrad0;
matXHOTGrad(:,1) = vecXHOTGrad;
if (0)
%% THIS CAN DIVERGE!
for n=1:numGradSteps
	[ vecFTemp, matJTemp ] = func1111_evalDeriv( vecXHOTGrad, funcPrm );
	vecGTemp = matJTemp'*vecFTemp;
	vecXHOTGrad -= 0.01*vecGTemp;
	matXHOTGrad(:,n+1) = vecXHOTGrad;
end
end
%
%
maxStepSize_HOTGC = 0.1;
minStepSize_HOTGC = 0.01;
vecX_HOTGC = vecXLev0;
[ vecF_HOTGC, matJ_HOTGC, vecG_HOTGC, matH_HOTGC ] = func1111_evalDeriv( vecX_HOTGC, funcPrm );
assert( 0.0 < norm(vecG_HOTGC) );
omega = (vecF_HOTGC'*vecF_HOTGC)/2.0;
matX_HOTGC(:,1) = vecX_HOTGC;
%
vecG0_HOTGC = vecG_HOTGC;
omega0 = omega;
%
vecDelta = -minStepSize_HOTGC*vecG_HOTGC/norm(vecG_HOTGC); % Unless...
[ matR_HOTGC, cholFlag ] = chol( matH_HOTGC );
if ( 0==cholFlag )
	vecDeltaN = -( matR_HOTGC\(matR_HOTGC'\vecG_HOTGC) );
	if (norm(vecDeltaN) <= maxStepSize_HOTGC )
		vecDelta = vecDeltaN;
	end
end
omegaPrev = omega;
%
n = 1;
while (1)
	%
	% Try new step.
	vecX_HOTGC += vecDelta;
	[ vecF_HOTGC, matJ_HOTGC, vecG_HOTGC, matH_HOTGC ] = func1111_evalDeriv( vecX_HOTGC, funcPrm );
	omega = (vecF_HOTGC'*vecF_HOTGC)/2.0;
	if ( omega >= omegaPrev )
		break;
	end
	n++;
	matX_HOTGC(:,n) = vecX_HOTGC;
	%
	if ( norm(vecG_HOTGC) <= eps050*norm(vecG0_HOTGC) )
		break;
	end
	if ( n >= 10000 )
		break;
	end
	%
	% Calc next step.
	vecDelta = -minStepSize_HOTGC*vecG_HOTGC/norm(vecG_HOTGC); % Unless...
	[ matR_HOTGC, cholFlag ] = chol( matH_HOTGC );
	if ( 0==cholFlag )
		vecDeltaN = -( matR_HOTGC\(matR_HOTGC'\vecG_HOTGC) );
		if (norm(vecDeltaN) <= maxStepSize_HOTGC )
			vecDelta = vecDeltaN;
		end
	end
	omegaPrev = omega;
end
msg( thisFile, __LINE__, sprintf( "HOTGC: %d steps to %10.3e.", n, omega ) );
%
%
%
% NOT ACTUALLY RK4!
stepSize = 0.03;
%
vecX0_RK4HOTGC = vecXLev0;
[ vecF0_RK4HOTGC, matJ0_RK4HOTGC ] = func1111_evalDeriv( vecX0_RK4HOTGC, funcPrm );
vecG0_RK4HOTGC = matJ0_RK4HOTGC'*vecF0_RK4HOTGC;
assert( 0.0 < norm(vecG0_RK4HOTGC) );
%
vecX_RK4HOTGC = vecX0_RK4HOTGC;
matX_RK4HOTGC(:,1) = vecX_RK4HOTGC;
omegaPrev = 0.5*vecF0_RK4HOTGC'*vecF0_RK4HOTGC;
n = 1;
while (true)
	n++;
	if ( n>=1000 )
		break;
	end
	[ vecFTemp, matJTemp ] = func1111_evalDeriv( vecX_RK4HOTGC, funcPrm );
	vecGTemp = matJTemp'*vecFTemp;
	omega = 0.5*vecFTemp'*vecFTemp;
	%msg( thisFile, __LINE__, sprintf( ...
	%  "RK4HOTGC: %5d,  %10.3e,  %10.3e,  %10.3e.", ...
	%  n, vecX_RK4HOTGC(1), vecX_RK4HOTGC(2), omega ) );
	if ( omega > omegaPrev )
		break;
	end
	if ( norm(vecGTemp) <= eps050*norm(vecG0_RK4HOTGC) )
		break;
	end
	vecX_RK4HOTGC -= stepSize*vecGTemp/norm(vecGTemp);
	matX_RK4HOTGC(:,n) = vecX_RK4HOTGC;
	omegaPrev = omega;
end
%%%matX_RK4HOTGC = matX_RK4HOTGC(:,1:end-1);
msg( thisFile, __LINE__, sprintf( "FAKE RK4HOTGC: %d steps", n ) );
%
%
%
numX1Vals = 51;
numX2Vals = 51;
%
x1Vals = linspace(x1Lo,x1Hi,numX1Vals);
x2Vals = linspace(x2Lo,x2Hi,numX2Vals);
[ x1Mesh, x2Mesh ] = meshgrid( x1Vals, x2Vals );
matX = [ reshape(x1Mesh,1,[]); reshape(x2Mesh,1,[]) ];
matF = func1111_eval(matX,funcPrm);
f1Mesh = reshape(matF(1,:),numX2Vals,numX1Vals);
f2Mesh = reshape(matF(2,:),numX2Vals,numX1Vals);
omegaMesh = 0.5*( f1Mesh.^2 + f2Mesh.^2 );
%
msg( thisFile, __LINE__, sprintf( "omega scale: %g to %g.", ...
  min(min((omegaMesh))), max(max((omegaMesh))) ) );
msg( thisFile, __LINE__, sprintf( "log(omega) scale: %g to %g.", ...
  min(min(log(omegaMesh))), max(max(log(omegaMesh))) ) );
msg( thisFile, __LINE__, sprintf( "sqrt(omega) scale: %g to %g.", ...
  min(min(sqrt(omegaMesh))), max(max(sqrt(omegaMesh))) ) );
msg( thisFile, __LINE__, sprintf( "|F1| scale: %g to %g.", ...
  min(min(abs(f1Mesh))), max(max(abs(f1Mesh))) ) );
msg( thisFile, __LINE__, sprintf( "|F2| scale: %g to %g.", ...
  min(min(abs(f2Mesh))), max(max(abs(f2Mesh))) ) );
msg( thisFile, __LINE__, sprintf( "F1 scale: %g to %g.", ...
  min(min((f1Mesh))), max(max((f1Mesh))) ) );
msg( thisFile, __LINE__, sprintf( "F2 scale: %g to %g.", ...
  min(min((f2Mesh))), max(max((f2Mesh))) ) );

%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, log(eps050*max(max(omegaMesh))+omegaMesh-min(min(omegaMesh))), 30 );
hold on;
plot( funcPrm.vecXE(1), funcPrm.vecXE(2), 'w+', 'markersize', 20 );
plot( funcPrm.vecXE(1), funcPrm.vecXE(2), 'ko', 'markersize', 20 );
plot( matXLevJTJ(1,:), matXLevJTJ(2,:), 'k-', 'linewidth', 2 );
plot( matXLevJTJB(1,:), matXLevJTJB(2,:), 'c-', 'linewidth', 2 );
plot( matXLevJTJC(1,:), matXLevJTJC(2,:), 'm-', 'linewidth', 2 );
plot( matXLev(1,:), matXLev(2,:), 'ro-' );
plot( matXLevB(1,:), matXLevB(2,:), 'gs-' );
plot( matXLevC(1,:), matXLevC(2,:), 'b^-' );
%plot( matXHOTGrad(1,:), matXHOTGrad(2,:), 'y*-' );
%plot( matX_RK4HOTGC(1,:), matX_RK4HOTGC(2,:), 'k*-', 'linewidth', 2 );
plot( matX_HOTGC(1,:), matX_HOTGC(2,:), 'y*-', 'linewidth', 2 );
grid on;
hold off;
axis equal;
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "log(omega-omegaMin) vs x1, x2" );
return
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, sqrt(omegaMesh), 30 );
hold on;
plot( funcPrm.vecXE(1), funcPrm.vecXE(2), 'w+', 'markersize', 20 );
plot( funcPrm.vecXE(1), funcPrm.vecXE(2), 'ko', 'markersize', 20 );
plot( matXLevJTJ(1,:), matXLevJTJ(2,:), 'k-', 'linewidth', 2 );
plot( matXLevJTJB(1,:), matXLevJTJB(2,:), 'c-', 'linewidth', 2 );
plot( matXLevJTJC(1,:), matXLevJTJC(2,:), 'm-', 'linewidth', 2 );
plot( matXLev(1,:), matXLev(2,:), 'ro-' );
plot( matXLevB(1,:), matXLevB(2,:), 'gs-' );
plot( matXLevC(1,:), matXLevC(2,:), 'b^-' );
%plot( matXHOTGrad(1,:), matXHOTGrad(2,:), 'y*-' );
%plot( matX_RK4HOTGC(1,:), matX_RK4HOTGC(2,:), 'k*-', 'linewidth', 2 );
plot( matX_HOTGC(1,:), matX_HOTGC(2,:), 'y*-', 'linewidth', 2 );
grid on;
hold off;
axis equal;
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "sqrt(omega) vs x1, x2" );
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, abs(f1Mesh), 30 );
hold on;
plot( funcPrm.vecXE(1), funcPrm.vecXE(2), 'w+', 'markersize', 20 );
plot( funcPrm.vecXE(1), funcPrm.vecXE(2), 'ko', 'markersize', 20 );
plot( matXLevJTJ(1,:), matXLevJTJ(2,:), 'k-', 'linewidth', 2 );
plot( matXLevJTJB(1,:), matXLevJTJB(2,:), 'c-', 'linewidth', 2 );
plot( matXLevJTJC(1,:), matXLevJTJC(2,:), 'm-', 'linewidth', 2 );
plot( matXLev(1,:), matXLev(2,:), 'ro-' );
plot( matXLevB(1,:), matXLevB(2,:), 'gs-' );
plot( matXLevC(1,:), matXLevC(2,:), 'b^-' );
%%%plot( matXHOTGrad(1,:), matXHOTGrad(2,:), 'y*-' );
plot( matX_HOTGC(1,:), matX_HOTGC(2,:), 'y*-', 'linewidth', 2 );
grid on;
hold off;
axis equal;
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "|F1| vs x1, x2" );
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, abs(f2Mesh), 30 );
hold on;
plot( matXLevJTJ(1,:), matXLevJTJ(2,:), 'k-', 'linewidth', 2 );
plot( matXLevJTJB(1,:), matXLevJTJB(2,:), 'c-', 'linewidth', 2 );
plot( matXLevJTJC(1,:), matXLevJTJC(2,:), 'm-', 'linewidth', 2 );
plot( matXLev(1,:), matXLev(2,:), 'ro-' );
plot( matXLevB(1,:), matXLevB(2,:), 'gs-' );
plot( matXLevC(1,:), matXLevC(2,:), 'b^-' );
%%%plot( matXHOTGrad(1,:), matXHOTGrad(2,:), 'y*-' );
plot( matX_HOTGC(1,:), matX_HOTGC(2,:), 'y*-', 'linewidth', 2 );
grid on;
hold off;
axis equal;
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "|F2| vs x1, x2" );
return;
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, sqrt(omegaMesh), 30 );
axis equal;
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "sqrt(omega) vs x1, x2" );
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, f1Mesh, 30 );
axis equal;
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "F1 vs x1, x2" );
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, f2Mesh, 30 );
axis equal;
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "F2 vs x1, x2" );
