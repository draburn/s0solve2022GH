clear;
thisFile = "func1111_viz";
commondefs;
tic();
setprngstates(0);
numFigs = 0;
%
msg( thisFile, __LINE__, "This simple case illustrates some issues with existing handling of 'bad points'." );
%
sizeX = 2;
sizeF = 2;
isBadMin = true;
%
%%% THIS STUFF ISN'T USED!!!
funcPrm.ary3K = randn(sizeX,sizeX,sizeF);
matJPre = randn(sizeF,sizeX);
vecFE = randn(sizeF,1);
funcPrm.vecXE = randn(sizeX,1);
%
vecFEHat = vecFE/norm(vecFE);
for n=1:sizeX
for m=1:sizeX
	vecT1 = randn(sizeF,1);
	vecT2 = vecT1 - vecFEHat*(vecFEHat'*vecT1);
	vecT3 = abs(randn)*vecFEHat + vecT2;
	funcPrm.ary3K(n,m,:) = vecT3;
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
%
%%% THIS IS WHAT MATTERS!
funcPrm.vecXE = [ 0.0; 0.0 ];
funcPrm.ary3K(:,:,:) = 0.0;
funcPrm.ary3K(1,1,1) = 1.0;
%funcPrm.ary3K(2,2,1) = 1.0;
matJ = [0,0;0,1];
vecFE = [1;0];
funcPrm.matJ = matJ;
funcPrm.vecFE = vecFE;
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
vecGLev0 = [ vecXLev0(1) + 0.5*vecXLev0(1)^3; vecXLev0(2) ];
matHLev0 = [ 1.0+1.5*vecXLev0(1)^2, 0.0; 0.0, 1.0 ];
vecXLevJTJ0 = [ 4.0; 4.0 ];
vecFLevJTJ0 = func1111_eval( vecXLevJTJ0, funcPrm );
matJLevJTJ0 = [ vecXLevJTJ0(1), 0.0; 0.0, 1.0 ];
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
vecGLevB0 = [ vecXLevB0(1) + 0.5*vecXLevB0(1)^3; vecXLevB0(2) ];
matHLevB0 = [ 1.0+1.5*vecXLevB0(1)^2, 0.0; 0.0, 1.0 ];
vecXLevJTJB0 = matXLevJTJ(:,end);
vecFLevJTJB0 = func1111_eval( vecXLevJTJB0, funcPrm );
matJLevJTJB0 = [ vecXLevJTJB0(1), 0.0; 0.0, 1.0 ];
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
vecXLevJTJC0 = matXLevJTJB(:,end);
vecFLevJTJC0 = func1111_eval( vecXLevJTJC0, funcPrm );
matJLevJTJC0 = [ vecXLevJTJC0(1), 0.0; 0.0, 1.0 ];
vecGLevJTJC0 = matJLevJTJC0'*vecFLevJTJC0; % Should be same as vecGLev0 if at same pt.
matHLevJTJC0 = matJLevJTJC0'*matJLevJTJC0;
vecXLevC0 = matXLevB(:,end);
vecGLevC0 = [ vecXLevC0(1) + 0.5*vecXLevC0(1)^3; vecXLevC0(2) ];
matHLevC0 = [ 1.0+1.5*vecXLevC0(1)^2, 0.0; 0.0, 1.0 ];
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
for n=1:numGradSteps
	[ vecFTemp, matJTemp ] = func1111_evalDeriv( vecXHOTGrad, funcPrm );
	vecGTemp = matJTemp'*vecFTemp;
	vecXHOTGrad -= 0.01*vecGTemp;
	matXHOTGrad(:,n+1) = vecXHOTGrad;
end
%
%
%
numX1Vals = 51;
numX2Vals = 51;
x1Lo = -5.0;
x1Hi = +5.0;
x2Lo = -5.0;
x2Hi = +5.0;
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
contourf( x1Mesh, x2Mesh, log(omegaMesh), 30 );
hold on;
plot( matXLevJTJ(1,:), matXLevJTJ(2,:), 'k-', 'linewidth', 2 );
plot( matXLevJTJB(1,:), matXLevJTJB(2,:), 'c-', 'linewidth', 2 );
plot( matXLevJTJC(1,:), matXLevJTJC(2,:), 'm-', 'linewidth', 2 );
plot( matXLev(1,:), matXLev(2,:), 'ro-' );
plot( matXLevB(1,:), matXLevB(2,:), 'gs-' );
plot( matXLevC(1,:), matXLevC(2,:), 'b^-' );
plot( matXHOTGrad(1,:), matXHOTGrad(2,:), 'y*-' );
grid on;
hold off;
axis equal;
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "log(omega) vs x1, x2" );
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, abs(f1Mesh), 30 );
hold on;
plot( matXLevJTJ(1,:), matXLevJTJ(2,:), 'k-', 'linewidth', 2 );
plot( matXLevJTJB(1,:), matXLevJTJB(2,:), 'c-', 'linewidth', 2 );
plot( matXLevJTJC(1,:), matXLevJTJC(2,:), 'm-', 'linewidth', 2 );
plot( matXLev(1,:), matXLev(2,:), 'ro-' );
plot( matXLevB(1,:), matXLevB(2,:), 'gs-' );
plot( matXLevC(1,:), matXLevC(2,:), 'b^-' );
plot( matXHOTGrad(1,:), matXHOTGrad(2,:), 'y*-' );
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
plot( matXHOTGrad(1,:), matXHOTGrad(2,:), 'y*-' );
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
