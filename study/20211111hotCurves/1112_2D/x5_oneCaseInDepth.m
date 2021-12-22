clear;
thisFile = "x5_oneCaseInDepth";
commondefs;
numFigs = 0;
tic();
%
useAxisEqual = false;
sizeX = 2;
sizeF = 2;
%testFuncPrm = testFunc_genPrm(sizeX,sizeF,97889488); % "Snakey"
%%%testFuncPrm = testFunc_genPrm(sizeX,sizeF,31039344); % "Both bad".
%testFuncPrm = testFunc_genPrm(sizeX,sizeF,4924368); % Another "green good".
testFuncPrm.sizeX = 2;
testFuncPrm.sizeF = 2;
testFuncPrm.vecXE = [ 0.0; 0.0 ];
testFuncPrm.vecFE = [ 1.0; 0.0 ];
testFuncPrm.matJ = [ 0.0, 0.0; 0.0, 1.0 ];
testFuncPrm.ary3K = zeros(sizeX,sizeX,sizeF);
testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.1 ];
testFuncPrm.ary3K(:,:,2) = [ 2.0, 0.0; 0.0, 0.0 ];
%
funchF = @(x)( testFunc_eval(x,testFuncPrm) );
%vecX0 = [ 1.0; 1.0 ];
%vecX0 = [ 0.01; 1.0 ];
vecX0 = [ 0.1; 0.1 ];
%vecX0 = [ -1; 0 ]
%ax = 1*[ -5.0, 5.0, -5.0, 5.0 ];
%ax = 1*[ -5.0, 5.0, -5.0, 5.0 ];
ax = [ -0.1, 3.0, -0.6, 0.2 ];
%ax = [ -1.0, 2.0, -1.0, 2.0 ];
%ax = [ -0.1, 0.2, -0.2, 0.1 ];
%
[ matPsi, matLam ] = eig(testFuncPrm.matJ);
[ lamAbsMin, indexOfAbsMin ] = min(abs(diag(matLam)));
lamAbsMax = max(abs(diag(matLam)));
msg( thisFile, __LINE__, sprintf( "Absolute eigenvalue range: %g ~ %g.", lamAbsMin, lamAbsMax ) );
assert( lamAbsMin < (eps^0.75)*lamAbsMax );
sizeU = 1;
matU = matPsi(:,indexOfAbsMin);
matU(:,1) /= norm( matU(:,1) );
sizeV = 1;
matV = randn(sizeX,sizeV);
matV(:,1) -= matU(:,1)*( matU(:,1)'*matV(:,1) );
matV(:,1) /= norm( matV(:,1) );

%%%msg( thisFile, __LINE__, "~~~ MODIFYING matJ AFTER CREATING matU, matV! ~~~" );
%%%testFuncPrm.matJ = [ 0.1, 0.0; 0.0, 1.0 ];
%%%funchF = @(x)( testFunc_eval(x,testFuncPrm) );
%
%%%matV = [ matU matV ] % Useful to check this is fullspace.
%
[ vecF0, matJ0, ary3K0 ] = calcDeriv( funchF, vecX0 );
%
[ matJU, matJV, ary3UTKU, ary3UTKV, ary3VTKU, ary3VTKV ] = calcSubSpaceDeriv( ...
  funchF, vecX0, matU, matV );
%
matJUUT = matJU*(matU');
matJVVT = matJV*(matV');
for n=1:sizeF
	ary3UUTKUUT(:,:,n) = matU*ary3UTKU(:,:,n)*(matU');
	ary3UUTKVVT(:,:,n) = matU*ary3UTKV(:,:,n)*(matV');
	ary3VVTKUUT(:,:,n) = matV*ary3VTKU(:,:,n)*(matU');
	ary3VVTKVVT(:,:,n) = matV*ary3VTKV(:,:,n)*(matV');
end
%
kSymNorm = 0.0;
for n=1:sizeF
	kSymNorm += sum(sum((ary3UTKV(:,:,n)'-ary3VTKU(:,:,n)).^2));
end
kSymNorm = sqrt(kSymNorm);
msg( thisFile, __LINE__, sprintf( "K sym norm: %10.3e.", kSymNorm ) );
msg( thisFile, __LINE__, sprintf( "J norms: %10.3e, %10.3e.", ...
  sqrt(sum(sum((matJUUT - matJ0).^2))), ...
  sqrt(sum(sum((matJVVT - matJ0).^2))) ) );
msg( thisFile, __LINE__, sprintf( "K norms: %10.3e, %10.3e, %10.3e, %10.3e.", ...
  sqrt(sum(sum(sum((ary3UUTKUUT - ary3K0).^2)))), ...
  sqrt(sum(sum(sum((ary3UUTKVVT - ary3K0).^2)))), ...
  sqrt(sum(sum(sum((ary3VVTKUUT - ary3K0).^2)))), ...
  sqrt(sum(sum(sum((ary3VVTKVVT - ary3K0).^2)))) ) );
%
%
%
[ vecF0, matJ0, ary3K0 ] = calcDeriv( funchF, vecX0 );
localFuncPrm_full.sizeX = sizeX;
localFuncPrm_full.sizeF = sizeF;
localFuncPrm_full.vecXE = vecX0;
localFuncPrm_full.vecFE = vecF0;
localFuncPrm_full.matJ = matJ0;
localFuncPrm_full.ary3K = ary3K0;
matYGrad_full = calcHOTGradCurve( localFuncPrm_full, vecX0 );
%
localFuncPrm_part = localFuncPrm_full;
localFuncPrm_part.ary3K = ary3UUTKUUT + ary3VVTKUUT + ary3UUTKVVT;
matYGrad_part = calcHOTGradCurve( localFuncPrm_part, vecX0 );
%
localFuncPrm_cold = localFuncPrm_full;
localFuncPrm_cold.ary3K(:,:,:) = 0.0;
matYGrad_cold = calcHOTGradCurve( localFuncPrm_cold, vecX0 );
%
numPts_local_cold = size(matYGrad_cold,2);
vecFPts_local_cold = testFunc_eval(matYGrad_cold,localFuncPrm_cold);
omegaPts_local_cold = 0.5*sum(vecFPts_local_cold.^2,1);
totDistPts_local_cold = sqrt(sum((matYGrad_cold-repmat(vecX0,1,numPts_local_cold)).^2));
%
%
%
%
%
%
numX1Vals = 51;
numX2Vals = 61;
x1Vals = linspace(ax(1),ax(2),numX1Vals);
x2Vals = linspace(ax(3),ax(4),numX2Vals);
[ x1Mesh, x2Mesh ] = meshgrid( x1Vals, x2Vals );
matX = [ reshape(x1Mesh,1,[]); reshape(x2Mesh,1,[]) ];
matF = funchF(matX,testFuncPrm);
f1Mesh = reshape(matF(1,:),numX2Vals,numX1Vals);
f2Mesh = reshape(matF(2,:),numX2Vals,numX1Vals);
omegaMesh = 0.5*( f1Mesh.^2 + f2Mesh.^2 );
%
matF_local_full = testFunc_eval(matX,localFuncPrm_full);
f1Mesh_local_full = reshape(matF_local_full(1,:),numX2Vals,numX1Vals);
f2Mesh_local_full = reshape(matF_local_full(2,:),numX2Vals,numX1Vals);
omegaMesh_local_full = 0.5*( f1Mesh_local_full.^2 + f2Mesh_local_full.^2 );
%
matF_local_part = testFunc_eval(matX,localFuncPrm_part);
f1Mesh_local_part = reshape(matF_local_part(1,:),numX2Vals,numX1Vals);
f2Mesh_local_part = reshape(matF_local_part(2,:),numX2Vals,numX1Vals);
omegaMesh_local_part = 0.5*( f1Mesh_local_part.^2 + f2Mesh_local_part.^2 );
%
matF_local_cold = testFunc_eval(matX,localFuncPrm_cold);
f1Mesh_local_cold = reshape(matF_local_cold(1,:),numX2Vals,numX1Vals);
f2Mesh_local_cold = reshape(matF_local_cold(2,:),numX2Vals,numX1Vals);
omegaMesh_local_cold = 0.5*( f1Mesh_local_cold.^2 + f2Mesh_local_cold.^2 );
%
msg( thisFile, __LINE__, sprintf( "F1 scale: %g to %g.", min(min((f1Mesh))), max(max((f1Mesh))) ) );
msg( thisFile, __LINE__, sprintf( "F2 scale: %g to %g.", min(min((f2Mesh))), max(max((f2Mesh))) ) );
msg( thisFile, __LINE__, sprintf( "omega scale: %g to %g.", min(min((omegaMesh))), max(max((omegaMesh))) ) );
msg( thisFile, __LINE__, sprintf( "LOCAL FULL omegaModel scale: %g to %g.", min(min(omegaMesh_local_full)), max(max(omegaMesh_local_full)) ) );
msg( thisFile, __LINE__, sprintf( "LOCAL PART omegaModel scale: %g to %g.", min(min(omegaMesh_local_part)), max(max(omegaMesh_local_part)) ) );
msg( thisFile, __LINE__, sprintf( "LOCAL COLD omegaModel scale: %g to %g.", min(min(omegaMesh_local_cold)), max(max(omegaMesh_local_cold)) ) );
%
funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
%
%
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, funchVizLog(omegaMesh), 30 );
hold on;
plot( matYGrad_full(1,:), matYGrad_full(2,:), 'r*-', 'linewidth', 2 );
plot( matYGrad_part(1,:), matYGrad_part(2,:), 'g*-', 'linewidth', 2 );
plot( matYGrad_cold(1,:), matYGrad_cold(2,:), 'b*-', 'linewidth', 2 );
hold off;
grid on;
if (useAxisEqual)
	axis equal;
end
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "log(omega-omegaMin) vs x1, x2" );
%
if (0)
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, funchVizLog(omegaMesh_local_full), 30 );
hold on;
plot( matYGrad_full(1,:), matYGrad_full(2,:), 'r*-', 'linewidth', 2 );
plot( matYGrad_part(1,:), matYGrad_part(2,:), 'g*-', 'linewidth', 2 );
plot( matYGrad_cold(1,:), matYGrad_cold(2,:), 'b*-', 'linewidth', 2 );
hold off;
grid on;
if (useAxisEqual)
	axis equal;
end
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "LOCAL FULL log(omega-omegaMin) vs x1, x2" );
end
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, funchVizLog(omegaMesh_local_part), 30 );
hold on;
plot( matYGrad_full(1,:), matYGrad_full(2,:), 'r*-', 'linewidth', 2 );
plot( matYGrad_part(1,:), matYGrad_part(2,:), 'g*-', 'linewidth', 2 );
plot( matYGrad_cold(1,:), matYGrad_cold(2,:), 'b*-', 'linewidth', 2 );
hold off;
grid on;
if (useAxisEqual)
	axis equal;
end
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "LOCAL PART log(omega-omegaMin) vs x1, x2" );
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, funchVizLog(omegaMesh_local_cold), 30 );
hold on;
plot( matYGrad_full(1,:), matYGrad_full(2,:), 'r*-', 'linewidth', 2 );
plot( matYGrad_part(1,:), matYGrad_part(2,:), 'g*-', 'linewidth', 2 );
plot( matYGrad_cold(1,:), matYGrad_cold(2,:), 'b*-', 'linewidth', 2 );
hold off;
grid on;
if (useAxisEqual)
	axis equal;
end
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "LOCAL COLD log(omega-omegaMin) vs x1, x2" );

return


%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, abs(f1Mesh), 30 );
hold on;
plot( matYGrad_full(1,:), matYGrad_full(2,:), 'r*-', 'linewidth', 2 );
plot( matYGrad_part(1,:), matYGrad_part(2,:), 'g*-', 'linewidth', 2 );
plot( matYGrad_cold(1,:), matYGrad_cold(2,:), 'b*-', 'linewidth', 2 );
hold off;
grid on;
if (useAxisEqual)
	axis equal;
end
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "|F1| vs x1, x2" );
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, abs(f1Mesh_local_part), 30 );
hold on;
plot( matYGrad_full(1,:), matYGrad_full(2,:), 'r*-', 'linewidth', 2 );
plot( matYGrad_part(1,:), matYGrad_part(2,:), 'g*-', 'linewidth', 2 );
plot( matYGrad_cold(1,:), matYGrad_cold(2,:), 'b*-', 'linewidth', 2 );
hold off;
grid on;
if (useAxisEqual)
	axis equal;
end
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "LOCAL PART |F1| vs x1, x2" );
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, abs(f1Mesh_local_cold), 30 );
hold on;
plot( matYGrad_full(1,:), matYGrad_full(2,:), 'r*-', 'linewidth', 2 );
plot( matYGrad_part(1,:), matYGrad_part(2,:), 'g*-', 'linewidth', 2 );
plot( matYGrad_cold(1,:), matYGrad_cold(2,:), 'b*-', 'linewidth', 2 );
hold off;
grid on;
if (useAxisEqual)
	axis equal;
end
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "LOCAL COLD |F1| vs x1, x2" );
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, abs(f2Mesh), 30 );
hold on;
plot( matYGrad_full(1,:), matYGrad_full(2,:), 'r*-', 'linewidth', 2 );
plot( matYGrad_part(1,:), matYGrad_part(2,:), 'g*-', 'linewidth', 2 );
plot( matYGrad_cold(1,:), matYGrad_cold(2,:), 'b*-', 'linewidth', 2 );
hold off;
grid on;
if (useAxisEqual)
	axis equal;
end
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "|F2| vs x1, x2" );
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, abs(f2Mesh_local_part), 30 );
hold on;
plot( matYGrad_full(1,:), matYGrad_full(2,:), 'r*-', 'linewidth', 2 );
plot( matYGrad_part(1,:), matYGrad_part(2,:), 'g*-', 'linewidth', 2 );
plot( matYGrad_cold(1,:), matYGrad_cold(2,:), 'b*-', 'linewidth', 2 );
hold off;
grid on;
if (useAxisEqual)
	axis equal;
end
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "LOCAL PART |F2| vs x1, x2" );
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, abs(f2Mesh_local_cold), 30 );
hold on;
plot( matYGrad_full(1,:), matYGrad_full(2,:), 'r*-', 'linewidth', 2 );
plot( matYGrad_part(1,:), matYGrad_part(2,:), 'g*-', 'linewidth', 2 );
plot( matYGrad_cold(1,:), matYGrad_cold(2,:), 'b*-', 'linewidth', 2 );
hold off;
grid on;
if (useAxisEqual)
	axis equal;
end
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "LOCAL COLD |F2| vs x1, x2" );
%
toc();
