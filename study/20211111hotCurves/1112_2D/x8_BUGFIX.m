clear;
thisFile = "x6_diakonal";
commondefs;
numFigs = 0;
tic();
%
useAxisEqual = false;
sizeX = 2;
sizeF = 2;
if (1)
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,97889488); % x5 "Snakey"
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,31039344); % x5 "Both bad".
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,4924368); % x5 Another "green good".
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,34123136); % x6 green has x-pt.
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,33394912,true,false,false);
	%%%testFuncPrm = testFunc_genPrm(sizeX,sizeF,60525840,true,false,false); % x8 Green bad!
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,45342256,true,false,false); % x8 Curves 2 are short.
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,60663152,true,false,false); % x8 red misses.
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,91475584,true,false,false); % Dramatic benefit.
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,72643728,true,false,false); %Hmm...
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,12186496,true,false,false);
	testFuncPrm = testFunc_genPrm(sizeX,sizeF,45553360,true,false,false); % Lev and Grad cnvg to different!
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,85541360,true,false,false); % Lev bonks.
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,81404512,true,false,false); % Lev bonks HARD.
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,82178352,true,false,false); % First Lev step is big.
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,96042480,true,false,false); % First Lev step is funky.
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,53162192,true,false,false); % Lev bonks.
	vecX0 = [ 3.0; 3.0 ];
	ax = [ -5.0, 5.0, -5.0, 5.0 ];
	%ax = [ -0.755, -0.753, -1.934, -1.931 ];
	%ax = [ 1.876, 1.878, -1.175, -1.174 ]
	%ax = [ 0.3, 0.34, 0.6, 0.64 ];
	%ax = [ 2.4, 2.8, 1.0, 1.5 ];
elseif (0)
	% Actual F is uniminal but (standard-basis-) diagonal K isn't.
	msg( thisFile, __LINE__, "*** CASE: Actual func is uniminal(?) but e1 'partial' is not." );
	msg( thisFile, __LINE__, "(Is fine when using deltaN, however.)" );
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	testFuncPrm.vecFE = [ 2.0; 2.0 ];
	testFuncPrm.matJ = [ 1.0, 0.0; 0.0, 1.0 ];
	testFuncPrm.ary3K = zeros(sizeX,sizeX,sizeF);
	testFuncPrm.ary3K(:,:,1) = [ 2.0, 2.0; 2.0, 2.0 ];
	testFuncPrm.ary3K(:,:,2) = [ 2.0, 2.0; 2.0, 2.0 ];
	vecX0 = [ 4.0; 4.0 ];
	ax = 5*[-1,1,-1,1];
elseif (0)
	msg( thisFile, __LINE__, "*** CASE: Make actual F non-uniminal." );
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	%testFuncPrm.vecFE = [ 70.0; 70.0 ];
	%testFuncPrm.matJ = [ 17.0, 16.0; 16.0, 17.0 ];
	testFuncPrm.vecFE = [ 2.0; 2.0 ];
	testFuncPrm.matJ = 10.0*[ 1.0, 0.0; 0.0, 1.0 ];
	testFuncPrm.ary3K = zeros(sizeX,sizeX,sizeF);
	testFuncPrm.ary3K(:,:,1) = [ 2.0, 0.0; 0.0, 2.0 ];
	testFuncPrm.ary3K(:,:,2) = [ 2.0, 0.0; 0.0, 2.0 ];
	vecX0 = [ 3.0; 0.0 ];
	ax = 5*[-1,1,-1,1];
elseif (1)
	msg( thisFile, __LINE__, "*** CASE: Illustration of bad local min original motivation." );
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	testFuncPrm.vecFE = [ 1.0; 0.0 ];
	testFuncPrm.matJ = [ 0.0, 0.0; 0.0, 1.0 ];
	testFuncPrm.ary3K = zeros(sizeX,sizeX,sizeF);
	testFuncPrm.ary3K(:,:,1) = [ 2.0, 0.5; 0.5, 0.0 ];
	testFuncPrm.ary3K(:,:,2) = [ 2.0, 0.5; 0.5, 0.0 ];
	vecX0 = [ 0.2; 0.2 ];
	if (0)
		testFuncPrm.ary3K(:,:,1) = [ 2.0, 0.5; 0.5, 2.0 ];
		testFuncPrm.ary3K(:,:,2) = [ 2.0, 0.5; 0.5, 2.0 ];
		vecX0 = [ 0.5; 0.5 ];
	end
	ax = [-3,1,-1,2];
elseif (1)
	error( "THIS CASE IS STALE!" );
	% Actual F is uniminal but eigenspace diagonal K isn't.
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	testFuncPrm.vecFE = [ 1.0; 0.0 ];
	%%%testFuncPrm.matJ = [ 0.0, 1.0; 1.0, 0.0 ];
	testFuncPrm.matJ = [ 1.0, 0.0; 0.0, 1.0 ];
	testFuncPrm.ary3K = zeros(sizeX,sizeX,sizeF);
	%%%testFuncPrm.ary3K(:,:,1) = [ -1.0, 0.0; 0.0, -0.5 ];
	%%%testFuncPrm.ary3K(:,:,2) = [  1.0, 0.0; 0.0, 0.5 ];
	testFuncPrm.ary3K(:,:,1) = [ 2.0, 0.0; 0.0, 2.0 ];
	testFuncPrm.ary3K(:,:,2) = [ 0.0, 0.0; 0.0, 0.0 ];
	vecX0 = [ 0.0; 4.0 ];
	ax = [ -5.0, 5.0, -5.0, 5.0 ];
elseif (1)
	error( "THIS CASE IS STALE!" );
	% Actual F is non-uniminal.
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	testFuncPrm.vecFE = [ 0.0; 0.0 ];
	testFuncPrm.matJ = [ 1.0, 0.0; 0.0, 1.0 ];
	testFuncPrm.ary3K = zeros(sizeX,sizeX,sizeF);
	testFuncPrm.ary3K(:,:,1) = [ 2.0, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,2) = [ 0.0, 0.0; 0.0, 0.0 ];
	vecX0 = [ 0.6; 1.0 ];
	ax = [ -1.0, 2.0, -1.0, 2.0 ];
else
	error( "THIS CASE IS STALE!" );
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	testFuncPrm.vecFE = [ 1.0; 0.0 ];
	testFuncPrm.matJ = [ 0.0, 0.0; 0.0, 1.0 ];
	testFuncPrm.ary3K = zeros(sizeX,sizeX,sizeF);
	testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.1 ];
	testFuncPrm.ary3K(:,:,2) = [ 2.0, 0.0; 0.0, 0.0 ];
end
%
funchF = @(x)( testFunc_eval(x,testFuncPrm) );
%
[ vecF0, matJ0, ary3K0 ] = calcDeriv( funchF, vecX0 );
%!!!!!
% PREVIOUS CODES HAD matJ0.*matJ0 instead!!!!!
[ matPsi0, matLam0 ] = eig(matJ0'*matJ0);
%!!!!!
%matPsi0'*matPsi0
%matPsi0*(matPsi0')
matU = matPsi0;
if (0)
	msg( thisFile, __LINE__, sprintf( "USING V = e1." ) );
	matV = [1;0];
elseif (0)
	msg( thisFile, __LINE__, sprintf( "USING V = g." ) );
	vecG0 = matJ0'*vecF0;
	matV = vecG0/norm(vecG0);
elseif (1)
	msg( thisFile, __LINE__, sprintf( "USING V = deltaN." ) );
	matV = -matJ0\vecF0;
	matV /= norm(matV);
elseif (0)
	msg( thisFile, __LINE__, sprintf( "USING V = psiCrit0." ) );
	[ lamAbsMin, indexOfAbsMin ] = min(abs(diag(matLam0)));
	lamAbsMax = max(abs(diag(matLam0)));
	msg( thisFile, __LINE__, sprintf( "Absolute eigenvalue range at X0: %g ~ %g.", lamAbsMin, lamAbsMax ) );
	matV = matPsi0(:,indexOfAbsMin);
else
	msg( thisFile, __LINE__, sprintf( "USING V = psiCritE." ) );
	[ matPsiE, matLamE ] = eig(testFuncPrm.matJ'*testFuncPrm.matJ);
	[ lamAbsMin, indexOfAbsMin ] = min(abs(diag(matLamE)));
	lamAbsMax = max(abs(diag(matLamE)));
	msg( thisFile, __LINE__, sprintf( "Absolute eigenvalue range at XE: %g ~ %g.", lamAbsMin, lamAbsMax ) );
	assert( lamAbsMin < (eps^0.75)*lamAbsMax );
	matV = matPsiE(:,indexOfAbsMin);
end
[ matJU, matJV, ary3UTKU, ary3UTKV, ary3VTKU, ary3VTKV ] = calcSubSpaceDeriv( ...
  funchF, vecX0, matU, matV );
matJUUT = matJU*(matU');
matJVVT = matJV*(matV');
for n=1:sizeF
	ary3UUTKUUT(:,:,n) = matU*ary3UTKU(:,:,n)*(matU');
	ary3UUTKVVT(:,:,n) = matU*ary3UTKV(:,:,n)*(matV');
	ary3VVTKUUT(:,:,n) = matV*ary3VTKU(:,:,n)*(matU');
	ary3VVTKVVT(:,:,n) = matV*ary3VTKV(:,:,n)*(matV');
	%
	ary3UTKUDiag(:,:,n) = diag(diag(ary3UTKU(:,:,n)));
	ary3VTKVDiag(:,:,n) = diag(diag(ary3VTKV(:,:,n)));
	ary3UUTKUUTDiag(:,:,n) = matU*diag(diag(ary3UTKU(:,:,n)))*(matU');
	ary3VVTKVVTDiag(:,:,n) = matV*diag(diag(ary3VTKV(:,:,n)))*(matV');
end
%
kSymNorm = 0.0;
for n=1:sizeF
	kSymNorm += sum(sum((ary3UTKV(:,:,n)'-ary3VTKU(:,:,n)).^2));
end
kSymNorm = sqrt(kSymNorm);
msg( thisFile, __LINE__, sprintf( "U norms: %10.3e, %10.3e.", ...
  sqrt(sum(sum((matU'*matU-eye(2,2)).^2))), ...
  sqrt(sum(sum((matU*(matU')-eye(2,2)).^2))) ) );
msg( thisFile, __LINE__, sprintf( "V norms: %10.3e, %10.3e.", ...
  sqrt(sum(sum((matV'*matV-eye(size(matV'*matV))).^2))), ...
  sqrt(sum(sum((matV*(matV')-eye(2,2)).^2))) ) );
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
localFuncPrm_full.sizeX = sizeX;
localFuncPrm_full.sizeF = sizeF;
localFuncPrm_full.vecXE = vecX0;
localFuncPrm_full.vecFE = vecF0;
localFuncPrm_full.matJ = matJ0;
localFuncPrm_full.ary3K = ary3K0;
matYGrad_full = calcHOTGradCurveRK4( localFuncPrm_full, vecX0 );
%
localFuncPrm_part = localFuncPrm_full;
%%%localFuncPrm_part.ary3K(:,:,:) = 0.0;
%%%localFuncPrm_part.ary3K = ary3VVTKVVTDiag;
[ matYGrad_part, vecSGrad_part ] = calcHOTLevCurveBound4( localFuncPrm_part, vecX0 );
numPts_part = size(matYGrad_part,2)
localFuncPrm_part.ary3K(:,:,:) = 0.0;
localFuncPrm_part.ary3K = ary3VVTKVVTDiag;
%
localFuncPrm_cold = localFuncPrm_full;
%%%localFuncPrm_cold.ary3K(:,:,:) = 0.0;
%%%matYGrad_cold = calcHOTGradCurveRK4( localFuncPrm_cold, vecX0 );
matYGrad_cold = calcHOTLevCurveRK4( localFuncPrm_cold, vecX0 );
if (0)
funch_gHat = @(x)( -testFunc_gHat(x,localFuncPrm_cold) );
matYGrad_cold = lsode( funch_gHat, vecX0, 10.0*linspace(0.0,1.0,101) );
matYGrad_cold = matYGrad_cold';
end
localFuncPrm_cold.ary3K(:,:,:) = 0.0;
%
%%%%%%%%%%%%%%%%%
% VARIOUS MISC WORK
numPts_local_cold = size(matYGrad_cold,2);
vecFPts_local_cold = testFunc_eval(matYGrad_cold,localFuncPrm_cold);
omegaPts_local_cold = 0.5*sum(vecFPts_local_cold.^2,1);
totDistPts_local_cold = sqrt(sum((matYGrad_cold-repmat(vecX0,1,numPts_local_cold)).^2));
%
matHdiakonal = zeros(sizeX,sizeX);
for n=1:sizeF
	matHdiakonal += localFuncPrm_part.vecFE(n) * localFuncPrm_part.ary3K(:,:,n);
end
%
%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
numX1Vals = 41;
numX2Vals = 41;
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
plot( matYGrad_full(1,:), matYGrad_full(2,:), 'rx-', 'linewidth', 2, 'markersize', 5 );
plot( matYGrad_part(1,:), matYGrad_part(2,:), 'go-', 'linewidth', 2, 'markersize', 5 );
plot( matYGrad_cold(1,:), matYGrad_cold(2,:), 'b+-', 'linewidth', 2, 'markersize', 5 );
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
plot( matYGrad_full(1,:), matYGrad_full(2,:), 'rx-', 'linewidth', 2, 'markersize', 5 );
plot( matYGrad_part(1,:), matYGrad_part(2,:), 'go-', 'linewidth', 2, 'markersize', 5 );
plot( matYGrad_cold(1,:), matYGrad_cold(2,:), 'b+-', 'linewidth', 2, 'markersize', 5 );
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
plot( matYGrad_full(1,:), matYGrad_full(2,:), 'rx-', 'linewidth', 2, 'markersize', 5 );
plot( matYGrad_part(1,:), matYGrad_part(2,:), 'go-', 'linewidth', 2, 'markersize', 5 );
plot( matYGrad_cold(1,:), matYGrad_cold(2,:), 'b+-', 'linewidth', 2, 'markersize', 5 );
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
plot( matYGrad_full(1,:), matYGrad_full(2,:), 'rx-', 'linewidth', 2, 'markersize', 5 );
plot( matYGrad_part(1,:), matYGrad_part(2,:), 'go-', 'linewidth', 2, 'markersize', 5 );
plot( matYGrad_cold(1,:), matYGrad_cold(2,:), 'b+-', 'linewidth', 2, 'markersize', 5 );
hold off;
grid on;
if (useAxisEqual)
	axis equal;
end
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "LOCAL COLD log(omega-omegaMin) vs x1, x2" );

%
numFigs++; figure(numFigs);
plot( vecSGrad_part, 'o-' );
grid on;
numFigs++; figure(numFigs);
plot( matYGrad_part(1,:), 'o-' );
grid on;
numFigs++; figure(numFigs);
plot( matYGrad_part(2,:), 'o-' );
grid on;

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
