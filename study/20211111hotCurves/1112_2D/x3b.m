thisFile = "x3b";
numFigs = 0;
tic();
%
if (1)
	[ matPsi, matLam ] = eig(testFuncPrm.matJ);
	[ foo, indexOfAbsMin ] = min(abs(diag(matLam)));
	clear foo;
	lamAbsMin = matLam(indexOfAbsMin,indexOfAbsMin);
	msg( thisFile, __LINE__, sprintf( "Absolute minimum eigenvalue: %g.", lamAbsMin ) );
	sizeU = 1;
	%matU = randn(sizeX,sizeU);
	matU = matPsi(:,indexOfAbsMin);
	matU(:,1) /= norm( matU(:,1) );
	%matU(:,2) -= matU(:,1)*(matU(:,1)'*matU(:,2));
	%matU(:,2) /= norm( matU(:,2) );
	%echo__matU = matU
	%echo__utu = matU'*matU
	%echo__uut = matU*(matU')
end
%
sizeV = sizeX;
matV = eye(sizeX,sizeX);
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
kSymNorm = sqrt(kSymNorm)
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
vecX0 = [ 4.0; 4.0 ];
ax = [ -5.0, 5.0, -5.0, 5.0 ];
%ax = [ -1.226, -1.224, 0.763, 0.765 ];
%vecX0 = [ 2.0; 1.0 ];
%ax = [ 0.0, 3.0, -1.0, 2.0 ];
%
%
%
funchF = @(x)( testFunc_eval(x,testFuncPrm) );
%
[ vecF0, matJ0, ary3K0 ] = calcDeriv( funchF, vecX0 );
modelFuncPrm.sizeX = sizeX;
modelFuncPrm.sizeF = sizeF;
modelFuncPrm.vecXE = vecX0;
modelFuncPrm.vecFE = vecF0;
modelFuncPrm.matJ = matJ0;
modelFuncPrm.ary3K = ary3K0;
modelFuncPrm0 = modelFuncPrm;
matYGrad_full = calcHOTGradCurve( modelFuncPrm, vecX0 );
%
modelFuncPrm.ary3K = ary3UUTKUUT;
matYGrad_part = calcHOTGradCurve( modelFuncPrm, vecX0 );
%
modelFuncPrm.ary3K(:,:,:) = 0.0;
matYGrad_cold = calcHOTGradCurve( modelFuncPrm, vecX0 );
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
matFModel = testFunc_eval(matX,modelFuncPrm);
f1ModelMesh = reshape(matFModel(1,:),numX2Vals,numX1Vals);
f2ModelMesh = reshape(matFModel(2,:),numX2Vals,numX1Vals);
omegaModelMesh = 0.5*( f1ModelMesh.^2 + f2ModelMesh.^2 );
%
msg( thisFile, __LINE__, sprintf( "F1 scale: %g to %g.", min(min((f1Mesh))), max(max((f1Mesh))) ) );
msg( thisFile, __LINE__, sprintf( "F2 scale: %g to %g.", min(min((f2Mesh))), max(max((f2Mesh))) ) );
msg( thisFile, __LINE__, sprintf( "omega scale: %g to %g.", min(min((omegaMesh))), max(max((omegaMesh))) ) );
msg( thisFile, __LINE__, sprintf( "omegaModel scale: %g to %g.", min(min((omegaModelMesh))), max(max((omegaModelMesh))) ) );
%
%
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, log(eps025*max(max(omegaMesh))+omegaMesh-min(min(omegaMesh))), 30 );
hold on;
plot( matYGrad_full(1,:), matYGrad_full(2,:), 'r*-', 'linewidth', 2 );
plot( matYGrad_part(1,:), matYGrad_part(2,:), 'g*-', 'linewidth', 2 );
plot( matYGrad_cold(1,:), matYGrad_cold(2,:), 'b*-', 'linewidth', 2 );
hold off;
grid on;
axis equal;
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "log(omega-omegaMin) vs x1, x2" );
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, log(eps025*max(max(omegaModelMesh))+omegaModelMesh-min(min(omegaModelMesh))), 30 );
hold on;
plot( matYGrad_full(1,:), matYGrad_full(2,:), 'r*-', 'linewidth', 2 );
plot( matYGrad_part(1,:), matYGrad_part(2,:), 'g*-', 'linewidth', 2 );
plot( matYGrad_cold(1,:), matYGrad_cold(2,:), 'b*-', 'linewidth', 2 );
hold off;
grid on;
axis equal;
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "MODEL log(omega-omegaMin) vs x1, x2" );
%
toc();
