clear;
thisFile = "calcMinfordCurve__test1";
commondefs;
setprngstates(0);
numFigs = 0;
startTime = time();
%
contourPlot_numColors = 51;
contourPlot_numLevels = 31;
numX1Vals = 71;
numX2Vals = 81;
numThetaVals = 1001;
%
%
switch (10)
case 1
	sizeX = 2;
	funchOmega = @(x)(0.5*(x'*x));
	funchG = @(x)(x);
	funchH = @(x)(eye(2,2));
	vecX1Hat = [1;0];
	vecX2Hat = [0;1];
	%
	vecX0 = [eps;eps]
	vecXC = zeros(sizeX,1)
	bigR = 1.0
	matS = [2.0,0.0;0.0,1.0]
case 2
	sizeX = 2;
	matH0 = [ 1.0, 0.0; 0.0, 16.0 ];
	funchOmega = @(x)(0.5*(x'*matH0*x));
	funchG = @(x)(matH0*x);
	funchH = @(x)(matH0);
	vecX1Hat = [1;0];
	vecX2Hat = [0;1];
	%
	vecX0 = [eps;eps]
	vecXC = zeros(sizeX,1)
	bigR = 1.0
	matS = [4.0,0.0;0.0,1.0]
case 3
	sizeX = 2;
	matH0 = [ 1.0, 0.0; 0.0, 4.0 ];
	funchOmega = @(x)(0.5*(x'*matH0*x));
	funchG = @(x)(matH0*x);
	funchH = @(x)(matH0);
	vecX1Hat = [1;0];
	vecX2Hat = [0;1];
	%
	vecX0 = [eps;eps]
	vecXC = zeros(sizeX,1)
	bigR = 1.0
	matS = [1.0,0.0;0.0,2.0]
case 4
	foo = 5.0
	sizeX = 2;
	matH0 = [ 1.0, 0.0; 0.0, foo^2 ];
	funchOmega = @(x)(0.5*(x'*matH0*x));
	funchG = @(x)(matH0*x);
	funchH = @(x)(matH0);
	vecX1Hat = [1;0];
	vecX2Hat = [0;1];
	%
	vecX0 = [eps;eps]
	vecXC = zeros(sizeX,1)
	bigR = 1.0
	matS = [foo,0.0;0.0,1.0]
case 10
	sizeF = 2;
	sizeX = 2;
	omega0 = abs(randn());
	vecG0 = randn(sizeX,1);
	matHGen = randn(sizeF,sizeX);
	matH0 = matHGen'*matHGen;
	funchOmega = @(x)( omega0 + x'*vecG0 + 0.5*(x'*matH0*x) );
	funchG = @(x)( vecG0 + matH0*x );
	funchH = @(x)( matH0 );
	%
	vecX1Hat = [1;0];
	vecX2Hat = [0;1];
	%
	vecX0 = zeros(sizeX,1)
	vecXC = zeros(sizeX,1)
	bigR = 1.0
	matS = [2.0,0.0;0.0,1.0]
case 11
	sizeF = 3;
	sizeX = 3;
	omega0 = randn();;
	vecG0 = randn(sizeX,1);
	matHGen = randn(sizeF,sizeX);
	matH0 = matHGen'*matHGen;
	funchOmega = @(x)( omega0 + x'*vecG0 + 0.5*(x'*matH0*x) );
	funchG = @(x)( vecG0 + matH0*x );
	funchH = @(x)( matH0 );
	vecX1Hat = randn(sizeX,1);
	vecX2Hat = randn(sizeX,1);
	%
	vecX1Hat = randn(sizeX,1);
	vecX2Hat = randn(sizeX,1);
	%
	vecX0 = randn(sizeX,1)
	vecXC = randn(sizeX,1)
	bigR = abs(randn)
	matS = abs(randn(sizeX,sizeX))
case 20
	sizeF = 10;
	sizeX = 10;
	omega0 = randn();;
	vecG0 = randn(sizeX,1);
	matHGen = randn(sizeF,sizeX);
	matH0 = matHGen'+matHGen;
	funchOmega = @(x)( omega0 + x'*vecG0 + 0.5*(x'*matH0*x) );
	funchG = @(x)( vecG0 + matH0*x );
	funchH = @(x)( matH0 );
	%
	vecX1Hat = randn(sizeX,1);
	vecX2Hat = randn(sizeX,1);
	%
	vecX0 = randn(sizeX,1)
	vecXC = randn(sizeX,1)
	bigR = abs(randn)
	matS = abs(randn(sizeX,sizeX))
otherwise
	error( "Invalid case." );
end
%
%
%
vecX1Hat /= norm(vecX1Hat);
vecX1Hat /= norm(vecX1Hat)
vecX2Hat -= vecX1Hat*(vecX1Hat'*vecX2Hat);
vecX2Hat /= norm(vecX2Hat);
vecX2Hat -= vecX1Hat*(vecX1Hat'*vecX2Hat);
vecX2Hat /= norm(vecX2Hat)
%
x1Vals = 3.0*linspace(-1.0,1.0,numX1Vals);
x2Vals = 3.0*linspace(-1.0,1.0,numX2Vals);
%x1Vals = linspace(1.0,1.4,numX1Vals);
%x2Vals = linspace(1.0,1.15,numX2Vals);
for i1=1:numX1Vals
for i2=1:numX2Vals
	vecX = vecX0 + vecX1Hat*x1Vals(i1) + vecX2Hat*x2Vals(i2);
	x1Mesh(i1,i2) = x1Vals(i1);
	x2Mesh(i1,i2) = x2Vals(i2);
	omegaMesh(i1,i2) = funchOmega(vecX);
	omegaSurfMesh(i1,i2) = calcMinfordCurve__evalOmegaSurf( funchOmega, vecXC, bigR, vecX, matS );
	altOmegaSurfMesh(i1,i2) = funchOmega(calcMinfordCurve__evalXSurf( vecXC, bigR, vecX, matS ));
	%
	vecD = vecX - vecXC;
	normD = norm(vecD);
	dMesh(i1,i2) = normD;
	if ( 0.0 == normD )
		dSurfMesh(i1,i2) = 0.0;
		distFromSurfMesh(i1,i2) = norm(vecXC);
		altDistFromSurfMesh(i1,i2) = norm(vecXC);
	else
		%%%dSurfMesh(i1,i2) = norm(bigR*matS*vecD/normD);
		%distFromSurfMesh(i1,i2) = norm( vecXC + bigR*matS*vecD/normD - vecX );
		%altDistFromSurfMesh(i1,i2) = norm( vecXC + bigR*vecD/norm(matS\vecD) - vecX );
		dSurfMesh(i1,i2) = norm(bigR*vecD/norm(matS*vecD));
		distFromSurfMesh(i1,i2) = norm( calcMinfordCurve__evalXSurf(vecXC,bigR,vecX,matS) - vecX );
		altDistFromSurfMesh(i1,i2) = norm( vecXC + bigR*(matS*vecD)/norm(vecD) - vecX );
	end
	%%%if ( norm(matS*vecD)^2 > norm(bigR*vecD) )
	if ( norm(vecD)*norm(matS*vecD) > norm(bigR*vecD) )
		omegaBallMesh(i1,i2) = omegaSurfMesh(i1,i2);
		dBallMesh(i1,i2) = dSurfMesh(i1,i2);
	else
		omegaBallMesh(i1,i2) = omegaMesh(i1,i2);
		dBallMesh(i1,i2) = dMesh(i1,i2);
	end
end
end
%
thetaVals = linspace( 0.0, 2.0*pi, numThetaVals );
for n=1:numThetaVals
	theta = thetaVals(n);
	vecX = vecX0 + vecX1Hat*cos(theta) + vecX2Hat*sin(theta);
	vecD = vecX - vecXC;
	normD = norm(vecD);
	if ( normD == 0.0 )
		vecXSurf = vecXC;
	else
		%%%vecXSurf = vecXC + (bigR*matS*vecD)/normD;
		vecXSurf = vecXC + (bigR*vecD)/norm(matS*vecD);
	end
	x1OfThetaVals(n) = vecX1Hat'*(vecXSurf-vecX0);
	x2OfThetaVals(n) = vecX2Hat'*(vecXSurf-vecX0);
	omegaOfThetaVals(n) = funchOmega( vecXSurf );
end
%
numFigs++; figure(numFigs);
plot( thetaVals, omegaOfThetaVals );
xlabel( "theta" );
ylabel( "omega" );
title( "omega vs theta" );
grid on;
%
numFigs++; figure(numFigs);
plot( thetaVals, x1OfThetaVals, 'o-', thetaVals, x2OfThetaVals, 'x-' );
xlabel( "theta" );
ylabel( "x1, x2" );
title( "x1, x2 vs theta" );
grid on;
%
%funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
funchVizLog = @(x)( sqrt(x) );
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, funchVizLog(omegaMesh), contourPlot_numLevels );
colormap( mycmap(contourPlot_numColors) );
hold on;
plot( x1OfThetaVals, x2OfThetaVals, 'r-' );
hold off;
xlabel( "x1" );
ylabel( "x2" );
title( "(omega) vs x1, x2" );
grid on;
omegaMeshScale = [ min(min(omegaMesh)), max(max(omegaMesh)) ]
funchVizLog = @(x)( sqrt(x) );
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, dMesh, contourPlot_numLevels );
colormap( mycmap(contourPlot_numColors) );
hold on;
plot( x1OfThetaVals, x2OfThetaVals, 'r-' );
hold off;
xlabel( "x1" );
ylabel( "x2" );
title( "d vs x1, x2" );
grid on;
dMeshScale = [ min(min(dMesh)), max(max(dMesh)) ]
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, dSurfMesh, contourPlot_numLevels );
colormap( mycmap(contourPlot_numColors) );
hold on;
plot( x1OfThetaVals, x2OfThetaVals, 'r-' );
hold off;
xlabel( "x1" );
ylabel( "x2" );
title( "dSurf vs x1, x2" );
grid on;
dSurfMeshScale = [ min(min(dSurfMesh)), max(max(dSurfMesh)) ]
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, distFromSurfMesh, contourPlot_numLevels );
colormap( mycmap(contourPlot_numColors) );
hold on;
plot( x1OfThetaVals, x2OfThetaVals, 'r-' );
hold off;
xlabel( "x1" );
ylabel( "x2" );
title( "distFromSurf vs x1, x2" );
grid on;
distFromSurfMeshScale = [ min(min(distFromSurfMesh)), max(max(distFromSurfMesh)) ]
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, altDistFromSurfMesh, contourPlot_numLevels );
colormap( mycmap(contourPlot_numColors) );
hold on;
plot( x1OfThetaVals, x2OfThetaVals, 'r-' );
hold off;
xlabel( "x1" );
ylabel( "x2" );
title( "altDistFromSurf vs x1, x2" );
grid on;
altDistFromSurfMeshScale = [ min(min(altDistFromSurfMesh)), max(max(altDistFromSurfMesh)) ]
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, funchVizLog(omegaSurfMesh), contourPlot_numLevels );
colormap( mycmap(contourPlot_numColors) );
hold on;
plot( x1OfThetaVals, x2OfThetaVals, 'r-' );
hold off;
xlabel( "x1" );
ylabel( "x2" );
title( "(omegaSurf) vs x1, x2" );
grid on;
omegaSurfMeshScale = [ min(min(omegaSurfMesh)), max(max(omegaSurfMesh)) ]
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, funchVizLog(altOmegaSurfMesh), contourPlot_numLevels );
colormap( mycmap(contourPlot_numColors) );
hold on;
plot( x1OfThetaVals, x2OfThetaVals, 'r-' );
hold off;
xlabel( "x1" );
ylabel( "x2" );
title( "(altOmegaSurfMesh) vs x1, x2" );
grid on;
altOmegaSurfMeshScale = [ min(min(altOmegaSurfMesh)), max(max(altOmegaSurfMesh)) ]
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, dBallMesh, contourPlot_numLevels );
colormap( mycmap(contourPlot_numColors) );
hold on;
plot( x1OfThetaVals, x2OfThetaVals, 'r-' );
hold off;
xlabel( "x1" );
ylabel( "x2" );
title( "dBall vs x1, x2" );
grid on;
dBallMeshScale = [ min(min(dBallMesh)), max(max(dBallMesh)) ]
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, funchVizLog(omegaBallMesh), contourPlot_numLevels );
colormap( mycmap(contourPlot_numColors) );
hold on;
plot( x1OfThetaVals, x2OfThetaVals, 'r-' );
hold off;
xlabel( "x1" );
ylabel( "x2" );
title( "(omegaBall) vs x1, x2" );
grid on;
omegaBallMeshScale = [ min(min(omegaBallMesh)), max(max(omegaBallMesh)) ]
