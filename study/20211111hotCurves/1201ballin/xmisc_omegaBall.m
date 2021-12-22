x1
%
thisFile = "xmisc_omegaBall";
numFigs = 5;
%
vecXC = [ 2.0; 5.0 ]
bigR = 2.0
ax = [ vecXC(1), vecXC(1), vecXC(2), vecXC(2) ] + 1.5*bigR*[-1,1,-1,1];
%ax = [ 2.05, 2.2, 0.15, 0.3 ];
%ax = [ 2.124, 2.130, 0.218, 0.224 ];
%ax = [ 2.4, 2.6, 0.1, 0.3 ]
xpost; thisFile = "xmisc_omegaBall";
bigC = 1.0
%
r1Mesh = x1Mesh - vecXC(1);
r2Mesh = x2Mesh - vecXC(2);
rNormMesh = sqrt( r1Mesh.^2 + r2Mesh.^2 );
rHat1Mesh = r1Mesh./(eps*max(max(rNormMesh))+rNormMesh);
rHat2Mesh = r2Mesh./(eps*max(max(rNormMesh))+rNormMesh);
y1Mesh = vecXC(1) + bigR*rHat1Mesh;
y2Mesh = vecXC(2) + bigR*rHat2Mesh;
matY = [ reshape(y1Mesh,1,[]); reshape(y2Mesh,1,[]) ];
matFPolar = funchF(matY,testFuncPrm);
fPolar1Mesh = reshape(matFPolar(1,:),contourPlot_numX2Vals,contourPlot_numX1Vals);
fPolar2Mesh = reshape(matFPolar(2,:),contourPlot_numX2Vals,contourPlot_numX1Vals);
omegaPolarMesh = 0.5*( fPolar1Mesh.^2 + fPolar2Mesh.^2 );
%%%betaMesh = 0.5*(rNormMesh-bigR).^2;
omegaRadialMesh = bigC*(rNormMesh.^2-bigR^2).^2;
omegaSphereMesh = omegaPolarMesh + omegaRadialMesh;
omegaBallMesh = omegaSphereMesh + (rNormMesh<bigR).*(omegaMesh-omegaSphereMesh);
%
numCirclePts = 1001;
vecTheta = 2.0*pi*linspace(0,1.0,numCirclePts);
matCirclePts = zeros(2,numCirclePts);
matCirclePts(1,:) = vecXC(1) + bigR*cos(vecTheta);
matCirclePts(2,:) = vecXC(2) + bigR*sin(vecTheta);

sr_showSRPts = true;
if (sr_showSRPts)
	sr_matRHat = zeros(2,numCirclePts);
	sr_matRHat(1,:) = cos(vecTheta);
	sr_matRHat(2,:) = sin(vecTheta);
	sr_matS = [ 3.0, 0.0; 0.0, 1.0 ]
	sr_matSRPts = zeros(2,numCirclePts);
	for n=1:numCirclePts
		sr_matSRPts(:,n) = vecXC + sr_matS*sr_matRHat(:,n);
	end
end
%
numFigs++; figure(numFigs);
funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
contourf( x1Mesh, x2Mesh, funchVizLog(omegaMesh), contourPlot_numLevels );
axis equal;
colormap( mycmap(contourPlot_numColors) );
hold on;
plot( vecXC(1), vecXC(2), 'ko', 'linewidth', 3, 'markersize', 10 );
plot( vecXC(1), vecXC(2), 'k+', 'linewidth', 3, 'markersize', 10 );
plot( matCirclePts(1,:), matCirclePts(2,:), 'k-', 'linewidth', 1, 'markersize', 1 );
if (sr_showSRPts)
	plot( sr_matSRPts(1,:), sr_matSRPts(2,:), 'g-', 'linewidth', 2, 'markersize', 1 );
end
grid on;
title( "omegaBall: omega vs x1, x2" );
xlabel( "x1" );
ylabel( "x2" );
hold off;
%
numFigs++; figure(numFigs);
funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
contourf( x1Mesh, x2Mesh, funchVizLog(omegaSphereMesh), contourPlot_numLevels );
axis equal;
colormap( mycmap(contourPlot_numColors) );
hold on;
plot( vecXC(1), vecXC(2), 'ko', 'linewidth', 3, 'markersize', 10 );
plot( vecXC(1), vecXC(2), 'k+', 'linewidth', 3, 'markersize', 10 );
plot( matCirclePts(1,:), matCirclePts(2,:), 'k-', 'linewidth', 1, 'markersize', 1 );
if (sr_showSRPts)
	plot( sr_matSRPts(1,:), sr_matSRPts(2,:), 'g-', 'linewidth', 2, 'markersize', 1 );
end
grid on;
title( "omegaBall: omegaSphere vs x1, x2" );
xlabel( "x1" );
ylabel( "x2" );
hold off;
%
numFigs++; figure(numFigs);
funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
contourf( x1Mesh, x2Mesh, funchVizLog(omegaBallMesh), contourPlot_numLevels );
axis equal;
colormap( mycmap(contourPlot_numColors) );
hold on;
plot( vecXC(1), vecXC(2), 'ko', 'linewidth', 3, 'markersize', 10 );
plot( vecXC(1), vecXC(2), 'k+', 'linewidth', 3, 'markersize', 10 );
plot( matCirclePts(1,:), matCirclePts(2,:), 'k-', 'linewidth', 1, 'markersize', 1 );
if (sr_showSRPts)
	plot( sr_matSRPts(1,:), sr_matSRPts(2,:), 'g-', 'linewidth', 2, 'markersize', 1 );
end
grid on;
title( "omegaBall: omegaBall vs x1, x2" );
xlabel( "x1" );
ylabel( "x2" );
hold off;

msg( thisFile, __LINE__, "Let's do a comparison..." );
index1 = 11
index2 = 20
vecX = [ x1Mesh(index1,index2); x2Mesh(index1,index2) ]

xmisc_omegaBall2
