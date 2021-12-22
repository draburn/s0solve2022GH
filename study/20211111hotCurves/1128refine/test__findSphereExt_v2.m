%funchG = from x1;
vecXC = [ 1.0; 1.0 ];
bigR = 1.0;
%theta0 = 2.0
%theta0 = 2.0
theta0 = pi/2.0;
vecX0 = vecXC+bigR*[cos(theta0);sin(theta0)]
prm = [];
findSphereExt_v2b;
prev__vecX = [
   0.536444155551432
   0.113932294303886 ]
echo__vecX = vecX
%
%
%
% Idea for next iteration...
%vecXC = [ 0.5; 0.5 ];
%bigR = sqrt(2.0)*0.5;
%theta0 = 0.0
%vecX0 = vecXC+bigR*[cos(theta0);sin(theta0)]
%
r1Mesh = x1Mesh - vecXC(1);
r2Mesh = x2Mesh - vecXC(2);
rNormMesh = sqrt( r1Mesh.^2 + r2Mesh.^2 );
rHat1Mesh = r1Mesh./(eps*max(max(rNormMesh))+rNormMesh);
rHat2Mesh = r2Mesh./(eps*max(max(rNormMesh))+rNormMesh);
y1Mesh = vecXC(1) + bigR*rHat1Mesh;
y2Mesh = vecXC(2) + bigR*rHat2Mesh;
matY = [ reshape(y1Mesh,1,[]); reshape(y2Mesh,1,[]) ];
matFAngular = funchF(matY,testFuncPrm);
fAngular1Mesh = reshape(matFAngular(1,:),contourPlot_numX2Vals,contourPlot_numX1Vals);
fAngular2Mesh = reshape(matFAngular(2,:),contourPlot_numX2Vals,contourPlot_numX1Vals);
omegaAngularMesh = 0.5*( fAngular1Mesh.^2 + fAngular2Mesh.^2 );
%%%betaMesh = 0.5*(rNormMesh-bigR).^2;
betaMesh = 0.5*(rNormMesh.^2-bigR^2).^2;
omegaSphereMesh = omegaAngularMesh + betaMesh;
omegaBallMesh = omegaSphereMesh + (rNormMesh<bigR).*(omegaMesh-omegaSphereMesh);
%
numFigs = 5;
%
numFigs++; figure(numFigs);
funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
contourf( x1Mesh, x2Mesh, funchVizLog(omegaSphereMesh), contourPlot_numLevels );
colormap( mycmap(contourPlot_numColors) );
grid on;
xlabel( "x1" );
ylabel( "x2" );
findSphereExt_v2b; % Draw markers on figure.
echo__vecX = vecX
hold off;
%
numFigs++; figure(numFigs);
funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
contourf( x1Mesh, x2Mesh, funchVizLog(omegaBallMesh), contourPlot_numLevels );
colormap( mycmap(contourPlot_numColors) );
grid on;
xlabel( "x1" );
ylabel( "x2" );
findSphereExt_v2b; % Draw markers on figure.
echo__vecX = vecX
hold off;



msg( "", __LINE__, "DOING NEW CALCULATION..." );
%
numFigs++; figure(numFigs);
funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
contourf( x1Mesh, x2Mesh, funchVizLog(omegaSphereMesh), contourPlot_numLevels );
colormap( mycmap(contourPlot_numColors) );
grid on;
xlabel( "x1" );
ylabel( "x2" );
findSphereExt_v3;
echo__vecX = vecX
hold off;
%
numFigs++; figure(numFigs);
funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
contourf( x1Mesh, x2Mesh, funchVizLog(omegaBallMesh), contourPlot_numLevels );
colormap( mycmap(contourPlot_numColors) );
grid on;
xlabel( "x1" );
ylabel( "x2" );
findSphereExt_v3;
echo__vecX = vecX
hold off;
