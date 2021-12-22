clear;
thisFile = "x1";
commondefs;
tic();
setprngstates(0); % "Good behavior" but not single ext.
numFigs = 0;
%
funcPrm = testFunc_genPrm();
%
x1Lo = -5.0; x1Hi = +5.0; x2Lo = -5.0; x2Hi = +5.0;
%x1Lo = -2.4; x1Hi = -1.8; x2Lo = -0.2; x2Hi = 0.4;
%
%
%
matI = eye(2,2);
numSVals = 100;
sVals = linspace(0.0,1.0,numSVals);
%
vecXLev0 = [ 4.0; 4.0 ];
[ omegaDemo, vecFDemo, matJDemo, vecGDemo, matHDemo ] = testFunc_evalDeriv( vecXLev0, funcPrm );
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
[ omegaDemo, vecFDemo, matJDemo, vecGDemo, matHDemo ] = testFunc_evalDeriv( vecXLevB0, funcPrm );
vecGLevB0 = vecGDemo;
matHLevB0 = matHDemo;
vecXLevJTJB0 = matXLevJTJ(:,end);
[ omegaDemo, vecFDemo, matJDemo, vecGDemo, matHDemo ] = testFunc_evalDeriv( vecXLevJTJB0, funcPrm );
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
[ omegaDeom, vecFDemo, matJDemo, vecGDemo, matHDemo ] = testFunc_evalDeriv( vecXLevC0, funcPrm );
vecGLevC0 = vecGDemo;
matHLevC0 = matHDemo;
vecXLevJTJC0 = matXLevJTJB(:,end);
[ omegaDemo, vecFDemo, matJDemo, vecGDemo, matHDemo ] = testFunc_evalDeriv( vecXLevJTJC0, funcPrm );
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
	[ omegaTemp, vecFTemp, matJTemp ] = testFunc_evalDeriv( vecXHOTGrad, funcPrm );
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
[ omega_HOTGC, vecF_HOTGC, matJ_HOTGC, vecG_HOTGC, matH_HOTGC ] = testFunc_evalDeriv( vecX_HOTGC, funcPrm );
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
	[ omega_HOTGC, vecF_HOTGC, matJ_HOTGC, vecG_HOTGC, matH_HOTGC ] = testFunc_evalDeriv( vecX_HOTGC, funcPrm );
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
[ omega0_RK4HOTGC, vecF0_RK4HOTGC, matJ0_RK4HOTGC ] = testFunc_evalDeriv( vecX0_RK4HOTGC, funcPrm );
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
	[ omegaTemp, vecFTemp, matJTemp ] = testFunc_evalDeriv( vecX_RK4HOTGC, funcPrm );
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
matF = testFunc_eval(matX,funcPrm);
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
contourf( x1Mesh, x2Mesh, log(eps025*max(max(omegaMesh))+omegaMesh-min(min(omegaMesh))), 30 );
hold on;
plot( funcPrm.vecXE(1), funcPrm.vecXE(2), 'w+', 'markersize', 20, 'linewidth', 3 );
plot( funcPrm.vecXE(1), funcPrm.vecXE(2), 'ko', 'markersize', 20, 'linewidth', 3 );
plot( matXLevJTJ(1,:), matXLevJTJ(2,:), 'y-', 'linewidth', 3 );
plot( matXLevJTJB(1,:), matXLevJTJB(2,:), 'm-', 'linewidth', 3 );
plot( matXLevJTJC(1,:), matXLevJTJC(2,:), 'c-', 'linewidth', 3 );
plot( matXLev(1,:), matXLev(2,:), 'ro-' );
plot( matXLevB(1,:), matXLevB(2,:), 'gs-' );
plot( matXLevC(1,:), matXLevC(2,:), 'b^-' );
plot( matX_HOTGC(1,:), matX_HOTGC(2,:), 'k*-', 'linewidth', 2 );
grid on;
hold off;
axis equal;
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "log(omega-omegaMin) vs x1, x2" );
