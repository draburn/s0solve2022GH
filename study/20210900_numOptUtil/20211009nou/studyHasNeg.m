clear;
setprngstates(0);
numFigs = 0;
%
% We'll make one direction barely negative.
sizeX = 1000;
sizeF = sizeX;
nu = eps^0.25; % Negativity coefficient.
%
matA0 = randn(sizeF,sizeX);
matH0 = matA0'*matA0;
echo__rcond_matH0 = rcond(matH0)
%
[ matPsi0, matLambda0 ] = eig(matH0);
vecLambda0 = diag(matLambda0);
lambdaAlgMax0 = max( vecLambda0 );
[ lambdaAlgMin0, nOfAlgMin0 ] = min( vecLambda0 );
vecPsiOfAlgMin0 = matPsi0(:,nOfAlgMin0);
%matH = matH0 - ( lambdaAlgMin0 + nu*lambdaAlgMax0 )*vecPsiOfAlgMin0*vecPsiOfAlgMin0'; % Not exactly symmetric!
matH = matH0 - ( lambdaAlgMin0 + nu*lambdaAlgMax0 )*(vecPsiOfAlgMin0*vecPsiOfAlgMin0');
vecG = randn(sizeX,1);
omega0 = 1.0+abs(randn());
return
%
%
%
[ matPsi, matLambda ] = eig(matH);
vecLambda = diag(matLambda);
muCrit = -min(vecLambda)
assert( muCrit > 0.0 );
%
matI = eye(sizeX,sizeX);
sizeMu = 1000;
vecMu = muCrit+10.0.^linspace(-3,3,sizeMu);
matDelta = NaN + zeros(sizeX,sizeMu);
vecD = NaN + zeros(sizeMu,1);
vecOmega = NaN + zeros(sizeMu,1);
for n=1:sizeMu
	mu = vecMu(n);
	%
	vecDelta = - ( matH + mu*matI ) \ vecG;
	d = norm(vecDelta);
	omega = omega0 + vecDelta'*vecG + 0.5*vecDelta'*matH*vecDelta;
	%
	matDelta(:,n) = vecDelta;
	vecD(n) = d;
	vecOmega(n) = omega;
end
vecOmegaCap = cap( vecOmega, -10.0*omega0, 10.0*omega0 );
vecDCap = cap( vecD, 0.0, 10.0 );
%
%
%
numFigs++; figure(numFigs);
plot( ...
  vecLambda0, 'o-', ...
  vecLambda, 'x-' );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  vecDCap, vecOmegaCap, 'o-' );
grid on;
grid on;
%
numFigs++; figure(numFigs);
semilogx( ...
  vecMu, vecOmegaCap, 'o-' );
grid on;
%
numFigs++; figure(numFigs);
semilogx( ...
  vecMu, vecDCap, 'o-' );
grid on;
%