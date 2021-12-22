clear;
thisFile = "studyHasNeg3 ?"
tic();
setprngstates(25848832);
numFigs = 0;
%
msg( thisFile, __LINE__, "" );
msg( thisFile, __LINE__, "THIS IS A CASE WHERE FIRST ITER MERELY GETS US TO A VALID FINITE MU." );
msg( thisFile, __LINE__, "" );
%
% We'll make one direction barely negative.
sizeX = 20;
sizeF = sizeX;
muCritApplied = exp(3.0*randn());
muCritGuess = 0.0;
%muCritGuess = 4.7971
omega0 = exp(3.0*randn());
muLo = 10.0^floor(log10(muCritApplied/10.0));
muHi = 1.0E11; % Adjust as necessary.
%
matA0 = randn(sizeF,sizeX).*exp(3.0*randn(sizeF,sizeX));
matH0 = matA0'*matA0;
echo__rcond_matH0 = rcond(matH0)
%
[ matEigVec0, matEigVal0 ] = eig(matH0);
[ eigValMin0, nOfMin0 ] = min(diag( matEigVal0 ));
eigVecMin0 = matEigVec0(:,nOfMin0);
matH = matH0 - ( eigValMin0 + muCritApplied )*(eigVecMin0*eigVecMin0');
vecG = randn(sizeX,1).*exp(3.0*randn(sizeX,1));
%
matI = eye(sizeX,sizeX);
sizeMu = 1000;
%vecMu = 10.0.^linspace(-3,3,sizeMu)';
vecMu = muLo * ( (muHi/muLo) .^ linspace(0.0,1.0,sizeMu)' );
matDelta = zeros(sizeX,sizeMu);
vecLambda = zeros(sizeMu,1);
vecOmega = omega0+zeros(sizeMu,1);
matDeltaP = zeros(sizeX,sizeMu);
matDeltaPP = zeros(sizeX,sizeMu);
vecDeltaPNorm = zeros(sizeMu,1);
vecDeltaPPNorm = zeros(sizeMu,1);
for n=1:sizeMu
	mu = vecMu(n);
	%
	[ matR, errFlag ] = chol( matH + mu*matI );
	if (~errFlag)
		vecDelta = -( matR \ (matR'\vecG) );
		lambda = norm(vecDelta);
		omega = omega0 + vecDelta'*vecG + 0.5*vecDelta'*matH*vecDelta;
		%
		matDelta(:,n) = vecDelta;
		vecLambda(n) = lambda;
		vecOmega(n) = omega;
		%
		%
		vecDeltaP = -( matR \ (matR'\vecDelta) );
		vecDeltaPP = -2.0*( matR \ (matR'\vecDeltaP) );
		matDeltaP(:,n) = vecDeltaP;
		matDeltaPP(:,n) = vecDeltaPP;
		vecDeltaPNorm(n) = norm(vecDeltaP);
		vecDeltaPPNorm(n) = norm(vecDeltaPP);
		clear vecDeltaP;
		clear vecDeltaPP;
		%
		%
		%
		clear omega;
		clear lambda;
		clear vecDelta;
	end
	%
	clear mu;
end
clear n;
%
%
%
sizeMuModel = 1000;
%%%vecMuModel = 10.0.^linspace(-3,3,sizeMuModel)';
vecMuModel = muLo * ( (muHi/muLo) .^ linspace(0.0,1.0,sizeMuModel)' );
matDeltaModel = zeros(sizeX,sizeMuModel);
vecLambdaModel = zeros(sizeMuModel,1);
vecOmegaModelPart = omega0+zeros(sizeMuModel,1);
vecOmegaModelFull = omega0+zeros(sizeMuModel,1);
vecG2 = (matH + muCritGuess*matI)*vecG;
bigADirect = omega0;
bigBDirect = vecG'*vecG;
vecOmegaDirect = omega0+zeros(sizeMuModel,1);
%
neo_bigA = -vecG'*vecG
neo_bigB = muCritGuess * (vecG'*vecG) + 1.5*(vecG'*matH*vecG)
neo_bigC = muCritGuess
neo_alpha = neo_bigA
neo_beta = neo_bigC + neo_bigB/neo_bigA
neo_vecOmega = zeros(sizeMuModel,1);
%
for n=1:sizeMuModel
	mu = vecMuModel(n);
	nu = mu - muCritGuess;
	%
	if (nu>0)
		vecDelta = -vecG/nu + vecG2/(nu^2);
		lambda = norm(vecDelta);
		omegaFull = omega0 + vecDelta'*vecG + 0.5*vecDelta'*matH*vecDelta;
		omegaPart = omega0 + vecDelta'*vecG + 0.5*vecG'*matH*vecG/(nu^2);
		%
		matDeltaModel(:,n) = vecDelta;
		vecLambdaModel(n) = lambda;
		vecOmegaModelPart(n) = omegaPart;
		vecOmegaModelFull(n) = omegaFull;
		%
		vecOmegaDirect(n) = bigADirect - bigBDirect/nu;
		%
		clear omegaFull;
		clear omegaPart;
		clear lambda;
		clear vecDelta;
	end
	%
	if ( mu>neo_beta )
		neo_vecOmega(n) = omega0 + neo_alpha / ( mu-neo_beta);
	end
	%
	clear nu
	clear mu;
end
clear n;
%
%
%
omegaTarget = 0.25*omega0
a = 2.0*muCritGuess + 3.0*(vecG'*matH*vecG)/(vecG'*vecG)
c = 2.0*( omega0 - omegaTarget )/(vecG'*vecG)
assert( c > 0.0 );
ac = a*c
if ( abs(a*c) < eps^0.75 )
	msg( thisFile, __LINE__, "Using a = 0 root form." );
	x = c/2.0
elseif ( a*c < 1.0 )
	msg( thisFile, __LINE__, "Using a != 0 root form." );
	x = ( 1.0 - sqrt( 1.0 - a*c ) ) / a
else
	msg( thisFile, __LINE__, "Using minimum form." );
	x = 1.0 / a
end
muSuggested = muCritGuess + 1.0/x
[ matR, errFlag ] = chol( matH + muSuggested*matI );
if (errFlag)
	msg( thisFile, __LINE__, "Cholesky decompsition for mu suggested failed." );
	lambdaSuggested = 0.0;
	omegaSuggested = omega0;
else
	vecDeltaSuggested = -( matR \ (matR'\vecG) );
	lambdaSuggested = norm(vecDeltaSuggested);
	omegaSuggested = omega0 + vecDeltaSuggested'*vecG + 0.5*vecDeltaSuggested'*matH*vecDeltaSuggested;
end
%
%
%
muDirect = muCritGuess + bigBDirect/(bigADirect-omegaTarget)
[ matR, errFlag ] = chol( matH + muDirect*matI );
if (errFlag)
	error( "Cholesky decompsition for mu direct failed." );
else
	vecDeltaDirect = -( matR \ (matR'\vecG) );
	lambdaDirect = norm(vecDeltaDirect);
	omegaDirect = omega0 + vecDeltaDirect'*vecG + 0.5*vecDeltaDirect'*matH*vecDeltaDirect;
	vecDeltaDirectP = -( matR \ (matR'\vecDeltaDirect) );
	omegaDirectP = ( vecG + matH*vecDeltaDirect )' * vecDeltaDirectP;
end
%
%
%
vecG2 = (matH + muCritGuess*matI)*vecG;
bigADirect2 = omega0;
bigBDirect2 = -( omega0 - omegaDirect )^2 / omegaDirectP;
bigCDirect2 = muDirect - ( omega0 - omegaDirect ) / omegaDirectP
vecOmegaDirect2 = omega0+zeros(sizeMuModel,1);
%
for n=1:sizeMuModel
	mu = vecMuModel(n);
	nu = mu - bigCDirect2;
	%
	if (nu>0)
		vecOmegaDirect2(n) = bigADirect2 + bigBDirect2/nu;
	end
	%
	clear nu
	clear mu;
end
clear n;
%
%
%
vecOmegaCap = cap( vecOmega, -2.0*omega0, 2.0*omega0 );
vecLambdaCap = cap( vecLambda, 0.0, 10.0 );
vecOmegaModelFullCap = cap( vecOmegaModelFull, -2.0*omega0, 2.0*omega0 );
vecOmegaModelPartCap = cap( vecOmegaModelPart, -2.0*omega0, 2.0*omega0 );
vecLambdaModelCap = cap( vecLambdaModel, 0.0, max(vecLambdaCap) );
vecOmegaDirectCap = cap( vecOmegaDirect, -2.0*omega0, 2.0*omega0 );
vecOmegaDirect2Cap = cap( vecOmegaDirect2, -2.0*omega0, 2.0*omega0 );
neo_vecOmegaCap = cap( neo_vecOmega, -2.0*omega0, omega0 );
%
numFigs++; figure(numFigs);
semilogx( ...
  vecMu, vecOmegaCap, 'o-', ...
  muSuggested, omegaSuggested, '+', 'markersize', 25, 'linewidth', 2, ...
  muSuggested, omegaSuggested, 'o', 'markersize', 25, 'linewidth', 2, ...
  muDirect, omegaDirect, '+', 'markersize', 25, 'linewidth', 2, ...
  muDirect, omegaDirect, 'o', 'markersize', 25, 'linewidth', 2, ...
  vecMuModel, neo_vecOmegaCap, '^-', 'linewidth', 2, ...
  vecMuModel, vecOmegaModelFullCap, '-', 'linewidth', 2, ...
  vecMuModel, vecOmegaModelPartCap, '-', 'linewidth', 2, ...
  vecMuModel, vecOmegaDirectCap, '-', 'linewidth', 2, ...
  vecMuModel, vecOmegaDirect2Cap, '^-', 'linewidth', 2, ...
  vecMu, omegaTarget+0*vecOmegaCap, 'c-', 'linewidth', 2, ...
  muSuggested+0*vecMu, vecOmegaCap, 'r-', ...
  vecMu, omega0+0*vecOmegaCap, 'k-', ...
  vecMu, omegaSuggested+0*vecOmegaCap, 'g-' );
xlabel( "mu" );
ylabel( "omega cap" );
title( "omega cap vs mu" );
grid on;
%
numFigs++; figure(numFigs);
semilogx( ...
  vecMu, vecOmegaDirect2Cap-vecOmegaCap, 'o-' );
xlabel( "mu" );
ylabel( "omega cap diff" );
title( "omega direct 2 cap minus omega cap vs mu" );
grid on;
%
numFigs++; figure(numFigs);
semilogx( ...
  vecMu, vecLambdaCap, 'o-', ...
  muSuggested, lambdaSuggested, '+', 'markersize', 25, 'linewidth', 2, ...
  muSuggested, lambdaSuggested, 'o', 'markersize', 25, 'linewidth', 2, ...
  muDirect, lambdaDirect, '+', 'markersize', 25, 'linewidth', 2, ...
  muDirect, lambdaDirect, 'o', 'markersize', 25, 'linewidth', 2, ...
  vecMuModel, vecLambdaModelCap, '-', 'linewidth', 2 );
xlabel( "mu" );
ylabel( "lambda cap" );
title( "lambda cap vs mu" );
grid on;
%
%
%
toc();
