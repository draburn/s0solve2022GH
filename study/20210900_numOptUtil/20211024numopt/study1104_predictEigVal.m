clear;
commondefs;
thisFile = "study1104_predictEigVal";
msg( thisFile, __LINE__, "Started." );
setprngstates(0);
numFigs = 0;
tic();
%
sizeX = 100;
sizeFA = sizeX;
sizeFB = sizeX;
%
msgnnl( thisFile, __LINE__, sprintf( ...
  "Generating matrices ( %d, %d, %d )...  ", ...
  sizeX, sizeFA, sizeFB ) );
tic();
matA = randn(sizeFA,sizeX);
matB = randn(sizeFB,sizeX);
matH = (matA'*matA)-(matB'*matB)/(sizeX^2);
%%%matH = (matA'*matA)-(matB'*matB)*0.1;
vecG = randn(sizeX,1);
matI = eye(sizeX,sizeX);
hFrobNorm = sqrt(sum(sum(matH.^2)));
toc();
msg( thisFile, __LINE__, sprintf( "Frob norm: %10.3e.", hFrobNorm ) );
%
msgnnl( thisFile, __LINE__, "Performing eig( matH )...  " );
tic();
[ matPsi, matLam ] = eig( matH ); vecLam = diag(matLam);
toc();
msg( thisFile, __LINE__, sprintf( ...
  "Eigenvaule range: %10.3e ~ %10.3e.", ...
  min(vecLam), max(vecLam) ) );
%
assert( min(vecLam) < 0.0 );
muCrit = -min(vecLam);
assert( max(vecLam) > 0.0 );
muSafe = max(vecLam);
%
numMuVals = 100;
muLo = muCrit*(1.0+eps050) + eps075*muSafe;
muHi = muSafe;
muVals = muLo + (muHi-muLo)*(linspace(0.0,1.0,numMuVals).^4);
%
cholFlagVals = 0*muVals - 1;
pVals = 0.0*muVals - 1.0;
mio_muPredVals = 0.0*muVals - 1.0;
%
alpha_opts.tol = 0.1;
alpha_opts.maxit = 1;
%alpha_opts.p = 20; % I think 20 is the default.
alpha_opts.v0 = vecG;
alpha_opts.cholB = true;
alpha_flagVals = 0*muVals - 1;
alpha_muPredVals = 0.0*muVals - 1.0;
%
beta_opts.tol = 0.1;
beta_opts.maxit = 1;
beta_opts.p = sizeX-1; % I think 20 is the default.
beta_opts.v0 = vecG;
beta_opts.cholB = true;
beta_flagVals = 0*muVals - 1;
beta_muPredVals = 0.0*muVals - 1.0;
%
msgnnl( thisFile, __LINE__, "Performing main loop...  " );
tic();
for n=1:numMuVals
	mu = muVals(n);
	matM = matH + (mu*matI);
	[ matR, cholFlag ] = chol( matM );
	cholFlagVals(n) = cholFlag;
	if ( 0 == cholFlag )
		vecDelta = -( matR \ (matR'\vecG) );
		vecDeltaPrime = -( matR \ (matR'\vecDelta) );
		kappa = vecDelta'*vecDelta;
		kappaPrime = vecDelta'*vecDeltaPrime;
		kappaPrimePrime = 3.0*(vecDeltaPrime'*vecDeltaPrime);
		mio_pVals(n) = (kappaPrime^2) / ( (kappa*kappaPrimePrime) - (kappaPrime^2) );
		mio_muPredVals(n) = mu + (kappa * kappaPrime / ( (kappa*kappaPrimePrime) - (kappaPrime^2) ));
		%
		[ v, d, eigsFlag ] = eigs( matI, matR, 1, 'lm', alpha_opts );
		alpha_flagVals(n) = eigsFlag;
		if ( 0 == eigsFlag )
			alpha_muPredVals(n) = mu - (1.0/d);
		end
		%
		[ v, d, eigsFlag ] = eigs( matI, matR, 1, 'lm', beta_opts );
		beta_flagVals(n) = eigsFlag;
		if ( 0 == eigsFlag )
			beta_muPredVals(n) = mu - (1.0/d);
		end
	end
end
toc();
%
mio_mask = (mio_muPredVals>0.0) & (mio_muPredVals<muVals);
alpha_mask = (alpha_muPredVals>0.0) & (alpha_muPredVals<muVals);
beta_mask = (beta_muPredVals>0.0) & (beta_muPredVals<muVals);
%
msgnnl( thisFile, __LINE__, "Generating plots...  " );
tic();
%
numFigs++; figure(numFigs);
loglog( ...
  muVals, 0.0*muVals + muCrit, 'k-', 'linewidth', 3, ...
  muVals(alpha_mask), alpha_muPredVals(alpha_mask), 'o-', 'markersize', 15, ...
  muVals(beta_mask), beta_muPredVals(beta_mask), 's-', 'markersize', 15, ...
  muVals(mio_mask), mio_muPredVals(mio_mask), 'x-', 'markersize', 15 );
grid on;
title( "muPred vs mu" );
xlabel( "mu" );
ylabel( "muPred" );
%
toc();
%
msg( thisFile, __LINE__, "" );
msg( thisFile, __LINE__, "Conclusion..." );
msg( thisFile, __LINE__, "This illustrates that even a crude eigs() is better than mio;" );
msg( thisFile, __LINE__, "where crude eigs() fails, mio is horribly wrong anyway." );
msg( thisFile, __LINE__, "But, maybe: mio at 2 vals of mu, make a straight line...?" );
