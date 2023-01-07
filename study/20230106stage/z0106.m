% DRaburn 2023-01-06:
%  This is a mock-up of the PyTorch NN code,
%  to help design before implementing in Python.
clear;
if ( stopsignalpresent() )
	msg( __FILE__, __LINE__, "ERROR: Stop signal already present." );
	return;
endif
setprngstates(0);
sizeX = 1E2
secret_sizeL = min([ sizeX, round(0.1*sqrt(sizeX)) ])
%secret_sizeL = min([ sizeX, round(sqrt(sqrt(1E6*sizeX))) ])
%secret_cVals = [ 1.0, 0.0, 0.0, 0.0 ] % Trivial.
secret_cVals = [ 1.0, 1.0e-2, 1.0e-2, 1.0e-2 ] % Easy?
%secret_cVals = [ 0.0, 1.0, 1.0e-2, 1.0e-2 ] % Moderate?
%secret_noisePrm = [ 0.0, 0.0; 0.0, 0.0; 0.0, 0.0 ] % Trivial
secret_noisePrm = [ 1.0e-12, 1.0e-2; 1.0e-2, 1.0e-2; 1.0e-2, 1.0e-2 ] % Moderate?
%
%
tic();
msgnnl( __FILE__, __LINE__, "Generating function... " );
secret_vecXCrit = randn(sizeX,1);
secret_fCrit = 1.0;
secret_matAS = ...
   secret_cVals(1)*sparse(eye(sizeX,sizeX)) ...
 + secret_cVals(2)*sparse(diag(randn(sizeX,1))) ...
 + secret_cVals(3)*sprandn(sizeX,sizeX,secret_sizeL*1.0/sizeX);
secret_matAW = secret_cVals(4)*randn(secret_sizeL,sizeX);
funchFG = @(x) funcQuad1230( x, secret_vecXCrit, secret_fCrit, secret_matAS, secret_matAW, secret_noisePrm );
secret_funchFG_noiseless = @(x) funcQuad1230( x, secret_vecXCrit, secret_fCrit, secret_matAS, secret_matAW, zeros(3,2) );
vecX0 = zeros(sizeX,1);
printf(sprintf( "(Mem: %0.1e vs %0.1e)... ", sizeof(secret_matAS) + sizeof(secret_matAW), (sizeX^2)*sizeof(secret_fCrit) ));
toc();
%
tic();
msgnnl( __FILE__, __LINE__, "Analyzing function... " );
numNLS = 100; % "Noise Level Samples"
[ rvecFNLS, matGNLS ] = funchFG( vecX0 + zeros(sizeX,numNLS) );
fAvg = sum(rvecFNLS)/numNLS;
fSqAvg = sum(rvecFNLS.^2)/numNLS;
fVar = sqrt(max([ fSqAvg - fAvg^2, 0.0 ]));
gAvg = sum(sqrt(sum(matGNLS.^2,1)))/numNLS;
gSqAvg = sum(sum(matGNLS.^2,1))/numNLS;
gVar = sqrt(max([ gSqAvg - gAvg^2, 0.0 ]));
toc();
%
msg( __FILE__, __LINE__, "Starting solver... " );
startTime = time();
prm = [];
prm.learningRate = 0.1;
prm.momentumFactor = 0.9;
prm.numFevalPerSuperPt = 100;
prm.xTol = 10.0*eps*norm(vecX0) + 10.0*eps*fAvg/sqrt(gSqAvg);
prm.fTol = (eps^0.5)*fVar + 10.0*eps*fAvg;
prm.gTol = (eps^0.5)*gVar + 10.0*eps*gAvg;
prm.fevalLimit = -1;
prm.iterLimit = -1;
prm.timeLimit = 600.0;
prm.stopSignalCheckInterval = 3.0;
prm.progressReportInterval = 1.0;
%prm.solverType = "qnj primordial";
prm.solverType = "sgd";
%
vecX = vecX0;
vecP = zeros(size(vecX));
%
iterCount = 0;
fevalCount = 0;
running_fevalCount = 0;
running_fTot = 0.0;
running_xtgTot = 0.0;
running_vecGTot = zeros(sizeX,1);
running_vecXTot = zeros(sizeX,1);
superPt_vecXPrev = vecX0; % Reasonable, but not entirely self-consistent.
proglog_lastTime = time();
stopsig_lastTime = time();
doMainLoop = true;
%
while (doMainLoop)
	%
	[ f, vecG ] = funchFG( vecX );
	fevalCount++;
	if ( prm.fevalLimit > 0 && fevalCount >= prm.fevalLimit )
		msg( __FILE__, __LINE__, "IMPOSED STOP: fevalCount >= prm.fevalLimit." );
		break;
	elseif ( prm.timeLimit >= 0.0 && time() - startTime >= prm.timeLimit )
		msg( __FILE__, __LINE__, "IMPOSED STOP: time() - startTime >= prm.timeLimit." );
		break;
	elseif ( prm.stopSignalCheckInterval >= 0.0 && time() - stopsig_lastTime >= prm.stopSignalCheckInterval )
		if ( stopsignalpresent() )
			msg( __FILE__, __LINE__, "IMPOSED STOP: Received stop signal." );
			break;
		endif
		stopsig_lastTime = time();
	endif
	%
	running_fevalCount++;
	running_fTot += f;
	running_xtgTot += vecX'*vecG;
	running_vecGTot += vecG;
	running_vecXTot += vecX;
	%
	vecP = ( prm.momentumFactor * vecP ) - ( prm.learningRate * vecG );
	vecX += vecP;
	%
	if ( running_fevalCount < prm.numFevalPerSuperPt )
		continue
	endif
	%
	% Okay, we have our latest superpoint.
	iterCount++;
	%
	superPt_vecX = running_vecXTot / running_fevalCount;
	superPt_vecG = running_vecGTot / running_fevalCount;
	temp_xtgAvg = running_xtgTot / running_fevalCount;
	temp_fAvg = running_fTot / running_fevalCount;
	superPt_f = temp_fAvg - (( temp_xtgAvg - (superPt_vecX'*superPt_vecG) )/2.0);
	%
	if ( norm( superPt_vecG ) <= prm.gTol )
		msg( __FILE__, __LINE__, "SUCCESS: norm( superPt_vecG ) <= prm.gTol." );
		break;
	elseif ( superPt_f <= prm.fTol )
		msg( __FILE__, __LINE__, "SUCCESS: superPt_f <= prm.fTol." );
		break;
	elseif ( norm( vecX - superPt_vecXPrev ) <= prm.xTol )
		% Note: This is NOT based on the super point!
		msg( __FILE__, __LINE__, "IMPOSED STOP: norm( vecX - prev_vecX ) <= prm.xTol." );
		break;
	elseif ( prm.iterLimit >= 0 && iterCount >= prm.iterLimit )
		msg( __FILE__, __LINE__, "IMPOSED STOP: iterCount >= prm.iterLimit." );
		break;
	endif
	%
	if ( time() > proglog_lastTime + prm.progressReportInterval )
		msg( __FILE__, __LINE__, sprintf( "   %10.3e (/%0.3e), %5d (/%d), %8d (/%d):  %10.3e (/%0.3e),  %10.3e (/%0.3e),  %10.3e (/%0.3e)", ...
		  time() - startTime, ...
		  prm.timeLimit, ...
		  iterCount,
		  prm.iterLimit, ...
		  fevalCount, ...
		  prm.fevalLimit, ...
		  norm( vecX - superPt_vecXPrev ), ...
		  prm.xTol, ...
		  temp_fAvg, ...
		  prm.fTol, ...
		  norm(superPt_vecG), ...
		  prm.gTol ) );
	 	proglog_lastTime = time();
	endif
	%
	switch (tolower(prm.solverType))
	case { "sgd" }
		% Nothing to do.
	case {"qnj primordial"}
		error( "prm.solverType qnj primordial is not yet supported." );
	otherwise
		echo__prm_solverType = prm.solverType
		error( "Invalid value of prm.solverType." );
	endswitch
	%
	% Prepare for next iteration.
	running_fevalCount = 0;
	running_fTot = 0.0;
	running_xtgTot = 0.0;
	running_vecGTot = zeros(sizeX,1);
	running_vecXTot = zeros(sizeX,1);
	superPt_vecXPrev = superPt_vecX;
endwhile
vecXF = vecX;
%
msg( __FILE__, __LINE__, sprintf( "   %10.3e (/%0.3e), %5d (/%d), %8d (/%d):  %10.3e (/%0.3e),  %10.3e (/%0.3e),  %10.3e (/%0.3e)", ...
  time() - startTime, ...
  prm.timeLimit, ...
  iterCount,
  prm.iterLimit, ...
  fevalCount, ...
  prm.fevalLimit, ...
  norm( vecX - superPt_vecXPrev ), ...
  prm.xTol, ...
  temp_fAvg, ...
  prm.fTol, ...
  norm(superPt_vecG), ...
  prm.gTol ) );
%
msg( __FILE__, __LINE__, "Results... " );
vecXC = secret_vecXCrit;
[ fN0, vecGN0 ] = secret_funchFG_noiseless( vecX0 );
[ fNF, vecGNF ] = secret_funchFG_noiseless( vecXF );
[ fNC, vecGNC ] = secret_funchFG_noiseless( vecXC );
msg( __FILE__, __LINE__, sprintf( " resX: %0.3e -> %0.3e (vs %0.3e)", norm( vecX0 - vecXC ), norm( vecXF - vecXC ), prm.xTol ) );
msg( __FILE__, __LINE__, sprintf( " resF: %0.3e -> %0.3e (vs %0.3e)", fN0 - fNC, fNF - fNC, prm.fTol ) );
msg( __FILE__, __LINE__, sprintf( " resG: %0.3e -> %0.3e (vs %0.3e)", norm( vecGN0 - vecGNC ), norm( vecGNF - vecGNC ), prm.gTol ) );
msg( __FILE__, __LINE__, sprintf( " elapsed time: %0.3es (/%0.3es)", time() - startTime, prm.timeLimit ) );
msg( __FILE__, __LINE__, sprintf( " fevals: %d (/%d)", fevalCount, prm.fevalLimit ) );
msg( __FILE__, __LINE__, sprintf( " iterations: %d (/%d)", iterCount, prm.iterLimit ) );
