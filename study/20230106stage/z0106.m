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
secret_sizeL = min([ sizeX, round(0.1*sqrt(sizeX)) ]) % "Small".
%secret_sizeL = min([ sizeX, round(sqrt(sqrt(1E6*sizeX))) ]) % "Large".
%secret_cVals = [ 1.0, 0.0, 0.0, 0.0 ] % Trivial.
%secret_cVals = [ 1.0, 1.0e-2, 1.0e-2, 1.0e-2 ] % Easy?
secret_cVals = [ 0.0, 1.0, 1.0e-2, 1.0e-2 ] % Moderate?
%secret_cVals = [ 0.0, 0.0, 1.0, 1.0 ] % Extra tricksy?
secret_noisePrm = [ 0.0, 0.0; 0.0, 0.0; 0.0, 0.0 ] % Trivial (except tolerances may be unreasonable!)
%secret_noisePrm = [ 1.0e-12, 1.0e-2; 1.0e-2, 1.0e-2; 1.0e-2, 1.0e-2 ] % Moderate?
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
%
prm.numFevalPerSuperPt = 100;
prm.xTol = 10.0*eps*norm(vecX0) + 10.0*eps*fAvg/sqrt(gSqAvg);
prm.fTol = (eps^0.5)*fVar + 10.0*eps*fAvg;
prm.gTol = (eps^0.5)*gVar + 10.0*eps*gAvg;
prm.fevalLimit = -1;
prm.iterLimit = -1;
prm.timeLimit = 600.0;
prm.stopSignalCheckInterval = 1.0;
prm.progressReportInterval = 1.0;
%
prm.useQNJ = true;
prm.qnj_basisDropThresh = sqrt(eps);
prm.qnj_fitType = "simple0106";
%prm.qnj_stepType = "placeholder0106";
prm.qnj_stepType = "basic0108";
prm.qnj_maxNumRecords = 50;
%
vecX = vecX0;
vecP = zeros(size(vecX));
iterCount = 0;
fevalCount = 0;
running_fevalCount = 0;
running_fTot = 0.0;
running_xtgTot = 0.0;
running_vecGTot = zeros(sizeX,1);
running_vecXTot = zeros(sizeX,1);
proglog_lastTime = time();
stopsig_lastTime = time();
%
superPt_vecXPrev = [];
superPt_vecGPrev = [];
superPt_fPrev = [];
superPt_wPrev = [];
% These initial superPt values are reasonable, albeit not entirely self-consistent.
superPt_vecX = vecX0;
[ foo, tempIndex ] = max(sum(matGNLS.^2,1));
superPt_vecG = matGNLS(:,tempIndex);
superPt_f = max(rvecFNLS);
superPt_w = prm.numFevalPerSuperPt;
%
record_matX = [];
record_matG = [];
record_rvecF = [];
record_rvecW = [];
%
doMainLoop = true;
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
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% SUPER-POINT ANALYSIS
	%
	% Okay, we have our latest superpoint.
	% Update "prev".
	superPt_vecXPrev = superPt_vecX;
	superPt_vecGPrev = superPt_vecG; % Not used?
	superPt_fPrev = superPt_f;
	superPt_wPrev = superPt_w; % Not Used?
	%
	% Calc new.
	iterCount++;
	%
	superPt_vecX = running_vecXTot / running_fevalCount;
	superPt_vecG = running_vecGTot / running_fevalCount;
	temp_xtgAvg = running_xtgTot / running_fevalCount;
	temp_fAvg = running_fTot / running_fevalCount;
	superPt_f = temp_fAvg - (( temp_xtgAvg - (superPt_vecX'*superPt_vecG) )/2.0);
	superPt_w = running_fevalCount;
	%
	if ( norm( superPt_vecG ) <= prm.gTol )
		msg( __FILE__, __LINE__, "SUCCESS: norm( superPt_vecG ) <= prm.gTol." );
		break;
	elseif ( superPt_f <= prm.fTol )
		msg( __FILE__, __LINE__, "SUCCESS: superPt_f <= prm.fTol." );
		break;
	elseif ( superPt_fPrev - superPt_f <= prm.fTol )
		msg( __FILE__, __LINE__, "SUCCESS: superPt_fPrev - superPt_f <= prm.fTol." );
		break;
	elseif ( norm( superPt_vecX - superPt_vecXPrev ) <= prm.xTol )
		msg( __FILE__, __LINE__, "IMPOSED STOP: norm( superPt_vecX - superPt_vecXPrev ) <= prm.xTol." );
		break;
	elseif ( prm.iterLimit >= 0 && iterCount >= prm.iterLimit )
		msg( __FILE__, __LINE__, "IMPOSED STOP: iterCount >= prm.iterLimit." );
		break;
	endif
	%
	if ( time() > proglog_lastTime + prm.progressReportInterval )
		msg( __FILE__, __LINE__, sprintf( "   %10.3e (/%0.3e); %5d (/%d); %8d (/%d):  %10.3e (/%0.3e);  %10.3e, %10.3e (/%0.3e);  %10.3e (/%0.3e)", ...
		  time() - startTime, ...
		  prm.timeLimit, ...
		  iterCount,
		  prm.iterLimit, ...
		  fevalCount, ...
		  prm.fevalLimit, ...
		  norm( superPt_vecX - superPt_vecXPrev ), ...
		  prm.xTol, ...
		  superPt_f, ...
		  superPt_fPrev - superPt_f, ...
		  prm.fTol, ...
		  norm(superPt_vecG), ...
		  prm.gTol ) );
	 	proglog_lastTime = time();
	endif
	%
	% Prepare for next iteration.
	running_fevalCount = 0;
	running_fTot = 0.0;
	running_xtgTot = 0.0;
	running_vecGTot = zeros(sizeX,1);
	running_vecXTot = zeros(sizeX,1);
	superPt_vecXPrev = superPt_vecX;
	if ( ~prm.useQNJ )
		continue;
	endif
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% QUASI-NEWTON-JUMP ANALYSIS
	%
	% Add latest super point to record.
	record_matX = [ superPt_vecX, record_matX ];
	record_matG = [ superPt_vecG, record_matG ];
	record_rvecF = [ superPt_f, record_rvecF ];
	record_rvecW = [ superPt_w, record_rvecW ]; % Not used?
	if ( size(record_matX,2) < 2 )
		% Can't make a jump yet.
		% Also, we won't drop any information.
		continue;
	endif
	%
	%
	% Determine anchor and basis.
	[ fAnchor, indexAnchor ] = min(record_rvecF);
	vecXAnchor = record_matX(:,indexAnchor);
	vecGAnchor = record_matG(:,indexAnchor);
	matD = record_matX - vecXAnchor;
	[ matV, rvecDrop ] = utorthdrop( matD, prm.qnj_basisDropThresh );
	% DRaburn 2023-01-07: My orthogonalization code seems faster than Octave's default ortho()!
	sizeK = size(matV,2);
	rvecDrop(indexAnchor) = false; % Never drop the anchor.
	%
	%
	% Calculate intermediate subspace-related stuff that could probably be determined by utorthdrop().
	matDSans = matD(:,~rvecDrop);
	matGSans = record_matG(:,~rvecDrop);
	rvecFSans = record_rvecF(~rvecDrop);
	matY = triu( matV'*matDSans ); % Could readily be returned by orthogonalization code.
	matGamma = matV'*matGSans;
	vecGammaAnchor = matV'*vecGAnchor;
	%
	%
	% Generate fit.
	switch ( tolower(prm.qnj_fitType) )
	case { "simple0106" }
		%  There are tons of alternatives, but this is almost certainly the simplest sensible method.
		fFit = fAnchor;
		vecGammaFit = vecGammaAnchor;
		matA = (matY') \ (( matGamma - vecGammaAnchor)');
		matHFit = (matA'+matA)/2.0;
	otherwise
		echo__prm_qnj_fitType = prm.qnj_fitType
		error( "Invalid value of prm.qnj_fitType." );
	endswitch
	%
	%
	% Generate step.
	switch ( tolower(prm.qnj_stepType) )
	case { "placeholder0106" }
		% This is a simple placeholder:
		%  no use of a trust region nor consideration that vecX is not superPt_vecX.
		% Testing indicates that 20221216sxsolve/hessfit.m works better.
		[ matR, cholFlag ] = chol( matHFit );
		epsChol = mygetfield( prm, "epsChol", sqrt(eps) );
		if ( 0 == cholFlag && min(diag(matR)) > epsChol*max(abs(diag(matR))) )
			newtStepCoeff = 1.0;
			vecDelta = matV * (matR\(  matR'  \  ((-newtStepCoeff)*vecGammaFit)  ));
			vecX += vecDelta;
			% We're not modifying the momentum.
		endif
		% Otherwise, do nothing.
	case { "basic0108" }
		% We'll use record-based TR scaling and external code to find the proper point on the Levenber curve.
		error( "Not implemented!" );
	otherwise
		echo__prm_qnj_stepType = prm.qnj_stepType
		error( "Invalid value of prm.qnj_stepType." );
	endswitch
	assert( isrealarray(vecX,[sizeX,1]) );
	assert( isrealarray(vecP,[sizeX,1]) );
	%
	%
	% Curate record.
	if ( size(record_matX,2) >= prm.qnj_maxNumRecords )
		record_matX = record_matX(:,1:prm.qnj_maxNumRecords);
		record_matG = record_matG(:,1:prm.qnj_maxNumRecords);
		record_rvecF = record_rvecF(1:prm.qnj_maxNumRecords);
		record_rvecW = record_rvecW(1:prm.qnj_maxNumRecords);
	endif
endwhile
vecXF = vecX;
%
msg( __FILE__, __LINE__, sprintf( "   %10.3e (/%0.3e); %5d (/%d); %8d (/%d):  %10.3e (/%0.3e);  %10.3e, %10.3e (/%0.3e);  %10.3e (/%0.3e)", ...
  time() - startTime, ...
  prm.timeLimit, ...
  iterCount,
  prm.iterLimit, ...
  fevalCount, ...
  prm.fevalLimit, ...
  norm( superPt_vecX - superPt_vecXPrev ), ...
  prm.xTol, ...
  superPt_f, ...
  superPt_fPrev - superPt_f, ...
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
