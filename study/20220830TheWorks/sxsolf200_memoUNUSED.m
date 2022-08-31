% ROOT FUNCTION
[ vecX, vecF, retCode, fevalCount, stepsCount, datOut ] = sxsolf200( funchF, vecX_initial, vecF_initial=[], prmIn=[] );

% MAJOR SEQUENTIAL FUNCTIONS
[ retCode, fevalIncr, vecF_initial, prm ] = __initPrm( funchF, vecX_initial, vecF_initial, prmIn );
[ retCode, fevalIncr, fModelDat ] = __initFModel( funchF, vecX, vecF, prm );
[ retCode, fevalIncr, fModelDat, vecX_next, vecF_next ] = __takeAction( funchF, fModelDat, studyDat, prm );

% "ACTION" FUNCTIONS
[ retCode, fevalIncr, fModelDat ] = __expandSubspace( vecV, funchF, fModelDat, prm );
[ retCode, fevalIncr, fModelDat ] = __reevalDirection( vecV, funchF, fModelDat, prm );
[ retCode, fevalIncr, fModelDat, vecX_next, vecF_next ] = __tryStep( vecY, funchF, fModelDat, studyDat, prm );
moveToCoast

% RANDOM ACCESS FUNCTIONS
vecU = __applyPrecon( vecRhoF, prm, vecX, vecF );
vecV = __calcOrthonorm( vecU, matV, prm );
vecW = __calcJV( vecV, funchF, vecX, vecF, prm );
vecY = __findCandStep( matH, vecG, matC, matB, bTrgt, prm );
vecY = __findCandStep_pow( vecG, matH, matC, matB, bTrgt, prm );
%%%vecY = __findCandStep_grad vecG, matH, matC, matB, bTrgt, prm );
[ retCode, studyDat ] = __studyFModel( funchF, prm );
[ retCode, fModelDat ] = __moveTo( vecY, vecF_next, fModelDat, prm );
[ retCode, fModelDat ] = __shrinkTR( vecY, fModelDat, prm );
[ retCode, fModelDat ] = __expandTR( vecY, fModelDat, prm );
[ retCode, fModelDat ] = __modifyB( vecY, fModelDat, prm );
moveToPlain
[ retCode, preconDat ] = __calcPrecons( fModelDat, prm );
[ retCode, candStepDat ] = __calcCandSteps( fModelDat, prm );

% DEV FUNCTIONS.
__validatePrm( prm );
__validateFModelDat( fModelDat, prm );
__validateStudyDat( fModelDat, studyDat, prm );


if ( 0~= retCode )
	msgretcodeif( true, __FILE__, __LINE__, retCode );
	return;
endif
