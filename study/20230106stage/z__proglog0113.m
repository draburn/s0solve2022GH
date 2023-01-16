msgnnl( __FILE__, __LINE__, "" );
printf( "%8.2E, %5d (%5dX), %9d:", ...
  time() - startTime, ...
  iterCount, ...
  badCount, ...
  fevalCount );
printf( "  %3d / %3d:", ...
  size(qnj_matV,2), ...
  size(record_matX,2) );
printf( "  %8.2E /%9.2E,  %8.2E /%9.2E, [%8.2E];", ...
  myternary( isempty(qnj_s), 0.0, qnj_s ), ...
  myternary( isempty(qnj_sMax), -1.0, qnj_sMax ), ...
  norm( qnj_vecDelta ), ...
  myternary( isempty(qnj_dMax), -1.0, qnj_dMax ), ...
  norm( vecPSeed ) );
printf( "  %8.2E +/- %8.2E (%8.2E);", ...
  norm( superPt_vecX - vecX0 ), ...
  norm( vecXHarvest - vecXSeed ), ...
  norm( superPt_vecX - minf_vecX ) );
printf( "  %8.2E +/- %8.2E (%8.2E, adj%9.2E);", ...
  superPt_f, ...
  real(sqrt(superPt_fSqVar)), ...
  minf_f, ...
  superPt_fAvg - superPt_f );
printf( "  %8.2E +/- %8.2E (/%0.2E)", ...
  norm( superPt_vecG ), ...
  real(sqrt(sum( superPt_vecGSqVar ))), ...
  prm.gTol );
if ( newIsMinf )
	printf( " ." );
elseif ( newIsBest )
	printf( "  " );
else
	printf( " X" );
endif
printf( "\n");
%
if (0)
	[ secret_f, secret_vecG ] = initDat.funchFG_noiseless(vecXSeed);
	secret_f
	secret_g = norm(secret_vecG)
endif
