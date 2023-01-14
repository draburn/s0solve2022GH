msgnnl( __FILE__, __LINE__, "" );
printf( "  %8.2E, %5d, %8d:", ...
  time() - startTime, ...
  iterCount,
  fevalCount );
printf( "  %8.2E +/- %8.2E;", ...
  norm( superPt_vecX - vecX0 ), ...
  norm( vecXHarvest - vecXSeed ) );
printf( "  %8.2E +/- %8.2E (%9.2E offset);", ...
  superPt_f, ...
  real(sqrt(superPt_fSqVar)), ...
  superPt_fAvg - superPt_f );
printf( "  %8.2E +/- %8.2E (/%0.2E)", ...
  norm( superPt_vecG ), ...
  real(sqrt(sum( superPt_vecGSqVar ))), ...
  prm.gTol );
printf( "\n");
return;

msg( __FILE__, __LINE__, sprintf( "   %10.3e (/%0.3e); %5d (/%d); %8d (/%d):  %10.3e (/%0.3e);  %10.3e, %10.3e (/%0.3e);  %10.3e (/%0.3e)", ...
  time() - startTime, ...
  prm.timeLimit, ...
  iterCount,
  prm.iterLimit, ...
  fevalCount, ...
  prm.fevalLimit, ...
  norm( vecXHarvest - vecXSeed ), ...
  prm.xTol, ...
  superPt_f, ...
  -1.0, %superPt_fPrev - superPt_f, ...
  prm.fTol, ...
  norm(superPt_vecG), ...
  prm.gTol ) );