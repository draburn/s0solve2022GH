function [ vecX, retCode, datOut ] = findLocMinLossG_simple( funchBigL, funchVecG, vecX0, prm=[] )
	commondefs;
	thisFile = "findLocMinLossG_simple";
	valdLev = mygetfield( prm, "valdLev", VALDLEV__MEDIUM );
	verbLev = mygetfield( prm, "verbLev", VERBLEV__WARN );
	%valdLev = mygetfield( prm, "valdLev", VALDLEV__HIGH );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	assert( isrealscalar(valdLev) );
	assert( isrealscalar(verbLev) );
	msg_copious( verbLev, thisFile, __LINE__, "Welcome." );
	%
	sizeX = size(vecX0,1);
	if ( valdLev >= VALDLEV__LOW )
		assert( isrealarray(vecX0,[sizeX,1]) );
	end
	%
	bigL0 = funchBigL(vecX0);
	vecG0 = funchVecG(vecX0);
	if ( valdLev >= VALDLEV__LOW )
		assert( isrealscalar(bigL0) );
		assert( 0.0 <= bigL0 );
		assert( isrealarray(vecG0,[sizeX,1]) );
	end
	%
	iterLimit = mygetfield( prm, "iterLimit", 10000 );
	bigLTol = mygetfield( prm, "bigLTol", 1e-6 * bigL0 + 1e-9  );
	magGTol = mygetfield( prm, "magGTol", 1e-9 * norm(vecG0) + 1e-12 * sizeX );
	if ( valdLev >= VALDLEV__LOW )
		assert( isrealscalar(iterLimit) )
		assert( isrealscalar(bigLTol) );
		assert( isrealscalar(magGTol) );
		assert( 0 <= iterLimit )
		assert( 0.0 <= bigLTol );
		assert( 0.0 <= magGTol );
	end
	%
	armijoCoeff = mygetfield( prm, "armijoCoeff", 1e-4 );
	magDTol = mygetfield( prm, "magDTol", 1e-12 * norm(vecX0) + 1e-15 );
	btLimit = mygetfield( prm, "btLimit", 20 );
	btCoeff = mygetfield( prm, "btCoeff", 0.1 );
	if ( valdLev >= VALDLEV__LOW )
		assert( 0.0 < armijoCoeff );
		assert( armijoCoeff < 1.0 );
		assert( 0.0 <= magDTol );
		assert( isrealscalar(btLimit) )
		assert( isrealscalar(btCoeff) );
		assert( 0 <= btLimit )
		assert( 0.0 < btCoeff );
		assert( btCoeff < 1.0 );
	end
	%
	bigLAbsFallTol = mygetfield( prm, "bigLAbsFallTol", 1e-12 * bigL0 + 1e-15 );
	bigLRelFallTol = mygetfield( prm, "bigLRelFallTol", 1e-6 );
	if ( valdLev >= VALDLEV__LOW )
		assert( isrealscalar(bigLAbsFallTol) );
		assert( isrealscalar(bigLRelFallTol) );
		assert( 0.0 <= bigLAbsFallTol );
		assert( 0.0 <= bigLRelFallTol );
	end
	%
	genDatOut = mygetfield( prm, "genDatOut", true );
	if (genDatOut)
		datOut.vecXVals = zeros(sizeX,iterLimit+1);
		datOut.bigLVals = zeros(1,iterLimit+1);
		datOut.vecGVals = zeros(sizeX,iterLimit+1);
	else
		datOut = [];
	end
	%
	%
	%
	vecX = vecX0;
	bigL = bigL0;
	vecG = vecG0;
	iterCount = 0;
	if (genDatOut)
		datOut.iterCount = iterCount;
		datOut.vecXVals(:,iterCount+1) = vecX;
		datOut.bigLVals(iterCount+1) = bigL;
		datOut.vecGVals(:,iterCount+1) = vecG;
	end
	while (1)
		% Check pre-iter stop crit.
		if ( bigL <= bigLTol )
			msg_main( verbLev, thisFile, __LINE__, "Success: bigLTol." );
			retCode = RETCODE__SUCCESS;
			return;
		end
		if ( norm(vecG) <= magGTol )
			msg_main( verbLev, thisFile, __LINE__, "Success: gNormTol." );
			retCode = RETCODE__SUCCESS;
			return;
		end
		%
		iterCount++;
		if ( iterCount > iterLimit )
			msg_warn( verbLev, thisFile, __LINE__, "Imposed stop: iterLimit." );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		end
		%
		%
		% Pick initial step;
		vecD = -bigL*vecG/(vecG'*vecG);
		%
		%
		% Do simple backtracking.
		btCount = 0;
		while (1)
			vecX_next = vecX + vecD;
			bigL_next = funchBigL( vecX_next );
			armijoDeltaL = armijoCoeff * (vecD' * vecG);
			if ( valdLev >= VALDLEV__HIGH )
				assert( isrealscalar(bigL) );
				assert( 0.0 <= bigL );
				assert( armijoDeltaL < 0.0 );
			end
			if ( bigL_next < bigL + armijoDeltaL )
				break;
			end
			if ( norm(vecX_next-vecX) <= magDTol )
				msg_warn( verbLev, thisFile, __LINE__, "Algorithm breakdown: magDTol." );
				retCode = RETCODE__ALGORITHM_BREAKDOWN;
				return;
			end
			%
			btCount++;
			if ( btCount > btLimit )
				msg_warn( verbLev, thisFile, __LINE__, "Algorithm breakdown: btLimit." );
				retCode = RETCODE__ALGORITHM_BREAKDOWN;
				return;
			end
			vecD *= btCoeff;
		end
		vecG_next = funchVecG(vecX_next);
		if ( valdLev >= VALDLEV__HIGH )
			assert( bigL_next < bigL );
			assert( isrealarray(vecG_next,[sizeX,1]) );
		end
		%
		%
		% Take step. 
		vecX_prev = vecX;
		bigL_prev = bigL;
		vecG_prev = vecG;
		vecX = vecX_next;
		bigL = bigL_next;
		vecG = vecG_next;
		if (genDatOut)
			datOut.iterCount = iterCount;
			datOut.vecXVals(:,iterCount+1) = vecX;
			datOut.bigLVals(iterCount+1) = bigL;
			datOut.vecGVals(:,iterCount+1) = vecG;
		end
		%
		%
		% Check post-iter stop crit.
		if ( abs(bigL-bigL_prev) <= bigLAbsFallTol )
			msg_notify( verbLev, thisFile, __LINE__, "Imposed stop: bigLAbsFallTol." );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		end
		if ( abs(bigL-bigL_prev) <= bigL*bigLRelFallTol )
			msg_notify( verbLev, thisFile, __LINE__, "Imposed stop: bigLRelFallTol." );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		end
	end
end

%!test
%!	clear;
%!	commondefs;
%!	thisFile = "findLocMinLossG_simple test 1";
%!	setprngstates(0);
%!	numFigs = 0;
%!	%
%!	sizeX = 10;
%!	sizeF = sizeX;
%!	vecXRoot = randn(sizeX,1);
%!	matA = randn(sizeF,sizeX);
%!	funchBigL = @(x)( 0.5*(norm(matA*(x-vecXRoot))^2) );
%!	funchVecG = @(x)( matA'*matA*(x-vecXRoot) );
%!	vecX0 = randn(sizeX,1);
%!	%
%!	msg( thisFile, __LINE__, "Calling findLocMinLossG_simple()..." );
%!	time0 = time();
%!	[ vecXF, retCode, datOut ] = findLocMinLossG_simple( funchBigL, funchVecG, vecX0 );
%!	timeF = time();
%!	%
%!	bigL0 = funchBigL(vecX0);
%!	bigLF = funchBigL(vecXF);
%!	magD0 = norm(vecX0-vecXRoot);
%!	magDF = norm(vecXF-vecXRoot);
%!	msg( thisFile, __LINE__, sprintf( "findLocMinLossG_simple() returned retCode %s.", retcode2str(retCode) ) );
%!	msg( thisFile, __LINE__, sprintf( "  Performed %d iterations", datOut.iterCount ) );
%!	msg( thisFile, __LINE__, sprintf( "  Loss went from %0.3e to %0.3e.", bigL0, bigLF ) );
%!	msg( thisFile, __LINE__, sprintf( "  ||x-xRoot|| went from %0.3e to %0.3e.",magD0, magDF ) );
%!	msg( thisFile, __LINE__, sprintf( "  Elapsed time was %gs.", time()-time0 ) );
%!	%
%!	numFigs++; figure(numFigs);
%!	semilogy( (0:datOut.iterCount), datOut.bigLVals(1:datOut.iterCount+1), 'o-' );
%!	grid on;
%!	xlabel( "iter" );
%!	ylabel( "Loss" );
%!	title( sprintf( "%s: Loss vs iter", thisFile ), "interpreter", "none" );
%!	%
%!	assert( bigLF < 1e-6*bigL0 + 1e-9 )
