function [ vecX, retCode, datOut ] = findLocMinLossG_bfgs( funchBigL, funchVecG, vecX0, prm=[] )
	commondefs;
	thisFile = "findLocMinLossG_bfgs";
	valdLev = mygetfield( prm, "valdLev", VALDLEV__MEDIUM );
	verbLev = mygetfield( prm, "verbLev", VERBLEV__WARN );
	assert( isrealscalar(valdLev) );
	assert( isrealscalar(verbLev) );
	msg_copious( verbLev, thisFile, __LINE__, "Welcome." );
	msg( thisFile, __LINE__, "TODO: Implement (slightly) more sophisticated backtracking." );
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
	iterLimit = mygetfield( prm, "iterLimit", 1000 );
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
	magDTol = mygetfield( prm, "magDTol", 1e-12 * norm(vecX0) + 1e-15 );
	btLimit = mygetfield( prm, "btLimit", 20 );
	btCoeff = mygetfield( prm, "btCoeff", 0.1 );
	if ( valdLev >= VALDLEV__LOW )
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
	if ( bigL0 <= bigLTol )
		msg_warn( verbLev, thisFile, __LINE__, "Input satisfied bigLTol." );
		retCode = RETCODE__SUCCESS;
		return;
	end
	if ( norm(vecG0) <= magGTol )
		msg_warn( verbLev, thisFile, __LINE__, "Input satisfied gNormTol." );
		retCode = RETCODE__SUCCESS;
		return;
	end
	%
	matA = eye(sizeX,sizeX)*((vecG0'*vecG0)/bigL0);
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
		% Pick initial step.
		[ matR, cholFlag ] = chol( matA );
		if ( 0 ~= cholFlag )
			msg_error( verbLev, thisFile, __LINE__, "Algorithm breakdown: Cholesky factorization." );
			retCode = RETCODE__ALGORITHM_BREAKDOWN;
			return;
		end
		vecD = -( matR \ (matR'\vecG) );
		%
		%
		% Do simple backtracking.
		btCount = 0;
		while (1)
			vecX_next = vecX + vecD;
			bigL_next = funchBigL( vecX_next );
			if ( valdLev >= VALDLEV__HIGH )
				assert( isrealscalar(bigL) );
				assert( 0.0 <= bigL );
			end
			if ( bigL_next < bigL )
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
		%
		%
		% Update approximate Hessian.
		vecDeltaG = vecG - vecG_prev;
		vecDeltaX = vecX - vecX_prev;
		vecU = vecDeltaG;
		vecV = matA * vecDeltaX;
		cuInv = vecU'*vecDeltaX;
		cvInv = vecV'*vecDeltaX;
		if ( abs(cuInv) <= eps050*norm(vecU)*norm(vecDeltaX) ...
		  || abs(cvInv) <= eps050*norm(vecV)*norm(vecDeltaX) );
			msg_error( verbLev, thisFile, __LINE__, "Algorithm breakdown: Updating approximate Hessian." );
			retCode = RETCODE__ALGORITHM_BREAKDOWN;
			return;
		end
		cu = 1.0./cvInv;
		cv = 1.0./cvInv;
		foo = matA + (cu*vecU)*(vecU') - (cv*vecV)*(vecV');
		matA = ( foo' + foo ) / 2.0;
		clear foo;
	end
end

%!test
%!	clear;
%!	commondefs;
%!	thisFile = "findLocMinLossG_bfgs test 1";
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
%!	msg( thisFile, __LINE__, "Calling findLocMinLossG_bfgs()..." );
%!	time0 = time();
%!	[ vecXF, retCode, datOut ] = findLocMinLossG_bfgs( funchBigL, funchVecG, vecX0 );
%!	timeF = time();
%!	%
%!	bigL0 = funchBigL(vecX0);
%!	bigLF = funchBigL(vecXF);
%!	magD0 = norm(vecX0-vecXRoot);
%!	magDF = norm(vecXF-vecXRoot);
%!	msg( thisFile, __LINE__, sprintf( "findLocMinLossG_bfgs() returned retCode %s.", retcode2str(retCode) ) );
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
