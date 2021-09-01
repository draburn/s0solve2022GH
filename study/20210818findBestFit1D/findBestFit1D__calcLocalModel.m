function [ omega, vecG, matH, retCode, datOut ] = findBestFit1D__calcLocalModel( funchRho, rhoArgs, vecZ, prm )
	%
	% Init
	commondefs;
	thisFile = "findBestFit1D__calcLocalModel";
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__NOTIFY );
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	valLev = mygetfield( prm, "valLev", VALLEV__HIGH );
	%
	omega = [];
	vecG = [];
	matH = [];
	retCode = RETCODE__NOT_SET;
	datOut = [];
	%
	%
	%
	% Calc stuff *at* vecZ.
	sizeZ = size(vecZ,1);
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealarray(vecZ,[sizeZ,1]) );
	end
	%
	[ errFlag, vecRho ] = funchRho( rhoArgs, vecZ );
	if (errFlag)
		msg_error( verbLev, thisFile, __LINE__, "funchRho() failed for current guess." );
		retCode = RETCODE__BAD_INPUT;
		return;
	end
	sizeRho = size(vecRho,1);
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealarray(vecRho,[sizeRho,1]) );
	end
	%
	if ( mygetfield(prm,"useCustomOmega",false) )
		funchOmega      = mygetfield(prm,"funchOmega");
		funchOmegaDRho  = mygetfield(prm,"funchOmegaDRho");
		funchOmegaDRho2 = mygetfield(prm,"funchOmegaDRho2");
	else
		funchOmega      = @(rho)( 0.5 * sum(rho.^2) );
		funchOmegaDRho  = @(rho)( rho );
		funchOmegaDRho2 = @(rho)( eye(sizeRho,sizeRho) );
	end
	omega         = funchOmega(      vecRho );
	vecOmegaDRho  = funchOmegaDRho(  vecRho );
	matOmegaDRho2 = funchOmegaDRho2( vecRho );
	if (~isrealscalar(omega))
		msg_error( verbLev, thisFile, __LINE__, "funchOmega() failed for current guess (not a real scalar)." );
		retCode = RETCODE__BAD_INPUT;
		return;
	end
	if ( omega < 0.0 )
		msg_error( verbLev, thisFile, __LINE__, "funchOmega() failed for current guess (negative)." );
		retCode = RETCODE__BAD_INPUT;
		return;
	end
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealscalar(omega) );
		assert( omega >= 0.0 );
		assert( isrealarray(vecOmegaDRho,[sizeRho,1]) );
		assert( isrealarray(matOmegaDRho2,[sizeRho,sizeRho]) );
		assert( sum(sum((matOmegaDRho2'-matOmegaDRho2).^2)) <= eps150*sum(sum(matOmegaDRho2.^2)) );
	end
	if ( issize(vecOmegaDRho,[1,sizeRho]) )
		vecOmegaDRho = vecOmegaDRho';
	end
	%
	%
	%
	% Calc stuff *near* vecZ.
	vecEpsZ = mygetfield( prm, "vecEpsZ", eps050*ones(size(vecZ)) );
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealarray(vecEpsZ,[sizeZ,1]) );
		for n=1:sizeZ
			assert( vecEpsZ(n) > eps*abs(vecZ(n)) );
		end
	end
	%
	matRho_plus  = zeros(sizeZ,sizeRho);
	matRho_minus = zeros(sizeZ,sizeRho);
	vecOmega_plus = zeros(sizeZ,1);
	vecOmega_minus = zeros(sizeZ,1);
	for n=1:sizeZ
		epsZ_temp = vecEpsZ(n);
		%
		vecZ_temp = vecZ;
		vecZ_temp(n) = vecZ(n)+epsZ_temp;
		%
		[ errFlag, vecRho_temp ] = funchRho( rhoArgs, vecZ_temp );
		if (errFlag)
			msg_error( verbLev, thisFile, __LINE__, sprintf( ...
			  "funchRho() failed for +%0.3g to element %d.", epsZ_temp, n ) );
			retCode = RETCODE__BAD_INPUT;
			return;
		end
		if ( valLev >= VALLEV__MEDIUM )
			assert( isrealarray(vecRho_temp,[sizeRho,1]) );
		end
		%
		omega_temp = funchOmega(vecRho_temp);
		if ( valLev >= VALLEV__MEDIUM )
			assert( isrealscalar(omega_temp) );
			assert( omega_temp >= 0.0 );
		end
		%
		matRho_plus(n,:) = vecRho_temp;
		vecOmega_plus(n) = omega_temp;
		%
		%
		%
		vecZ_temp = vecZ;
		vecZ_temp(n) = vecZ(n)-epsZ_temp;
		%
		[ errFlag, vecRho_temp ] = funchRho( rhoArgs, vecZ_temp );
		if (errFlag)
			msg_error( verbLev, thisFile, __LINE__, sprintf( ...
			  "funchRho() failed for -%0.3g to element %d.", epsZ_temp, n ) );
			retCode = RETCODE__BAD_INPUT;
			return;
		end
		if ( valLev >= VALLEV__MEDIUM )
			assert( isrealarray(vecRho_temp,[sizeRho,1]) );
		end
		%
		omega_temp = funchOmega(vecRho_temp);
		if ( valLev >= VALLEV__MEDIUM )
			assert( isrealscalar(omega_temp) );
			assert( omega_temp >= 0.0 );
		end
		%
		matRho_minus(n,:) = vecRho_temp;
		vecOmega_minus(n) = omega_temp;
		%
		%
		%
		clear vecRho_temp;
		clear omega_temp;
		clear vecZ_temp;
		clear epsZ_temp;
	end
	%
	%
	%
	% Calc finite-differences.
	matRhoDZ = zeros(sizeZ,sizeRho);
	matRhoDZ2Diag = zeros(sizeZ,sizeRho);
	vecOmegaDZ = zeros(sizeZ,1);
	vecOmegaDZ2Diag = zeros(sizeZ,1);
	for n=1:sizeZ
		matRhoDZ(n,:) = (matRho_plus(n,:) - matRho_minus(n,:)) / ( 2.0*vecEpsZ(n) );
		matRhoDZ2Diag(n,:) = ...
		  ( matRho_plus(n,:) + matRho_minus(n,:) - 2.0*vecRho' ) / ( vecEpsZ(n)^2 );
		vecOmegaDZ(n) = ( vecOmega_plus(n) - vecOmega_minus(n) )/( 2.0*vecEpsZ(n) );
		vecOmegaDZ2Diag(n) = ...
		  ( vecOmega_plus(n) + vecOmega_minus(n) - 2.0*omega )/( vecEpsZ(n)^2 );
	end
	%
	%
	%
	% Calc combination stuff.
	% Note that even though H1+H2 would be more accurate locally,
	% H1 is likely(?) to be more accurate near the actual minimum.
	vecG_a = vecOmegaDZ;
	vecG_b = matRhoDZ * vecOmegaDRho;
	matH1 = matRhoDZ * matOmegaDRho2 * (matRhoDZ');
	matH2Diag = diag( matRhoDZ2Diag * vecOmegaDRho );
	matH_inexact = matH1 + matH2Diag;
	matHDiag_a = diag(vecOmegaDZ2Diag);
	matHDiag_b = diag(diag(matH_inexact));
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealarray(vecG_a,[sizeZ,1]) );
		assert( isrealarray(vecG_b,[sizeZ,1]) );
		assert( isrealarray(matH1,[sizeZ,sizeZ]) );
		assert( isrealarray(matH2Diag,[sizeZ,sizeZ]) );
		assert( isrealarray(matH_inexact,[sizeZ,sizeZ]) );
		assert( isrealarray(matHDiag_a,[sizeZ,sizeZ]) );
		assert( isrealarray(matHDiag_b,[sizeZ,sizeZ]) );
		%
		assert( sum(sum((matH1'-matH1).^2))               <= eps150*sum(sum(matH1.^2)) );
		assert( sum(sum((matH2Diag'-matH2Diag).^2))       <= eps150*sum(sum(matH2Diag.^2)) );
		assert( sum(sum((matH_inexact'-matH_inexact).^2)) <= eps150*sum(sum(matH_inexact.^2)) );
		assert( sum(sum((matHDiag_a'-matHDiag_a).^2))     <= eps150*sum(sum(matHDiag_a.^2)) );
		assert( sum(sum((matHDiag_b'-matHDiag_b).^2))     <= eps150*sum(sum(matHDiag_b.^2)) );
	end
	%
	% Set output.
	vecG = vecG_a; % May be overwritten by "__tweak".
	matH = matH1;  % May be overwritten by "__tweak".
	%
	datOut.vecZ = vecZ;
	%
	datOut.vecRho = vecRho;
	datOut.omega = omega;
	datOut.vecOmegaDRho = vecOmegaDRho;
	datOut.matOmegaDRho2 = matOmegaDRho2;
	%
	datOut.vecEpsZ = vecEpsZ;
	datOut.matRho_plus  = matRho_plus;
	datOut.matRho_minus = matRho_minus;
	datOut.vecOmega_plus  = vecOmega_plus;
	datOut.vecOmega_minus = vecOmega_minus;
	%
	datOut.matRhoDZ = matRhoDZ;
	datOut.matRhoDZ2Diag = matRhoDZ2Diag;
	datOut.vecOmegaDZ = vecOmegaDZ;
	datOut.vecOmegaDZ2Diag = vecOmegaDZ2Diag;
	%
	datOut.vecG_a = vecG_a;
	datOut.vecG_b = vecG_b;
	datOut.matH1 = matH1;
	datOut.matH2Diag = matH2Diag;
	datOut.matH_inexact = matH_inexact;
	datOut.matHDiag_a = matHDiag_a;
	datOut.matHDiag_b = matHDiag_b;
	%
	%
	%
	prm_tweak = mygetfield( prm, "prm_tweak", prm );
	[ vecG, matH, retCode, datOut_tweak ] = findBestFit1D__calcLocalModel__tweak( datOut, prm_tweak );
	datOut.datOut_tweak = datOut_tweak;
	if ( RETCODE__SUCCESS ~= retCode )
		msg_error( verbLev, thisFile, __LINE__, sprintf(
		  "__calcLocalModel__tweak() returned %s.", retcode2str(retCode) ) );
		retCode = RETCODE__BAD_INPUT;
		return;
	end
	%
	%
	%
	retCode = RETCODE__SUCCESS;
	return;
end
