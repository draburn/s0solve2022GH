function [ vecX, datOut ] = calcMinfordCurve__findNextPt_mybfgs( funchOmega, funchG, onSurf0, vecX0, vecXC, bigR, matS=[], prm=[] )
	thisFile = "calcMinfordCurve__findNextPt_mybfgs";
	msg( thisFile, __LINE__, "WARNING: THIS APPROACH DOESN'T WORK BECAUSE THE GRADEINT IS DISCONTINUOUS!" );
	%
	%vecXfoo = [ pi; e ]
	%
	%%%pullCoeff = 0.01;
	pullCoeff = 10.0;
	%msg( thisFile, __LINE__, "Hey?" );
	function [ omega ] = funcOmegaBall( vecX )
		%thisSubfile = "calcMinfordCurve__findNextPt_mybfgs::funcOmegaBall";
		%msg( thisSubfile, __LINE__, "Welcome!" );
		%echo__vecXC = vecXC
		%echo__bigR = bigR
		%echo__matS = matS
		%echo__vecX = vecX
		vecD = vecX - vecXC;
		if ( isempty(matS) )
			s = norm(vecD);
		else
			s = norm(matS*vecD);
		end
		if ( s <= bigR )
			omega = funchOmega( vecX );
			%msg( thisSubfile, __LINE__, "Goodbye!" );
			return;
		end
		vecXSurf = vecXC + (vecD*bigR/s);
		omega = funchOmega( vecXSurf ) + 0.5*pullCoeff*(norm(vecX - vecXSurf)^2);
		%msg( thisSubfile, __LINE__, "Goodbye!" );
		return;
	endfunction
	%msg( thisFile, __LINE__, "You listen?" );
	%
	%
	%%
	function [ omega ] = funcOmega_forceSurf( vecX )
		vecD = vecX - vecXC;
		if ( isempty(matS) )
			s = norm(vecD);
		else
			s = norm(matS*vecD);
		end
		if ( s <= eps*bigR )
			omega = funchOmega( vecX );
			%msg( thisSubfile, __LINE__, "Goodbye!" );
			return;
		end
		vecXSurf = vecXC + (vecD*bigR/s);
		omega = funchOmega( vecXSurf ) + 0.5*pullCoeff*(norm(vecX - vecXSurf)^2);
		return;
	endfunction
	%%
	%
	%
	function [ vecG ] = funcGBall( vecX )
		%thisSubfile = "calcMinfordCurve__findNextPt_mybfgs::funcGBall";
		%msg( thisSubfile, __LINE__, "Welcome!" );
		%echo__vecXC = vecXC
		%echo__bigR = bigR
		%echo__matS = matS
		%echo__vecX = vecX
		vecD = vecX - vecXC;
		if ( isempty(matS) )
			s = norm(vecD);
		else
			s = norm(matS*vecD);
		end
		if ( s <= bigR )
			%omega = funchOmega( vecX );
			vecG = funchG( vecX );
			%msg( thisSubfile, __LINE__, "Goodbye!" );
			return;
		end
		%
		vecXSurf = vecXC + (vecD*bigR/s);
		%omega = funchOmega( vecXSurf ) + 0.5*pullCoeff*(norm(vecX - vecXSurf)^2);
		vecG0 = funchG( vecXSurf );
		if (isempty(matS))
			vecT = vecD/(s^2);
		else
			vecT = matS'*(matS*vecD)/(s^2);
		end
		%
		switch (2)
		case 1
		vecG = ( vecG0 - vecT*(vecD'*vecG0) )*(bigR/s) + pullCoeff*( vecX - vecXSurf );
		case 2
		vecG_omega = ( vecG0 - vecT*(vecD'*vecG0) )*(bigR/s);
		dSq = vecD'*vecD;
		vecG_pull = pullCoeff*(1.0-bigR/s)*( (1.0-bigR/s)*vecD + (bigR*dSq/s)*vecT );
		vecG = vecG_omega + vecG_pull;
		otherwise
		error("Invalid case.");
		end
		%
		doFullTest = false;
		if (doFullTest)
			sizeX = size(vecX,1);
			if ( isempty(matS) )
				matS_nonEmpty = eye(sizeX,sizeX);
			else
				matS_nonEmpty = matS;
			end
			%
			vecG_test = zeros(sizeX,1);
			for n=1:sizeX
				epsFD = 1e-5;
				vecXP = vecX;
				vecXM = vecX;
				vecXP(n) += epsFD;
				vecXM(n) -= epsFD;
				switch (2)
				case 1
				vecXP_surf = vecXC + bigR*(vecXP-vecXC)/norm(matS_nonEmpty*(vecXP-vecXC));
				vecXM_surf = vecXC + bigR*(vecXM-vecXC)/norm(matS_nonEmpty*(vecXM-vecXC));
				omegaP_surf = funchOmega(vecXP_surf);
				omegaM_surf = funchOmega(vecXM_surf);
				vecG_test(n) = ( omegaP_surf - omegaM_surf ) / (2.0*epsFD);
				case 2
				omegaP = funcOmega_forceSurf(vecXP);
				omegaM = funcOmega_forceSurf(vecXM);
				vecG_test(n) = ( omegaP - omegaM ) / (2.0*epsFD);
				end
				clear epsFD;
				clear vecXP;
				clear vecXM;
				clear vecXP_surf;
				clear vecXM_surf;
				clear omegaM_surf;
				clear omegaP_surf;
			end
			clear n;
			if ( norm(vecG-vecG_test) > (1e-3)*( norm(vecG) + norm(vecG_test) ) )
			%if (1)
				echo__bigR = bigR
				echo__vecX = vecX
				echo__vecG_old =  ( vecG0 - vecT*(vecD'*vecG0) )*(bigR/s) + pullCoeff*( vecX - vecXSurf )
				echo__vecG = vecG
				echo__vecG_test = vecG_test
				echo__norm_vecG = norm(vecG)
				echo__norm_vecG_test = norm(vecG_test)
				echo__norm_diff = norm(vecG-vecG_test)
				echo__norm_vsOld = norm(vecG-echo__vecG_old)
				%error( "HALT." );
			end
			assert( norm(vecG-vecG_test) <= (1e-3)*( norm(vecG) + norm(vecG_test) ) );
		end
		%msg( thisSubfile, __LINE__, "Goodbye!" );
		return;
	endfunction
	%
	funchOmegaBall = @(x)( funcOmegaBall(x) );
	funchGBall = @(x)( funcGBall(x) );
	
	if (0)
	funchGBall( [3.0; 3.0] )
	funchGBall( [2.9001; 3.0] )
	funchGBall( [2.8999; 3.0] )
	theta = 0.01021287022984
	funchGBall( [3.0;3.0]-bigR*(1.0+(eps^0.75))*[cos(theta);sin(theta)] )
	echo__theta = theta
	error( "Halt!" );
	end
	
	
	%
	%msg( thisFile, __LINE__, "Before findMin_bfgs()." );
	vecX = findMin_bfgs( vecX0, funchOmegaBall, funchGBall );
	%%%error( "Halt!" )
	%msg( thisFile, __LINE__, "Back from findMin_bfgs()." );
	datOut = [];
	%
	return;
end
