function [ vecX, datOut ] = calcMinfordCurve__findNextPt_mybfgs( funchOmega, funchG, onSurf0, vecX0, vecXC, bigR, matS=[], prm=[] )
	thisFile = "calcMinfordCurve__findNextPt_mybfgs";
	%
	%vecXfoo = [ pi; e ]
	%
	pullCoeff = 0.01;
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
		vecG = ( vecG0 - vecT*(vecD'*vecG0) )*(bigR/s) + pullCoeff*( vecX - vecXSurf );
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
