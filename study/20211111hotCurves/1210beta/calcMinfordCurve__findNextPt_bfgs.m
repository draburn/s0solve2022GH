function [ vecX, datOut ] = calcMinfordCurve__findNextPt_bfgs( funchOmega, funchG, onSurf0, vecX0, vecXC, bigR, matS=[], prm=[] )
	thisFile = "calcMinfordCurve__findNextPt_bfgs";
	msg( thisFile, __LINE__, "WARNING: THIS APPROACH DOESN'T WORK BECAUSE THE GRADEINT IS DISCONTINUOUS!" );
	%msg( thisFile, __LINE__, "HACK-SACK!" );
	%
	iterLimit = 100;
	verbLev = 0;
	%%%in_args = { vecX0, vecXC, bigR, matS, funchOmega };
	%%%in_ctrl = { iterLimit, verbLev };
	%%%[ out_x, out_obj_value, out_cnvg, out_iters ] = bfgsmin( "bfgsmin_funchOmegaBall", in_args, in_ctrl )
	%
	cnvgCrit = 1; % Default is 1 "strong", alt is 0 "weak"
	funchArgIndex = 1;
	memLim = 0; % none, use regular bfgs
	fTol = 1e-12; % Default is 1e-12
	dTol = 1e-2; % Default is 1e-6
	gTol = 1e-5; % Default is 1e-5
	in_args = { vecX0, vecXC, bigR, matS, funchOmega, funchG };
	in_ctrl = { iterLimit, verbLev, cnvgCrit, funchArgIndex, memLim, fTol, dTol, gTol };
	[ out_x, out_obj_value, out_cnvg, out_iters ] = bfgsmin( "bfgsmin_funchOmegaGBall", in_args, in_ctrl );
	%
	%vecX = vecX0;
	%return
	%
	vecX = out_x;
	%
	return;
end
