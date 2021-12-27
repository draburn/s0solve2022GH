clear;
numFigs = 0;
%
%funchOmega = @(x)( (x-pi).^e );
%funchG = @(x)( e*(x-pi).^(e-1) );
%function [ omega, vecG ] = evalOmegaG( vecX )
function [ omega ] = evalOmegaG( vecX )
	omegaRoot = 0.1;
	vecXRoot = [ 3.0; 1.0 ];
	matH = [ 1.0, 2.0; 2.0, 5.0 ];
	vecY = vecX - vecXRoot;
	omega = omegaRoot + 0.5*(vecY'*matH*vecY);
	%vecG = matH*vecY;
endfunction
%
%
if (0)
xVals = linspace(0.0,5.0,1001);
numFigs++; figure(numFigs);
plot( xVals, evalOmega(xVals), 'o-' );
hold on;
plot( xVals, 0*xVals, 'k-' );
hold off;
grid on;
end
%
vecX0 = [ 0.0; 0.0 ];
iterLimit = 100;
verbLev = 1;
%
in_args = { vecX0 };
in_ctrl = { iterLimit, verbLev };
[ out_x, out_obj_value, out_cnvg, out_iters ] = bfgsmin( "evalOmegaG", in_args, in_ctrl );
