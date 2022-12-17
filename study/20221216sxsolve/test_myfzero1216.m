clear;
%funchFC_arrayForm = @(x)( [ x.^3 - 2.0, 1.0, 3.0*x.^2, 0.0 ] );
funchFC_arrayForm = @(x)( [ x.^3 - 2.0, 1.0-x, 3.0*x.^2, -1.0+(0*x) ] );
funchFC = @(x)( myfzero1216_repackFunc( x, funchFC_arrayForm ) );
myfzero1216( funchFC, [0.0,2.0], prm )

