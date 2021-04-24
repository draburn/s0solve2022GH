if (0)
clear;
assert( 0.0 == pickNextPoint_1D_0423 )
clear;
%
clear;
assert( abs(pickNextPoint_1D_0423(0,0,0)-sqrt(sqrt(eps))) < sqrt(eps) );
clear;
%
clear;
funch_ef = @(x)([ (0>=x)+(1<=x); x./((0<=x).*(1>=x)) ]);
fa_x = linspace(-2,3,50);
foobar = funch_ef(fa_x);
ia_exceptionFlag = foobar(1,:);
fa_f = foobar(2,:);
clear;
end
%
clear;
setprngstates(0);
funch_exceptionFlag = @(x)( (0>=x)+(1<=x) );
funch_f = @(x)( x./((0<=x).*(1>=x)) );
fa_x = [ 0.4, -0.2, 0.6, 10.0, 1.3 ];
ia_exceptionFlag = funch_exceptionFlag(fa_x);
fa_f = funch_f(fa_x);
pickNextPoint_1D_0423( fa_x, ia_exceptionFlag, fa_f );
