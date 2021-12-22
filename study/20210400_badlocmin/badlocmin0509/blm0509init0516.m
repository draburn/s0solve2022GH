% Setup to looking at following "root-ish curve".
myclear;
setprngstates(0);
%
c1 = 1.0;
cx = 0.1;
c0 = (c1^2)*(1.0+cx)/4.0;
%
funchRot1 = @(x1,x2,t)( cos(t)*x1 - sin(t)*x2 );
funchRot2 = @(x1,x2,t)( sin(t)*x1 + cos(t)*x2 );
if (1)
	theta1 = 15.0*pi/180.0;
	funchF1A = @(x1,x2)( x1.*(c0+(c1*x1)+(x1.^2)) );
	funchF1B = @(x1,x2)( funchF1A(funchRot1(x1,x2,theta1),funchRot2(x1,x2,theta1)) );
	theta2 = 30.0*pi/180.0;
	funchF2A = @(x1,x2)( (-x1-5*(x2.^2))./(1.0+(x1.^2)+(x2.^2)) );
	funchF2B = @(x1,x2)( funchF2A(funchRot1(x1,x2,theta2),funchRot2(x1,x2,theta2)) );
	%%matA = [ 1.0, -0.01; 10.0, 1.0 ];
	%matA = [ 1.0, 0.0; 0.0, 1.0 ];
	matA = [ 1.0, 1.0; -2.0, 1.0 ];
	funchF = @(x)(  matA * [ funchF1B(x(1,:),x(2,:)); funchF2B(x(1,:),x(2,:)) ]);
elseif (1)
	theta1 = 0.0*pi/180.0;
	funchF1A = @(x1,x2)( x1.*(c0+(c1*x1)+(x1.^2)) );
	funchF1B = @(x1,x2)( funchF1A(funchRot1(x1,x2,theta1),funchRot2(x1,x2,theta1)) );
	theta2 = 0.0*pi/180.0;
	funchF2A = @(x1,x2)( (-x1-5*(x2.^2))./(1.0+(x1.^2)+(x2.^2)) );
	funchF2B = @(x1,x2)( funchF2A(funchRot1(x1,x2,theta2),funchRot2(x1,x2,theta2)) );
	matA = [ 1.0, 0.0; 10.0, 1.0 ];
	funchF = @(x)(  matA * [ funchF1B(x(1,:),x(2,:)); funchF2B(x(1,:),x(2,:)) ]);
elseif (1)
	theta1 = 15.0*pi/180.0;
	funchF1A = @(x1,x2)( x1.*(c0+(c1*x1)+(x1.^2)) );
	funchF1B = @(x1,x2)( funchF1A(funchRot1(x1,x2,theta1),funchRot2(x1,x2,theta1)) );
	theta2 = 0.0*pi/180.0;
	funchF2A = @(x1,x2)( x1+5.0*(x2.^2) );
	funchF2B = @(x1,x2)( funchF2A(funchRot1(x1,x2,theta2),funchRot2(x1,x2,theta2)) );
	matA = [ 1.0, 0.0; 0.0, 1.0 ];
	funchF = @(x)(  matA * [ funchF1B(x(1,:),x(2,:)); funchF2B(x(1,:),x(2,:)) ] );
else
	funchF1 = @(x1,x2)( x1.*(c0+(c1*x1)+(x1.^2)) );
	funchF2 = @(x1,x2)( x2 );
	funchF = @(x)([ funchF1(x(1,:),x(2,:)); funchF2(x(1,:),x(2,:)) ]);
end
