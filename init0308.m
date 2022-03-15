page_screen_output(0)
pkg load optim;
addpath( "myutil" );
%
%my_editor = "/usr/bin/codeblocks";
my_editor = "/usr/bin/mousepad";
EDITOR( [my_editor ' %s'] ); % Moot if not USING_GUI; see options within GUI.
%
switch 20
case 10
	msg(__FILE__,__LINE__,"Adding path \"study/20220222findLocMin/testfunc2021/\".");
	addpath("study/20220222findLocMin/testfunc2021/");
	msg(__FILE__,__LINE__,"Adding path \"study/20220222findLocMin/0308baselineSolvers/\".");
	addpath("study/20220222findLocMin/0308baselineSolvers/");
case 20
	msg(__FILE__,__LINE__,"Adding path \"study/20220222findLocMin/0314pieIsBased/\".");
	addpath("study/20220222findLocMin/0314pieIsBased/");
otherwise
	error( "Invalid case." );
end
