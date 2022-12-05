page_screen_output(0)
%pkg load optim;
addpath( "myutil" );
%
%my_editor = "/usr/bin/codeblocks";
my_editor = "/usr/bin/mousepad";
EDITOR( [my_editor ' %s'] ); % Moot if not USING_GUI; see options within GUI.
%
switch 310
case 10
	msg(__FILE__,__LINE__,"Adding path \"study/20220222findLocMin/testfunc2021/\".");
	addpath("study/20220222findLocMin/testfunc2021/");
	msg(__FILE__,__LINE__,"Adding path \"study/20220222findLocMin/0308baselineSolvers/\".");
	addpath("study/20220222findLocMin/0308baselineSolvers/");
case 20
	msg(__FILE__,__LINE__,"Adding path \"study/20220222findLocMin/0314pieIsBased/\".");
	addpath("study/20220222findLocMin/0314pieIsBased/");
case 30
	msg(__FILE__,__LINE__,"Adding path \"study/20220222findLocMin/0317compactify/\".");
	addpath("study/20220222findLocMin/0317compactify/");
case 40
	msg(__FILE__,__LINE__,"Adding path \"study/20220222findLocMin/0320zeroingIn/\".");
	addpath("study/20220222findLocMin/0320zeroingIn/");
case 100
	msg(__FILE__,__LINE__,"Adding path \"study/20220326BeyondZero/0326ExamineMyAssumptions/\".");
	addpath("study/20220326BeyondZero/0326ExamineMyAssumptions/");
case 110
	msg(__FILE__,__LINE__,"Adding path \"study/20220326BeyondZero/0327jfnk/\".");
	addpath("study/20220326BeyondZero/0327jfnk/");
case 120
	msg(__FILE__,__LINE__,"Adding path \"study/20220326BeyondZero/0328ReturnToZero/\".");
	addpath("study/20220326BeyondZero/0328ReturnToZero/");
case 130
	msg(__FILE__,__LINE__,"Adding path \"study/20220326BeyondZero/0505ReExamine/\".");
	addpath("study/20220326BeyondZero/0505ReExamine/");
case 140
	msg(__FILE__,__LINE__,"Adding path \"study/20220326BeyondZero/0507zlinsolf/\".");
	addpath("study/20220326BeyondZero/0507zlinsolf/");
case 150
	msg(__FILE__,__LINE__,"Adding path \"study/20220326BeyondZero/0520NonPolyExpansion/\".");
	addpath("study/20220326BeyondZero/0520NonPolyExpansion/");
case 160
	msg(__FILE__,__LINE__,"Adding path \"study/20220530sxsolve/\".");
	addpath("study/20220530sxsolve/");
	msg(__FILE__,__LINE__,"Adding path \"study/20220530sxsolve/demo/\".");
	addpath("study/20220530sxsolve/demo"); % Not needed???
case 170
	msg(__FILE__,__LINE__,"Adding path \"study/20220710AndOther/\".");
	addpath("study/20220710AndOther");
case 180
	msg(__FILE__,__LINE__,"Adding path \"study/20220830TheWorks/\".");
	addpath("study/20220830TheWorks");
case 181
	msg(__FILE__,__LINE__,"Adding path \"cxcpp/20220624forPIES/post/\".");
	addpath("cxcpp/20220624forPIES/post/");
case 200
	msg(__FILE__,__LINE__,"Adding path \"study/20220909APS2022a/\".");
	addpath("study/20220909APS2022a");
case 300
	msg(__FILE__,__LINE__,"Adding path \"study/20221111hessfit/\".");
	addpath("study/20221111hessfit");
case 310
	msg(__FILE__,__LINE__,"Adding path \"study/20221204hessfit/\".");
	addpath("study/20221204hessfit");
otherwise
	error( "Invalid case." );
end
