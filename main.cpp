// pppar01.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
//20220417 增加INS功能
#include "pch.h"
#include <iostream>
#include "rtklib.h"
#include <string>

#include <sys/types.h>
#include "unistd.h"

using namespace std;
//int paid_main(void)
//{
//	int     var;
//	pid_t   pid;
//
//	var = 88;
//	if (write(STDOUT_FILENO, buf, sizeof(buf) - 1) != sizeof(buf) - 1)//write函数将buf中的内容写到标准输出流中，此处制作输出使用。
//		err_sys("write error");
//	printf("before fork\n");
//
//	if ((pid = fork()) < 0) {
//		err_sys("fork error");
//	}
//	else if (pid == 0) {
//		glob++;
//		var++;
//	}
//	else {
//		sleepms(2);  //学过Linux的都知道，在fork()之后，是父进程先执行还是子进程先执行取决于内核的调度算法。所以在这里为了让子进程先执行，我们先让父进程sleep 2秒，但是2秒钟不一定够，所以不一定能够保证子进程先执行。
//	}
//}

int main(int argc, char** argv)
{
	char** in_str = NULL;

	std::cout << "Hello World!\n";
	//char *file="C://Users//win7//source//repos//ppparmaster//GNSS_DATA//GNSS_DATA//atx//igs14_2164.atx";
	//pcvs_t pcvs = { 0 };
	//readantex(file, &pcvs);
	//cout << argc << ":  " << argv << endl;
	//int argc = 8;
	//char* argv[] = { "-C","H:/winmove/project/pppar/pppar_main/conf/PPP/ppp_mgex_2021_3f.conf",
	//		"-S", "GREC", "-M", "PPP-KINE", "-A", "0", "L", "0" };
	if (1) {
		//pppar_main(/*argc, &argv*/);
		pppar_main(argc, argv);
	}
	else {
		rtknavi_main();
		while (true)
		{
			TimerTimer(in_str);
			// cout << in_str << endl;
			sleepms(1000);
		}
	}
	return 1;
}
// 入门提示:
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件