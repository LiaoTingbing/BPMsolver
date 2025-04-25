
#include "common.h"

//#include "BPM.h"
#include "thomas.h"



// 测试托马斯算法


int main()
{
	cx_vec a(10000000, fill::randu) ;
	cx_vec b(10000000 , fill::randu) ;
	cx_vec c(10000000 , fill::randu) ;
	cx_vec d(10000000 , fill::randu);

	clock_t t1 = clock();
	cx_vec ans = thomas_algorithm(a, b, c, d);
	//ans.print();
	clock_t t2 = clock();
	cout << "总时间：\t" << (double)(t2 - t1) / CLOCKS_PER_SEC << "s" << endl;


	  t1 = clock();
	  ans = thomas_algorithm(a, b, c, 1,d);
	//ans.print();
	  t2 = clock();

	cout << "总时间：\t" << (double)(t2 - t1) / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}