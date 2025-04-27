#include "../include/thomas.h"

cx_vec thomasAlgorithm(const cx_vec& a, const cx_vec& b, const cx_vec& c, const cx_vec& d)
{
	int n = d.size();
	//% a, b, c, d 同列
	cx_vec  g(n);
	cx_vec  r(n); //% 存放中间d

	g(0) = c(0) / b(0);
	r(0) = d(0) / b(0);

	//% 正向
	for (int i = 1; i < n; i++)
	{
		g(i) = c(i) / (b(i) - a(i) * g(i - 1));
		r(i) = (d(i) - a(i) * r(i - 1)) / (b(i) - a(i) * g(i - 1));
	}

	//% 反向
	cx_vec x(n);
	x(n - 1) = r(n - 1);

	for (int i = n - 2;i > -1;i--)
	{
		x(i) = r(i) - g(i) * x(i + 1);
	}

	return x;
}

cx_vec thomasAlgorithm(const cx_vec& a, const cx_vec& b, const cx_vec& c, int p, const cx_vec& d)
{

	//n = length(d);
	int n = d.size();
	//% a, b, c, d 同列
	cx_vec	g(n,fill::zeros);
	cx_vec r(n, fill::zeros);     //% 存放中间d

	g.rows(0, p - 1) = c.rows(0, p - 1) / b.rows(0, p - 1);
	r.rows(0, p - 1) = d.rows(0, p - 1) / b.rows(0, p - 1);

	//% 正向

	for (int i = p;i < n;i++)
	{
		g(i) = c(i) / (b(i) - a(i) * g(i - p));
		r(i) = (d(i) - a(i) * r(i - p)) / (b(i) - a(i) * g(i - p));
	}


	//% 反向
	cx_vec x(n);
	x.rows(n - p, n - 1) = r.rows(n - p, n - 1);

	for (int i = n - p - 1;i > -1;i--)
	{
		x(i) = r(i) - g(i) * x(i + p);
	}

	return x;
}
