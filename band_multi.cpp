#include "band_multi.h"

cx_vec triMulti(const cx_vec& a, const cx_vec& b, const cx_vec& c, int q, const cx_vec& x)
{
	int n = x.n_elem;
	cx_vec y(n,fill::zeros);

	for (int i = 0; i < q;i++)
	{
		y(i) = b(i) * x(i) + c(i) * x(i + q);
	}
	for (int i = q; i < n - q; i++)
	{
		y(i) = a(i) * x(i - q) + b(i) * x(i) + c(i) * x(i + q);
	}
	for (int i = n - q;i < n;i++)
	{
		y(i) = a(i) * x(i - q) + b(i) * x(i);
	}

	return y;
}

cx_vec fiveMulti(const cx_vec& a, const cx_vec& b, const cx_vec& c, const cx_vec& d, const cx_vec& e, int q, const cx_vec& x)
{
	return triMulti(a, c, e, q + 1, x) +
		triMulti(b, cx_vec(size(c),fill::zeros), d, q - 1, x);
}
