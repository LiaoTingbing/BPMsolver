#include "BPM.h"

void BPM::init()
{
	nx = (*dev)["x"].n_elem;
	ny = (*dev)["y"].n_elem;
	nz = (*dev)["z"].n_elem;
	nt = nx * ny;

	cout << "nx = " << nx << endl;
	cout << "ny = " << ny << endl;
	cout << "nz = " << nz << endl;
	cout << "nt = " << nt << endl;

	dx = (*dev)["x"](1) - (*dev)["x"](0);
	dy = (*dev)["y"](1) - (*dev)["y"](0);
	dz = (*dev)["z"](1) - (*dev)["z"](0);
	lambda = (*dev)["lambda"](0);
	k0 = 2 * pi / lambda;
	n0 = (*dev)["neff"](0);
	alpha = 0.5; // CN 差分控制参数 

	cout << "dx = " << dx << endl;
	cout << "dy = " << dy << endl;
	cout << "dz = " << dz << endl;
	cout << "lambda = " << lambda << endl;
	cout << "k0 = " << k0 << endl;
	cout << "n0 = " << n0 << endl;
	cout << "alpha = " << alpha << endl;


	Ein = join_cols(vectorise((*dev)["Exin"]), vectorise((*dev)["Eyin"]));

	Eout.set_size( 2* nt , nz );
}

sp_mat BPM::Calculate_P(int i)
{
	// 计算P矩阵
	//field<sp_cx_mat> Ptrans(nz);

	// 实数折射率

	//int i = 0;
	vec Iv(nt, fill::ones);
	//sp_cx_mat Idiv = spdiags(cx_vec(nt, fill::ones), ivec{ 0 }, nt, nt);

	vec ERX, ERY, ERZ, ERXY(nt), ERYX(nt);  // 如各向异性对ERXY，ERYX赋值

	ERX =pow( vectorise((*dev)["index_x"].slice(i)) , 2.0);
	ERY =pow( vectorise((*dev)["index_y"].slice(i)) ,2);
	ERZ =pow( vectorise((*dev)["index_z"].slice(i)) ,2 );

	sp_mat Pxx = dxdx(Iv, ERZ, ERX, nx, ny, dx, dy)
		+ dydy(Iv, Iv, Iv, nx, ny, dx, dy)
		+  spdiags(k0 * k0 * (ERX / eps0 - n0 * n0), ivec{ 0 }, nt, nt);

	sp_mat Pyy = dxdx(Iv, Iv, Iv, nx, ny, dx, dy)
		+ dydy(Iv, ERZ, ERY, nx, ny, dx, dy)
		+ spdiags(k0 * k0 * (ERY / eps0 - n0 * n0), ivec{ 0 }, nt, nt);

	sp_mat Pxy = dxdy(Iv, ERZ, ERY, nx, ny, dx, dy)
		- dxdy(Iv, Iv, Iv, nx, ny, dx, dy)
		+ spdiags(ERXY / eps0 * k0 * k0, ivec{ 0 }, nt, nt);

	sp_mat Pyx = dydx(Iv, ERZ, ERX, nx, ny, dx, dy)
		- dydx(Iv, Iv, Iv, nx, ny, dx, dy)
		+ spdiags(ERYX / eps0 * k0 * k0, ivec{ 0 }, nt, nt);

	sp_mat P(2 * nt, 2 * nt);
	P(0, 0, size(nt, nt)) = Pxx;
	P(0, nt, size(nt, nt)) = Pxy;
	P(nt, 0, size(nt, nt)) = Pyx;
	P(nt, nt, size(nt, nt)) = Pyy;

	return P;
}

void BPM::propagate()
{
	Eout.col(0) = Ein + 0.0 * iu;   // 第一个场值
	cx_vec Iv(nt, fill::ones);
	sp_cx_mat II = spdiags(Iv * (2.0 * iu * n0 * k0 / dz), ivec{ 0 }, 2 * nt, 2 * nt);

	// 沿着层传播
	sp_mat P1, P2;
	P1 = Calculate_P(0);

	for (int i = 0; i < nz - 1; i++) {

		cout << "传输\t:\t" << i+1 <<"~"<<i+2 << "层" << endl;
		
		P2 = Calculate_P(i + 1);

		// solve CN 差分 2in0k0/dz[ u2 - u1] = alpha*P2*u2 + (1 - alpha)*p1*u1
		Eout.col(i + 1) = spsolve(II - alpha * P2, (II + (1 - alpha) * P1) * Eout.col(i));

		P1 = P2; 
		//Eout(i + 1).print();
	}
}

void BPM::postData()
{

	//cx_cube Eout_cube(2 * nt, nz, 1);
	//Eout_cube.slice(0) = Eout;
	//Eout_cube.resize(nx, ny, nz);

	string filePath = "lumerical/bpmOut.h5";

	mat Eout_real = real(Eout);
	mat Eout_imag = imag(Eout);
	Eout_real.save(hdf5_name(filePath, "Eout_real"));
	Eout_imag.save(hdf5_name(filePath, "Eout_imag", hdf5_opts::append));

	cx_cube Ex(nt, nz, 1);
	Ex.slice(0) = Eout(0, 0, size(nt, nz));
	Ex.resize(nx, ny, nz);
	cube Ex_real = real(Ex);

	Ex_real.save(hdf5_name(filePath, "Ex_real", hdf5_opts::append));

	(*dev)["x"].save(hdf5_name(filePath, "x", hdf5_opts::append));
	(*dev)["y"].save(hdf5_name(filePath, "y", hdf5_opts::append));
	(*dev)["z"].save(hdf5_name(filePath, "z", hdf5_opts::append));

	(*dev)["index_x"].slice(nz-5).save(hdf5_name(filePath, "index_x", hdf5_opts::append));


}