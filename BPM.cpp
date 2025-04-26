#include "BPM.h"

void BPM::init()
{
	layersPML = 10;
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

	sx = ones(nt) + 0.0 * iu;
	sy = ones(nt) + 0.0 * iu;
	isx = 1 / sx;
	isy = 1 / sy;

	Ex.set_size(nt, nz);
	Ey.set_size(nt, nz);

}

void BPM::compute_PML()
{
	cout << "\t初始化PML参数\n";
	cx_mat Sx(nx, ny, fill::ones), Sy(nx, ny, fill::ones);
	int m = 3;
	double omega = c0 * k0;
	double R0 = 1e-6;
	double sigmaXmax = -(m + 1) * log10(R0) / (2 * 373 * dx * layersPML);
	double sigmaYmax = -(m + 1) * log10(R0) / (2 * 373 * dy * layersPML);


	cx_double value = 0;
	for (int i = 0;i < layersPML + 1;i++) {

		value = 1.0 + sigmaXmax * pow((layersPML - i + 0.0) / layersPML, m) / (iu * omega * eps0);
		Sx.row(i) = cx_rowvec(ny, fill::value(value));
		Sx.row(nx - 1 - i) = Sx.row(i);

		value = 1.0 + sigmaYmax * pow((layersPML - i + 0.0) / layersPML, m) / (iu * omega * eps0);
		Sy.col(i) = cx_vec(nx, fill::value(value));
		Sy.col(ny - 1 - i) = Sy.col(i);

	}
	sx = Sx.as_col();
	sy = Sy.as_col();
	isx = 1 / sx;
	isy = 1 / sy;
}

void BPM::compute_Matrix()
{
	cout << "\t计算Ax Ay Bx By C D";

	Ax_V.set_size(nz);
	Ay_V.set_size(nz);
	Bx_V.set_size(nz);
	By_V.set_size(nz);
	C_V.set_size(nz);
	D_V.set_size(nz);
	clock_t t1 = clock();
#pragma omp parallel for
	for (int i = 0; i < nz;i++)
	{
		cx_vec erx = (*dev)["Epsx"].slice(i).as_col() + 0.0 * iu;
		cx_vec erz = (*dev)["Epsx"].slice(i).as_col() + 0.0 * iu;
		cx_vec ery = (*dev)["Epsy"].slice(i).as_col() + 0.0 * iu;
		cx_vec erxy = (*dev)["Epsxy"].slice(i).as_col() + 0.0 * iu;
		cx_vec eryx = (*dev)["Epsyx"].slice(i).as_col() + 0.0 * iu;


		dxdxfunc(isx, sx % erz, erx, nx, ny, dx, dy, Ax_V(i));
		Ax_V(i).diagArr.col(1) += 0.5 * k0 * k0 * (erx - n0 * n0);

		dydyfunc(isy, sy, cx_vec(nx * ny, fill::ones), nx, ny, dx, dy, Ay_V(i));
		Ay_V(i).diagArr.col(1) += 0.5 * k0 * k0 * (erx - n0 * n0);

		dxdxfunc(isx, sx, cx_vec(nx * ny, fill::ones), nx, ny, dx, dy, Bx_V(i));
		Bx_V(i).diagArr.col(1) += 0.5 * k0 * k0 * (ery - n0 * n0);

		dydyfunc(isy, sy % erz, ery, nx, ny, dx, dy, By_V(i));
		By_V(i).diagArr.col(1) += 0.5 * k0 * k0 * (ery - n0 * n0);

		DiagStruct C1, C2;
		dxdyfunc(isx, sy % erz, ery, nx, ny, dx, dy, C1);
		dxdyfunc(isx, sy, cx_vec(nx * ny, fill::ones), nx, ny, dx, dy, C2);
		C_V(i).diagArr = C1.diagArr - C2.diagArr;
		C_V(i).diagArr.col(2) = erxy * k0 * k0;
		C_V(i).diagIndex = C1.diagIndex;

		DiagStruct C3, C4;
		dydxfunc(isy, sx % erz, erx, nx, ny, dx, dy, C3);
		dydxfunc(isy, sx, cx_vec(nx * ny, fill::ones), nx, ny, dx, dy, C4);
		D_V(i).diagArr = C3.diagArr - C4.diagArr;
		D_V(i).diagArr.col(2) = erxy * k0 * k0;
		D_V(i).diagIndex = C3.diagIndex;

		// d b c a排序的
	}
	clock_t t2 = clock();
	cout << "\t:\t" << (double)(t2 - t1) / CLOCKS_PER_SEC << "s" << endl;


}










void BPM::Qusi_TM_Propagate()
{



}

void BPM::Qusi_TE_Propagate()
{

}

void BPM::FullVector_propagate()
{
	cout << "\t全矢量传播\n";

	Ex.col(0) = vectorise((*dev)["Exin"]) + 0.0 * iu;
	Ey.col(0) = vectorise((*dev)["Eyin"]) + 0.0 * iu;

	cx_double a_ = (1 - alpha) * dz / 2 / 1i / n0 / k0;
	cx_double b_ = -alpha * dz / 2 / 1i / n0 / k0;

	clock_t t1 = clock();

	for (int i = 0; i < nz - 1; i++) {

		cout << "\t第" << i + 1 << "/" << nz - 1 << "层";
		// CN差分1.1
		cx_vec d = spdiags(join_rows(a_*Ay_V(i).diagArr.col(0) ,	 1 + a_ * Ay_V(i).diagArr.col(1), a_ * Ay_V(i).diagArr.col(2)),
			Ay_V(i).diagIndex, nt, nt).st() * Ex.col(i);
		cx_vec c = b_*Ay_V(i + 1).diagArr.col(0);
		cx_vec b = 1.0 + b_ * Ay_V(i + 1).diagArr.col(1); //对角
		cx_vec a = b_*Ay_V(i + 1).diagArr.col(2);
		cx_vec utmp = thomas_algorithm(a, b, c,nx, d);

		//// CN差分1.2
		d = spdiags(join_rows(a_ * By_V(i).diagArr.col(0), 1 + a_ * By_V(i).diagArr.col(1), a_ * By_V(i).diagArr.col(2)),
			By_V(i).diagIndex, nt, nt).st() * Ey.col(i)
		+spdiags( a_*D_V(i).diagArr, D_V(i).diagIndex, nt, nt).st() * Ex.col(i)
		- spdiags( b_*D_V(i+1).diagArr, D_V(i+1).diagIndex, nt, nt).st() * utmp;
		c = b_* By_V(i + 1).diagArr.col(0);
		b = 1.0 + b_ * By_V(i + 1).diagArr.col(1); //对角
		a = b_* By_V(i + 1).diagArr.col(2);
		cx_vec vtmp = thomas_algorithm(a, b, c,nx, d);


		// CN差分2.1
		d = spdiags(
			join_rows(a_ * Bx_V(i).diagArr.col(0), 1 + a_ * Bx_V(i).diagArr.col(1), a_ * Bx_V(i).diagArr.col(2)), Bx_V(i).diagIndex, nt, nt).st() * vtmp;

		c = b_ * Bx_V(i + 1).diagArr.col(0);
		b = 1.0+b_ * Bx_V(i + 1).diagArr.col(1);
		a = b_ * Bx_V(i + 1).diagArr.col(2);
		
		Ey.col(i + 1) = thomas_algorithm(a, b, c, d);

	 
		// CN差分
		d = spdiags(join_rows(a_ * Ax_V(i).diagArr.col(0), 1 + a_ * Ax_V(i).diagArr.col(1), a_ * Ax_V(i).diagArr.col(2)),
			Ax_V(i).diagIndex, nt, nt).st()* utmp
			+ spdiags(a_ * C_V(i).diagArr, C_V(i).diagIndex, nt, nt).st() * vtmp
			- spdiags(b_ * C_V(i + 1).diagArr, C_V(i + 1).diagIndex, nt, nt).st() * Ex.col(i + 1);

		c = b_ * Ax_V(i + 1).diagArr.col(0);
		b = 1.0 + b_ * Ax_V(i + 1).diagArr.col(1);
		a = b_ * Ax_V(i + 1).diagArr.col(2);

		Ex.col(i + 1) = thomas_algorithm(a, b, c, d);

		cout << "\r";
		cout.flush();
	}
	clock_t t2 = clock();
	cout << "\t:\t" << (double)(t2 - t1) / CLOCKS_PER_SEC << "s" << endl;


}

void BPM::postData()
{

	system("del output/bpmOut.h5");

	cout << "\n写入文件中：" << endl;

	string filePath = "output/bpmOut.h5";

	(*dev)["x"].save(hdf5_name(filePath, "x"));
	(*dev)["y"].save(hdf5_name(filePath, "y", hdf5_opts::append));
	(*dev)["z"].save(hdf5_name(filePath, "z", hdf5_opts::append));


	
	mat Ex_abs = abs(Ex);
	Ex_abs.save(hdf5_name(filePath, "Ex_abs", hdf5_opts::append));
	mat Ey_abs = abs(Ey);
	Ey_abs.save(hdf5_name(filePath, "Ey_abs", hdf5_opts::append));


	mat Ex_real = real(Ex);
	Ex_real.save(hdf5_name(filePath, "Ex_real", hdf5_opts::append));

	mat Ey_real = real(Ey);
	Ey_real.save(hdf5_name(filePath, "Ey_real", hdf5_opts::append));

}