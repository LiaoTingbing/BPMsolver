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
	//cout << "k0 = " << k0 << endl;
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
		C_V(i).pos = C1.pos;


		DiagStruct C3, C4;
		dydxfunc(isy, sx % erz, erx, nx, ny, dx, dy, C3);
		dydxfunc(isy, sx, cx_vec(nx * ny, fill::ones), nx, ny, dx, dy, C4);
		D_V(i).diagArr = C3.diagArr - C4.diagArr;
		D_V(i).diagArr.col(2) = erxy * k0 * k0;
		D_V(i).pos = C3.pos;


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


void BPM::FullVector_propagate_simple()
{

	cout << "\t全矢量传播\n";

	Ex.col(0) = vectorise((*dev)["Exin"]) + 0.0 * iu;
	Ey.col(0) = vectorise((*dev)["Eyin"]) + 0.0 * iu;

	clock_t t1 = clock();

	for (int i = 0; i < nz - 1; i++) {
		cout.flush();
		cout << "\r\t" << i + 1 << "/" << nz - 1 << "层";

		//求解CN差分方程
		cx_double a_ = (1 - alpha) * dz / 2 / 1i / n0 / k0;
		cx_double b_ = -alpha * dz / 2 / 1i / n0 / k0;
		cx_vec uout, vout;
		CNsolve(
			a_, b_,
			Ay_V(i), Ay_V(i + 1),
			By_V(i), By_V(i + 1),
			Ax_V(i), Ax_V(i + 1),
			Bx_V(i), Bx_V(i + 1),
			C_V(i), C_V(i + 1),
			D_V(i), D_V(i + 1),
			Ex.col(i), Ey.col(i),
			uout, vout
		);
		Ex.col(i + 1) = uout;
		Ey.col(i + 1) = vout;
	}
	clock_t t2 = clock();
	cout << "\t\t\t" << (double)(t2 - t1) / CLOCKS_PER_SEC << "s" << endl;

}

void BPM::FullVector_WideAngle_propagate_simple(int order)
{
	cout << "\t全矢量传播\t广角阶数\t"<<order<<"\n";

	Ex.col(0) = vectorise((*dev)["Exin"]) + 0.0 * iu;
	Ey.col(0) = vectorise((*dev)["Eyin"]) + 0.0 * iu;

	clock_t t1 = clock();

	// 第i层传播
	for (int i = 0; i < nz - 1; i++) {
		cout.flush();
		cout << "\r\t" << i + 1 << "/" << nz - 1 << "层";

		// dz放里面可能不均匀的dz和变化的n0
		cx_vec coffUp = Nn(order) - iu * n0 * k0 * dz * (1 - alpha) * Mn(order);
		cx_vec coffDown = Nn(order) + iu * n0 * k0 * dz *  alpha * Mn(order);

		cx_vec a = -1.0 / roots(coffUp)/n0/n0/k0/k0 ;
		cx_vec b = -1.0 / roots(coffDown) / n0 / n0 / k0 / k0;

		cx_vec uout, vout,uin,vin;

		//广角循环
		uin = Ex.col(i) ;
		vin = Ey.col(i) ;
		for (int j = 0;j < Nn(order).size()-1;j++) {
			//cx_double a_ = (1 - alpha) * dz / 2 / 1i / n0 / k0;
			//cx_double b_ = -alpha * dz / 2 / 1i / n0 / k0;
			//求解CN差分方程
			CNsolve(
				a(j), b(j),
				Ay_V(i), Ay_V(i),
				By_V(i), By_V(i),
				Ax_V(i), Ax_V(i),
				Bx_V(i), Bx_V(i),
				C_V(i), C_V(i),
				D_V(i), D_V(i),
				uin, vin,
				uout, vout
			);
			uin = uout;
			vin = vout;
		}
		Ex.col(i + 1) = uout;
		Ey.col(i + 1) = vout;
	}
	clock_t t2 = clock();
	cout << "\t\t\t" << (double)(t2 - t1) / CLOCKS_PER_SEC << "s" << endl;

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