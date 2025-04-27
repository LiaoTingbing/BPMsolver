#include "bpm_class.h"

void Bpm::init()
{
 
	nx_ = (*dev_)["x"].n_elem;
	ny_ = (*dev_)["y"].n_elem;
	nz_ = (*dev_)["z"].n_elem;
	nt_ = nx_ * ny_;

	cout << "nx = " << nx_ << endl;
	cout << "ny = " << ny_ << endl;
	cout << "nz = " << nz_ << endl;
	cout << "nt = " << nt_ << endl;


	dx_ = (*dev_)["x"](1) - (*dev_)["x"](0);
	dy_ = (*dev_)["y"](1) - (*dev_)["y"](0);
	dz_ = (*dev_)["z"](1) - (*dev_)["z"](0);
	lambda_ = (*dev_)["lambda"](0);
	k0_ = 2 * PI / lambda_;
	n0_ = (*dev_)["neff"](0);
	alpha_ = 0.5; // CN 差分控制参数 

	cout << "dx = " << dx_ << endl;
	cout << "dy = " << dy_ << endl;
	cout << "dz = " << dz_ << endl;
	cout << "lambda = " << lambda_ << endl;
	//cout << "k0 = " << k0 << endl;
	cout << "n0 = " << n0_ << endl;
	cout << "alpha = " << alpha_ << endl;

	ex_.set_size(nt_, nz_);
	ey_.set_size(nt_, nz_);

}

void Bpm::computePML(int layersPML)
{
	//int layersPML = 10 ;
	cout << "\t初始化PML参数\tPML层数："<<layersPML<<"\n";
	cx_mat sx(nx_, ny_, fill::ones);
	cx_mat sy(nx_, ny_, fill::ones);
	int m = 3;
	double omega = C0 * k0_;
	const double R0 = 1e-6;
	double sigmaXMax = -(m + 1) * log10(R0) / (2 * 373 * dx_ * layersPML);
	double sigmaYMax = -(m + 1) * log10(R0) / (2 * 373 * dy_ * layersPML);


	cx_double value ;
	for (int i = 0;i < layersPML + 1;i++) {

		value = 1.0 + sigmaXMax * pow((layersPML - i + 0.0) / layersPML, m) / (IU * omega * EPS0);
		sx.row(i) = cx_rowvec(ny_, fill::value(value));
		sx.row(nx_ - 1 - i) = sx.row(i);

		value = 1.0 + sigmaYMax * pow((layersPML - i + 0.0) / layersPML, m) / (IU * omega * EPS0);
		sy.col(i) = cx_vec(nx_, fill::value(value));
		sy.col(ny_ - 1 - i) = sy.col(i);

	}
	sx_ = sx.as_col();
	sy_ = sy.as_col();
	isx_ = 1 / sx_;
	isy_ = 1 / sy_;
}

void Bpm::computeMatrix()
{
	cout << "\t计算Ax Ay Bx By C D";

	Ax_.set_size(nz_);
	Ay_.set_size(nz_);
	Bx_.set_size(nz_);
	By_.set_size(nz_);
	C_.set_size(nz_);
	D_.set_size(nz_);
	clock_t t1 = clock();
#pragma omp parallel for
	for (int i = 0; i < nz_;i++)
	{
		cx_vec erx = (*dev_)["Epsx"].slice(i).as_col() + 0.0 * IU;
		cx_vec erz = (*dev_)["Epsx"].slice(i).as_col() + 0.0 * IU;
		cx_vec ery = (*dev_)["Epsy"].slice(i).as_col() + 0.0 * IU;
		cx_vec erxy = (*dev_)["Epsxy"].slice(i).as_col() + 0.0 * IU;
		cx_vec eryx = (*dev_)["Epsyx"].slice(i).as_col() + 0.0 * IU;


		dxdxFunc(isx_, sx_ % erz, erx, nx_, ny_, dx_, dy_, Ax_(i));
		Ax_(i).diagArr.col(1) += 0.5 * k0_ * k0_ * (erx - n0_ * n0_);

		dydyFunc(isy_, sy_, cx_vec(nx_ * ny_, fill::ones), nx_, ny_, dx_, dy_, Ay_(i));
		Ay_(i).diagArr.col(1) += 0.5 * k0_ * k0_ * (erx - n0_ * n0_);

		dxdxFunc(isx_, sx_, cx_vec(nx_ * ny_, fill::ones), nx_, ny_, dx_, dy_, Bx_(i));
		Bx_(i).diagArr.col(1) += 0.5 * k0_ * k0_ * (ery - n0_ * n0_);

		dydyFunc(isy_, sy_ % erz, ery, nx_, ny_, dx_, dy_, By_(i));
		By_(i).diagArr.col(1) += 0.5 * k0_ * k0_ * (ery - n0_ * n0_);

		DiagStruct c1, c2;
		dxdyFunc(isx_, sy_ % erz, ery, nx_, ny_, dx_, dy_, c1);
		dxdyFunc(isx_, sy_, cx_vec(nx_ * ny_, fill::ones), nx_, ny_, dx_, dy_, c2);
		C_(i).diagArr = c1.diagArr - c2.diagArr;
		C_(i).diagArr.col(2) = erxy * k0_ * k0_;
		C_(i).pos = c1.pos;


		DiagStruct c3, c4;
		dydxFunc(isy_, sx_ % erz, erx, nx_, ny_, dx_, dy_, c3);
		dydxFunc(isy_, sx_, cx_vec(nx_ * ny_, fill::ones), nx_, ny_, dx_, dy_, c4);
		D_(i).diagArr = c3.diagArr - c4.diagArr;
		D_(i).diagArr.col(2) = erxy * k0_ * k0_;
		D_(i).pos = c3.pos;
 
	}
	clock_t t2 = clock();
	cout << "\t:\t" << (double)(t2 - t1) / CLOCKS_PER_SEC << "s" << endl;


}










void Bpm::qusiTmPropagate()
{



}

void Bpm::qusiTePropagate()
{

}


void Bpm::fullVectorPropagateSimple()
{

	cout << "\t全矢量传播\n";

	ex_.col(0) = vectorise((*dev_)["Exin"]) + 0.0 * IU;
	ey_.col(0) = vectorise((*dev_)["Eyin"]) + 0.0 * IU;

	clock_t t1 = clock();

	for (int i = 0; i < nz_ - 1; i++) {
		cout.flush();
		cout << "\r\t" << i + 1 << "/" << nz_ - 1 << "层";

		//求解CN差分方程
		cx_double a = (1 - alpha_) * dz_ / 2 / 1i / n0_ / k0_;
		cx_double b = -alpha_ * dz_ / 2 / 1i / n0_ / k0_;
		cx_vec uout, vout;
		cnSolve(
			a, b,
			Ay_(i), Ay_(i + 1),
			By_(i), By_(i + 1),
			Ax_(i), Ax_(i + 1),
			Bx_(i), Bx_(i + 1),
			C_(i), C_(i + 1),
			D_(i), D_(i + 1),
			ex_.col(i), ey_.col(i),
			uout, vout
		);
		ex_.col(i + 1) = uout;
		ey_.col(i + 1) = vout;
	}
	clock_t t2 = clock();
	cout << "\t\t\t" << (double)(t2 - t1) / CLOCKS_PER_SEC << "s" << endl;

}

void Bpm::fullVectorWideAnglePropagateSimple(int order)
{
	cout << "\t全矢量传播\t";
	switch (order)
	{
	case 0:
		cout << "广角阶数\t(1,0)\n";
		break;
	case 1:
		cout << "广角阶数\t(1,1)\n";
		break;
	case 2:
		cout << "广角阶数\t(2,2)\n";
		break;
	case 3:
		cout << "广角阶数\t(3,3)\n";
		break;
	case 4:
		cout << "广角阶数\t(4,4)\n";
		break;
	case 5:
		cout << "广角阶数\t(5,5)\n";
		break;
	default:
		break;
	}

	ex_.col(0) = vectorise((*dev_)["Exin"]) + 0.0 * IU;
	ey_.col(0) = vectorise((*dev_)["Eyin"]) + 0.0 * IU;

	clock_t t1 = clock();
	cx_vec coffUp = NN(order) - IU * n0_ * k0_ * dz_ * (1 - alpha_) * MN(order);
	cx_vec coffDown = NN(order) + IU * n0_ * k0_ * dz_ * alpha_ * MN(order);
	cx_vec a = -1.0 / roots(coffUp) / n0_ / n0_ / k0_ / k0_;
	cx_vec b = -1.0 / roots(coffDown) / n0_ / n0_ / k0_ / k0_;
	cx_vec uout, vout, uin, vin;

	// 第i层传播
	for (int i = 0; i < nz_ - 1; i++) {
		cout.flush();
		cout << "\r\t" << i + 1 << "/" << nz_ - 1 << "层";

		// dz放里面可能不均匀的dz和变化的n0

		//广角循环
		uin = ex_.col(i) ;
		vin = ey_.col(i) ;
		for (int j = 0;j < NN(order).size()-1;j++) {
			//cx_double a_ = (1 - alpha) * dz / 2 / 1i / n0 / k0;
			//cx_double b_ = -alpha * dz / 2 / 1i / n0 / k0;
			//求解CN差分方程
			cnSolve(
				a(j), b(j),
				Ay_(i), Ay_(i + 1),
				By_(i), By_(i + 1),
				Ax_(i), Ax_(i + 1),
				Bx_(i), Bx_(i + 1),
				C_(i), C_(i + 1),
				D_(i), D_(i + 1),
				uin, vin,
				uout, vout
			);
			uin = uout;
			vin = vout;
		}
		ex_.col(i + 1) = uout;
		ey_.col(i + 1) = vout;
	}
	clock_t t2 = clock();
	cout << "\t\t\t" << (double)(t2 - t1) / CLOCKS_PER_SEC << "s" << endl;

}

void Bpm::postData()
{

	system("del output/bpmOut.h5");

	cout << "\n写入文件中：" << endl;

	string filePath = "output/bpmOut.h5";

	(*dev_)["x"].save(hdf5_name(filePath, "x"));
	(*dev_)["y"].save(hdf5_name(filePath, "y", hdf5_opts::append));
	(*dev_)["z"].save(hdf5_name(filePath, "z", hdf5_opts::append));

	mat exAbs = abs(ex_);
	exAbs.save(hdf5_name(filePath, "Ex_abs", hdf5_opts::append));
	mat eyAbs = abs(ey_);
	eyAbs.save(hdf5_name(filePath, "Ey_abs", hdf5_opts::append));

	mat exReal = real(ex_);
	exReal.save(hdf5_name(filePath, "Ex_real", hdf5_opts::append));

	mat eyReal = real(ey_);
	eyReal.save(hdf5_name(filePath, "Ey_real", hdf5_opts::append));

}