#include "load_hdf5.h"


void loadData(map<string, cube>& dev,
	string filepath,
	field<string>& s) {
	// map<string,cube> device;
	// device["z"]=aload()
	cout << "加载数据	：\n" ;
	cube tmp;

	for (int i = 0; i < s.size(); i++) {
		tmp.load(hdf5_name(filepath, s(i)));
		 cout <<"\t"<< s(i) << "\t\t" <<
			 tmp.n_rows << " " << tmp.n_cols << " " << tmp.n_slices << endl;
		dev[s(i)] = tmp;
		// tmp.print();
	}
	//dev["x"].print();
}