#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

void generate(std::string fname, size_t K, size_t n) {
  
	std::ofstream file;

  file.open (fname);
  file << "%%MatrixMarket matrix coordinate real diagonal" << std::endl;
  file << n << " " << n << " " <<  n << std::endl;
	for (size_t i = 0; i < n; ++i) {
		file << i + 1 << " " << i + 1 << " " << std::scientific << std::exp(std::log(1./K)*i/(n - 1)) << std::endl;
	}
	
	file.close();	
}

void generate2(std::string fname, size_t K, size_t n) {
  
	std::ofstream file;

  file.open (fname);
  file << "%%MatrixMarket matrix coordinate real diagonal" << std::endl;
  file << n << " " << n << " " <<  n << std::endl;
	for (size_t i = 0; i < n; ++i) {
		file << i + 1 << " " << i + 1 << " " << std::scientific << std::log(std::exp(1. / K) + (std::exp(1) - std::exp(1./K))*i/(n - 1)) << std::endl;
	}
	
	file.close();	
}

int main() {
	
  size_t n = 10000;

	generate("dmat1.mtx", 1e+05, n);
	generate("dmat2.mtx", 1e+10, n);
	generate("dmat3.mtx", 1e+15, n);

	generate2("dmat4.mtx", 1e+05, n);
	generate2("dmat5.mtx", 1e+10, n);
	generate2("dmat6.mtx", 1e+15, n);
	
	return 0;
}