#include <iostream>
#include "hyp.h"

using namespace mpfr;
using namespace std;

int main() {

  unsigned int prec = 9999;
  mpreal::set_default_prec(prec);
  try {
    cout.precision(128);
    complexd a(3.3, 0.0);
    complexd b(3.3 + 24.5, 0.0);
    complexd z(2513.53 * 0.301, 0.0);
    complexd result = conhyp(a, b, z);
    cout << "1F1(" << a << ", " << b << ", " << z << ") = ";
    cout << result.real() << " + " << result.imag() << "i" << endl;
    cout << "with " << mpfr::mpreal_min(prec) << " minimum and ";
    cout << mpfr::mpreal_max(prec) << " maximum precision." << endl;
    return 0;
  } catch (const char *errmsg) {
    std::cerr << "error: " << errmsg << endl;
    return 1;
  }
}
