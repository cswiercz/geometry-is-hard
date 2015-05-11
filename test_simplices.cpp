
#include <iostream>

#include "simplices.hpp"



int main(int argc, char **argv) {

  Simplices<2, 1> S;
  std::array<double, 2> p = {1.0, 2.0};

  S.add(p);

  Simplices<2, 0>& T = S.sub();

  std::cout << S.size() << std::endl;
  std::cout << T.size() << std::endl;

  Simplex<2, 0> s = T[0];
  std::cout << s[0] << ", " << s[1] << std::endl;

  return 0;
}
