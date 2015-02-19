
#include <iostream>

#include "primitives.hpp"


int main()
{

  Point<2> x;

  std::cout << x[0] << ", " << x[1] << std::endl;

  Point<2> y = {1.0, 2.0};

  std::cout << y[0] << ", " << y[1] << std::endl;

  return 0;
}
