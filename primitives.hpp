
#ifndef PRIMITIVES_HPP
#define PRIMITIVES_HPP

template <int dim,
          class T,
          class Real = double>
class PointBase
{
  Real operator[](const unsigned int &i) const
  {
    return static_cast<T*>(this)[i];
  }
};


template <int dim,
	  class Real = double>
class Point : PointBase<dim, Point<dim> >
{
protected:
  Real coords[dim];

public:
  Point() {
    for (unsigned int i = 0; i < dim; i++) {
      coords[i] = 0.0;
    }
  }

  Real operator[](const unsigned int &i) const
  {
    return coords[i];
  }
};

#endif
