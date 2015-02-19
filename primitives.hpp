
#ifndef PRIMITIVES_HPP
#define PRIMITIVES_HPP

template <size_t dim,
          class T,
          class Real = double>
class PointBase
{
  Real operator[](const size_t& i) const
  {
    return static_cast<T*>(this)[i];
  }
};


template <size_t dim,
	  class Real = double>
class Point : PointBase<dim, Point<dim> >
{
protected:
  Real coords[dim];

public:
  Point() : coords() { }

  Real operator[](const size_t& i) const
  {
    return coords[i];
  }
};

#endif
