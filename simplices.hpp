
#ifndef SIMPLICES_HPP
#define SIMPLICES_HPP

#include <vector>
#include <memory>
#include <array>

using std::size_t;

template <size_t d, size_t k, typename Real = double> class Simplices;
template <size_t d, size_t k, typename Real = double> class Simplex;

/**
 * This class represents k-dimensional simplices in d-dimensional space.
 * Simplices are hierarchical objects; each simplex of the collection
 * consists of the set of indices in the collection of k-1 - simplices of
 * its faces.
 */
template <size_t d, size_t k, typename Real>
class Simplices
{
public:
  Simplices() {};

  size_t size() const { return simplices.size(); }

  void add(const std::array<size_t, k+1>& simplex) {
    simplices.push_back(simplex);
  }

  void add(const std::array<Real, d>& x) {
    sub_simplices.add(x);
  }

  Simplices<d, k-1, Real>& sub() {
    return sub_simplices;
  }

  Simplex<d, k, Real> operator[](const size_t index) {
    return Simplex<d, k, Real>(*this, index);
  }

private:
  std::vector<std::array<size_t, k+1> > simplices;
  Simplices<d, k-1, Real> sub_simplices;

  friend class Simplex<d, k, Real>;
};


/**
 * A simplex object, which internally stores a reference to its parent
 * collection and its index in the collection.
 */
template <size_t d, size_t k, typename Real>
class Simplex
{
public:
  size_t operator[](const size_t i) const {
    return simplices->simplices[index][i];
  }

private:
  Simplex(const Simplices<d, k, Real>& _simplices, const size_t _index)
    : simplices(&_simplices), index(_index) {}
  const Simplices<d, k, Real>* simplices;
  const size_t index;

  friend class Simplices<d, k, Real>;
};


/**
 * Specialization of the Simplices/Simplex classes for k = 0, i.e., sets of
 * interned d-dimensional points.
 */
template <size_t d, typename Real>
class Simplices <d, 0, Real>
{
public:
  Simplices() {};

  size_t size() const { return X.size(); }

  void add(const std::array<Real, d>& x) {
    X.push_back(x);
  }

  Simplex<d, 0, Real> operator[](const size_t index) const {
    return Simplex<d, 0, Real>(*this, index);
  }

private:
  std::vector<std::array<Real, d> > X;

  friend class Simplex<d, 0, Real>;
};


template <size_t d, typename Real>
class Simplex <d, 0, Real>
{
public:
  Real operator[](const size_t i) const {
    return simplices->X[index][i];
  }

private:
  Simplex(const Simplices<d, 0, Real>& _simplices, const size_t _index)
    : simplices(&_simplices), index(_index) {}
  const Simplices<d, 0, Real>* simplices;
  const size_t index;

  friend class Simplices<d, 0, Real>;
};



#endif
