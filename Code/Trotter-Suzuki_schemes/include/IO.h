#ifndef _IO_H_
#define _IO_H_

#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "Scheme.h"
#include "NScheme.h"
#include "Minimization.h"
#include "NMinimization.h"

using Eigen::VectorXcd;
using Eigen::VectorXd;
using VectorXld = Eigen::Matrix<long double, Dynamic, 1>;
using VectorXcld = Eigen::Matrix<complex<long double>, Dynamic, 1>;
using VectorXQ = Eigen::Matrix<quad, Dynamic, 1>;
using VectorXcQ = Eigen::Matrix<complex<quad>, Dynamic, 1>;

using namespace std;

string trim(const string &str);

// Helper function to parse a boolean value
bool parse_bool(const string &str);

// Helper function to parse a real or complex value
template <typename Scalar>
Scalar parse_value(const string &str);

// Helper function to parse a list to extract vectors or arrays
template <typename List>
List parse_list(const string &str);

template <typename RealT>
void print_num_min(ostream &stream, const int &order, const int &no_cycles,
                   const array<RealT, 4> &eps1, const array<RealT, 4> &eps2, const vector<RealT> &wi,
                   const RealT step, const int n_iter, const array<RealT, 2> &Ls, const RealT lambda, const int steps);

template <typename RealT>
void print_num_min_origin(ostream &stream, const int &order, const int &no_cycles,
                          const array<RealT, 4> &eps1, const array<RealT, 4> &eps2, const vector<RealT> &wi,
                          const RealT step, const int n_iter, const array<RealT, 2> &Ls, const RealT lambda, const RealT ratio, const int steps);

template <typename RealT>
void print_sym_min(ostream &stream, const int &order, const int &no_cycles,
                   const array<RealT, 4> &eps1, const array<RealT, 4> &eps2, const vector<RealT> &wi,
                   const int n_iter, const array<RealT, 2> &Ls, const RealT lambda, const int steps);

template <typename Vec>
void print_num_min1(ostream &stream, NScheme<Vec> &scheme, const Vec &a_init, const Vec &b_init);

template <typename Vec>
void print_sym_min1(ostream &stream, Scheme<typename Eigen::NumTraits<typename Vec::Scalar>::Real> &scheme,
                    const Vec &a_init, const Vec &b_init, const Vec &a_vec, const Vec &b_vec);

// Helper functions for printing in the minimization routines
template <typename Scalar, typename RealT>
void print_num_find(ostream &stream, const int &order, const int &no_cycles,
                    const array<RealT, 4> &eps1, const array<RealT, 4> &eps2, const vector<RealT> &wi,
                    const RealT step, const int n_iter, const array<RealT, 2> &Ls, const RealT lambda, const int N, const int steps,
                    Scalar mu, const vector<Scalar> &mus, Scalar sigma, int sum_convs, const RealT tol);

template <typename Scalar, typename RealT>
void print_num_find_origin(ostream &stream, const int &order, const int &no_cycles,
                           const array<RealT, 4> &eps1, const array<RealT, 4> &eps2, const vector<RealT> &wi,
                           const RealT step, const int n_iter, const array<RealT, 2> &Ls, const RealT lambda, const RealT &ratio, const int N, const int steps,
                           Scalar mu, const vector<Scalar> &mus, Scalar sigma, int sum_convs, const RealT tol);

template <typename Scalar, typename RealT>
void print_sym_find(ostream &stream, const int &order, const int &no_cycles,
                    const array<RealT, 4> &eps1, const array<RealT, 4> &eps2, const vector<RealT> &wi,
                    const int n_iter, const array<RealT, 2> &Ls, const RealT lambda, const int N, const int steps,
                    Scalar mu, const vector<Scalar> &mus, Scalar sigma, int sum_convs, const RealT tol);

template <typename RealT>
void print_find1(ostream &stream, const int &tracker, const double &ratio,
                 const RealT &eff2, const RealT &eff4, const RealT &eff6, const RealT &eff8);

template <typename Vec>
void print_find2(ostream &stream, const Vec &a_vec, const Vec &b_vec);

template <typename Vec>
void print_sym_scheme(ostream &stream, Scheme<double> &scheme, const Vec &a_vec, const Vec &b_vec);

template <typename Vec>
void print_sym_scheme(ostream &stream, Scheme<long double> &scheme, const Vec &a_vec, const Vec &b_vec);

template <typename Vec>
void print_num_scheme(ostream &stream, NScheme<Vec> &scheme);


// Class responsible for reading the input file
class IO
{
public:
  // Read file stream
  ifstream rfile;

  int routine;         // Which routine to use (Minimization: 0, Scheme: 1)
  int mode;            // Which mode to use (Numeric: 0, Symbolic: 1)
  string scalar_type;  // Scalar type to use (double, complex<double>, long_double, complex<long_double>)
  string save_dir{""}; // File to save to (if not provided no writing occurs)

  // General parameters read from the read_file
  int order;
  int no_cycles;

  // IO constructor
  IO(const string read_file);

  // IO deconstructor
  ~IO();

  // *** Routines ***

  // Numeric Minimization
  // Methods: - minimize: minimize the computed scheme manifold once using initial vector a_vec, b_vec
  //          - min_twostep: minimize the computed scheme manifold twice (constraint minimization in the 2. step)
  //          - find: find as many minima for the computed scheme manifold
  template <typename Scalar>
  void num_minim();

  // Numeric Scheme
  // Computes the errors and efficiencies of the scheme for the provided a_vec and b_vec
  template <typename Scalar>
  void num_scheme();


  // Symbolic Minimization
  // Methods: - minimize: minimize the computed/loaded scheme manifold once using initial vector a_vec, b_vec
  //          - min_twostep: minimize the computed/loaded scheme manifold twice (constraint minimization in the 2. step)
  //          - find: find as many minima for the computed/loaded scheme manifold
  template <typename Scalar>
  void sym_minim();

  // Symbolic Scheme
  // Methods: - scratch: compute scheme from scratch
  //          - load: load scheme from file
  //          - load_iterate: load scheme from file and iterate for 2 cycles
  // Optional methods (Execute if provided): - Compute errors and efficiencies (needs a_vec and b_vec)
  //                                         - Save the scheme (needs a save folder)
  template <typename Scalar>
  void sym_scheme();
};

#endif // _IO_H_
