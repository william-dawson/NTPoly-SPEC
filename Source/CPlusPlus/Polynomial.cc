#include "Polynomial.h"
#include "SolverParameters.h"
using namespace NTPoly;

////////////////////////////////////////////////////////////////////////////////
extern "C" {
#include "Polynomial_c.h"
}

////////////////////////////////////////////////////////////////////////////////
namespace NTPoly {
//////////////////////////////////////////////////////////////////////////////
Polynomial::Polynomial(int degree) {
  ConstructPolynomial_wrp(this->ih_this, &degree);
}
//////////////////////////////////////////////////////////////////////////////
Polynomial::~Polynomial() { DestructPolynomial_wrp(this->ih_this); }
//////////////////////////////////////////////////////////////////////////////
void Polynomial::SetCoefficient(int degree, double coefficient) {
  int temp_degree = degree + 1;
  SetCoefficient_wrp(this->ih_this, &temp_degree, &coefficient);
}
//////////////////////////////////////////////////////////////////////////////
void Polynomial::HornerCompute(
    const Matrix_ps &InputMat, Matrix_ps &OutputMat,
    const SolverParameters &solver_parameters) const {
  HornerCompute_wrp(GetIH(InputMat), GetIH(OutputMat), this->ih_this,
                    GetIH(solver_parameters));
}
//////////////////////////////////////////////////////////////////////////////
void Polynomial::PatersonStockmeyerCompute(
    const Matrix_ps &InputMat, Matrix_ps &OutputMat,
    const SolverParameters &solver_parameters) const {
  PatersonStockmeyerCompute_wrp(GetIH(InputMat), GetIH(OutputMat),
                                this->ih_this, GetIH(solver_parameters));
}
} // namespace NTPoly
