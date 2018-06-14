%module NTPolySwig
%include "typemaps.i"
%apply double& OUTPUT { double& chemical_potential_out };
%apply double *OUTPUT { double *max_power_eig };
%apply double *OUTPUT { double *min_ger_eig };
%apply double *OUTPUT { double *max_ger_eig };
%{
#include "SolverBase.h"
#include "ChebyshevSolvers.h"
#include "DensityMatrixSolvers.h"
#include "DistributedSparseMatrix.h"
#include "DistributedMatrixMemoryPool.h"
#include "EigenBounds.h"
#include "EigenSolvers.h"
#include "ExponentialSolvers.h"
#include "FixedSolversParameters.h"
#include "GeometryOptimization.h"
#include "HermiteSolvers.h"
#include "InverseSolvers.h"
#include "IterativeSolversParameters.h"
#include "LinearSolvers.h"
#include "LoadBalancer.h"
#include "MatrixMemoryPool.h"
#include "MinimizerSolvers.h"
#include "Permutation.h"
#include "Polynomial.h"
#include "ProcessGrid.h"
#include "RootSolvers.h"
#include "SignSolvers.h"
#include "SparseMatrix.h"
#include "SquareRootSolvers.h"
#include "TrigonometrySolvers.h"
#include "Triplet.h"
#include "TripletList.h"
using namespace NTPoly;
%}
%include "std_string.i"

%include "SolverBase.h"
%include "ChebyshevSolvers.h"
%include "DensityMatrixSolvers.h"
%include "DistributedSparseMatrix.h"
%include "DistributedMatrixMemoryPool.h"
%include "EigenBounds.h"
%include "EigenSolvers.h"
%include "ExponentialSolvers.h"
%include "FixedSolversParameters.h"
%include "GeometryOptimization.h"
%include "HermiteSolvers.h"
%include "InverseSolvers.h"
%include "IterativeSolversParameters.h"
%include "LinearSolvers.h"
%include "LoadBalancer.h"
%include "MatrixMemoryPool.h"
%include "MinimizerSolvers.h"
%include "Permutation.h"
%include "Polynomial.h"
%include "ProcessGrid.h"
%include "RootSolvers.h"
%include "SignSolvers.h"
%include "SparseMatrix.h"
%include "SquareRootSolvers.h"
%include "TrigonometrySolvers.h"
%include "Triplet.h"
%include "TripletList.h"

%template(MatrixMemoryPool_r) NTPoly::MatrixMemoryPool<double>;
%template(SparseMatrix_r) NTPoly::SparseMatrix<double>;
%template(TripletList_r) NTPoly::TripletList<double>;
%template(Triplet_r) NTPoly::Triplet<double>;

%template(MatrixMemoryPool_c) NTPoly::MatrixMemoryPool<double _Complex>;
%template(SparseMatrix_c) NTPoly::SparseMatrix<double _Complex>;
%template(TripletList_c) NTPoly::TripletList<double _Complex>;
%template(Triplet_c) NTPoly::Triplet<double _Complex>;
