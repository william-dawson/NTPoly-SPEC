''' @package testDistributedSparseMatrix
A test suite for the Distributed Sparse Matrix module.'''
import unittest
import NTPolySwig as nt

import warnings
warnings.filterwarnings(action="ignore", module="scipy",
                        message="^internal gelsd")

import scipy
from scipy.linalg import pinv, funm, polar, cholesky
from scipy.sparse import csr_matrix, csc_matrix, rand, identity
from scipy.io import mmread, mmwrite
from scipy.sparse.linalg import norm, inv, eigsh
from numpy import zeros, sqrt, power, \
    sign, exp, log, sin, cos, linspace, diag, dot, sort
from numpy.linalg import eigh, svd
from numpy.linalg import norm as densenorm
import numpy
from random import random
import os
from numpy.polynomial.chebyshev import chebfit, chebval
from numpy.polynomial.hermite import hermfit, hermval
from mpi4py import MPI
# MPI global communicator.
comm = MPI.COMM_WORLD

from Helpers import THRESHOLD
from Helpers import result_file
from Helpers import scratch_dir


class TestSolvers(unittest.TestCase):
    '''A test class for the different kinds of solvers.'''
    # First input file.
    input_file = scratch_dir + "/input.mtx"
    # Second input file.
    input_file2 = scratch_dir + "/input2.mtx"
    # Matrix to compare against.
    CheckMat = 0
    # Rank of the current process.
    my_rank = 0
    # Dimension of the matrices to test.
    matrix_dimension = 17

    @classmethod
    def setUpClass(self):
        '''Set up all of the tests.'''
        rows = int(os.environ['PROCESS_ROWS'])
        columns = int(os.environ['PROCESS_COLUMNS'])
        slices = int(os.environ['PROCESS_SLICES'])
        nt.ConstructProcessGrid(rows, columns, slices)

    def setUp(self):
        '''Set up all of the tests.'''
        # Rank of the current process.
        self.my_rank = comm.Get_rank()
        # Parameters for iterative solvers.
        self.iterative_solver_parameters = nt.IterativeSolverParameters()
        # Parameters for fixed solvers.
        self.fixed_solver_parameters = nt.FixedSolverParameters()
        self.fixed_solver_parameters.SetVerbosity(True)
        self.iterative_solver_parameters.SetVerbosity(True)

    def check_result(self):
        '''Compare two computed matrices.'''
        normval = 0
        relative_error = 0
        if (self.my_rank == 0):
            ResultMat = mmread(result_file)
            normval = abs(norm(self.CheckMat - ResultMat))
            relative_error = normval / norm(self.CheckMat)
            print("\nNorm:", normval)
            print("Relative_Error:", relative_error)
        global_norm = comm.bcast(normval, root=0)
        global_error = comm.bcast(relative_error, root=0)
        self.assertLessEqual(global_error, THRESHOLD)

    def check_diag(self):
        '''Compare two diagonal matrices.'''
        normval = 0
        relative_error = 0
        if (self.my_rank == 0):
            ResultMat = sort(diag(mmread(result_file).todense()))
            CheckDiag = sort(diag(self.CheckMat.todense()))
            normval = abs(densenorm(CheckDiag - ResultMat))
            relative_error = normval / norm(self.CheckMat)
            print("\nNorm:", normval)
            print("Relative_Error:", relative_error)
        global_norm = comm.bcast(normval, root=0)
        global_error = comm.bcast(relative_error, root=0)
        self.assertLessEqual(global_error, THRESHOLD)

    def test_invert(self):
        '''Test our ability to invert matrices.'''
        # Starting Matrix
        temp_mat = zeros((self.matrix_dimension, self.matrix_dimension))
        for j in range(0, self.matrix_dimension):
            for i in range(0, self.matrix_dimension):
                if i == j:
                    temp_mat[j, i] = 1.0
                else:
                    temp_mat[j, i] = random() / (float(i - j)**2)
        temp_mat = (temp_mat.T + temp_mat) * 0.5
        matrix1 = csr_matrix(temp_mat)

        # Check Matrix
        self.CheckMat = inv(csc_matrix(matrix1))
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        overlap_matrix = nt.DistributedSparseMatrix(
            self.input_file, False)
        inverse_matrix = nt.DistributedSparseMatrix(
            overlap_matrix.GetActualDimension())
        permutation = nt.Permutation(overlap_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.iterative_solver_parameters.SetLoadBalance(permutation)

        nt.InverseSolvers.Invert(overlap_matrix, inverse_matrix,
                                 self.iterative_solver_parameters)

        inverse_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_pseudoinverse(self):
        '''Test our ability to compute the pseudoinverse of matrices.'''
        # Starting Matrix.
        temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                        density=1.0)
        matrix1 = csr_matrix(temp_mat)
        # Make it rank deficient
        k = int(self.matrix_dimension / 2)
        matrix1 = matrix1[k:].dot(matrix1[k:].T)

        # Check Matrix
        self.CheckMat = csr_matrix(pinv(matrix1.todense()))
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        overlap_matrix = nt.DistributedSparseMatrix(
            self.input_file, False)
        inverse_matrix = nt.DistributedSparseMatrix(
            overlap_matrix.GetActualDimension())
        permutation = nt.Permutation(overlap_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.iterative_solver_parameters.SetLoadBalance(permutation)

        nt.InverseSolvers.PseudoInverse(overlap_matrix, inverse_matrix,
                                        self.iterative_solver_parameters)

        inverse_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_inversesquareroot(self):
        '''Test our ability to compute the inverse square root of matrices.'''
        # Starting Matrix. Care taken to make sure eigenvalues are positive.
        temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                        density=1.0)
        temp_mat = (temp_mat.T + temp_mat)
        identity_matrix = identity(self.matrix_dimension)
        temp_mat = temp_mat + identity_matrix * self.matrix_dimension
        matrix1 = csr_matrix(temp_mat)

        # Check Matrix
        dense_check = funm(temp_mat.todense(), lambda x: 1.0 / sqrt(x))
        self.CheckMat = csr_matrix(dense_check)
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        overlap_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        inverse_matrix = nt.DistributedSparseMatrix(
            overlap_matrix.GetActualDimension())
        permutation = nt.Permutation(overlap_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.iterative_solver_parameters.SetLoadBalance(permutation)
        nt.SquareRootSolvers.InverseSquareRoot(overlap_matrix, inverse_matrix,
                                               self.iterative_solver_parameters)
        inverse_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_squareroot(self):
        '''Test our ability to compute the square root of matrices.'''
        # Starting Matrix. Care taken to make sure eigenvalues are positive.
        temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                        density=1.0)
        temp_mat = (temp_mat.T + temp_mat)
        identity_matrix = identity(self.matrix_dimension)
        temp_mat = temp_mat + identity_matrix * self.matrix_dimension
        matrix1 = csr_matrix(temp_mat)

        # Check Matrix
        dense_check = funm(temp_mat.todense(), lambda x: sqrt(x))
        self.CheckMat = csr_matrix(dense_check)
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        root_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.iterative_solver_parameters.SetLoadBalance(permutation)
        nt.SquareRootSolvers.SquareRoot(input_matrix, root_matrix,
                                        self.iterative_solver_parameters)
        root_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_inverseroot(self):
        '''Test our ability to compute  general matrix inverse root.'''
        roots = [1, 2, 3, 4, 5, 6, 7, 8]
        for root in roots:
            print("Root:", root)
            # Starting Matrix. Care taken to make sure eigenvalues are
            # positive.
            temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                            density=1.0)
            temp_mat = (temp_mat.T + temp_mat)
            identity_matrix = identity(self.matrix_dimension)
            temp_mat = temp_mat + identity_matrix * self.matrix_dimension
            matrix1 = csr_matrix(temp_mat)

            # Check Matrix
            dense_check = funm(temp_mat.todense(),
                               lambda x: power(x, -1.0 / root))
            self.CheckMat = csr_matrix(dense_check)
            if self.my_rank == 0:
                mmwrite(self.input_file, csr_matrix(matrix1))
            comm.barrier()

            # Result Matrix
            input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
            inverse_matrix = nt.DistributedSparseMatrix(
                input_matrix.GetActualDimension())
            permutation = nt.Permutation(input_matrix.GetLogicalDimension())
            permutation.SetRandomPermutation()
            self.iterative_solver_parameters.SetLoadBalance(permutation)
            nt.RootSolvers.ComputeInverseRoot(input_matrix, inverse_matrix,
                                              root,
                                              self.iterative_solver_parameters)
            inverse_matrix.WriteToMatrixMarket(result_file)
            comm.barrier()

            self.check_result()

    def test_root(self):
        '''Test our ability to compute  general matrix root.'''
        roots = [1, 2, 3, 4, 5, 6, 7, 8]
        for root in roots:
            print("Root", root)
            # Starting Matrix. Care taken to make sure eigenvalues are
            # positive.
            temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                            density=1.0)
            temp_mat = (temp_mat.T + temp_mat)
            identity_matrix = identity(self.matrix_dimension)
            temp_mat = temp_mat + identity_matrix * self.matrix_dimension
            matrix1 = csr_matrix(temp_mat)

            # Check Matrix
            dense_check = funm(temp_mat.todense(),
                               lambda x: power(x, 1.0 / root))
            self.CheckMat = csr_matrix(dense_check)
            if self.my_rank == 0:
                mmwrite(self.input_file, csr_matrix(matrix1))
            comm.barrier()

            # Result Matrix
            input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
            inverse_matrix = nt.DistributedSparseMatrix(
                input_matrix.GetActualDimension())
            permutation = nt.Permutation(input_matrix.GetLogicalDimension())
            permutation.SetRandomPermutation()
            self.iterative_solver_parameters.SetLoadBalance(permutation)
            nt.RootSolvers.ComputeRoot(input_matrix, inverse_matrix,
                                       root, self.iterative_solver_parameters)
            inverse_matrix.WriteToMatrixMarket(result_file)
            comm.barrier()
            self.check_result()

    def test_signfunction(self):
        '''Test our ability to compute the matrix sign function.'''
        # Starting Matrix
        temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                        density=1.0)
        temp_mat = (temp_mat.T + temp_mat)
        matrix1 = csr_matrix(temp_mat)

        # Check Matrix
        dense_check = funm(temp_mat.todense(), lambda x: sign(x))
        self.CheckMat = csr_matrix(dense_check)
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        sign_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.iterative_solver_parameters.SetLoadBalance(permutation)
        nt.SignSolvers.ComputeSign(input_matrix, sign_matrix,
                                   self.iterative_solver_parameters)
        sign_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_exponentialfunction(self):
        '''Test our ability to compute the matrix exponential.'''
        # Starting Matrix
        temp_mat = rand(
            self.matrix_dimension, self.matrix_dimension, density=1.0)
        temp_mat = 8 * (1.0 / self.matrix_dimension) * (temp_mat.T + temp_mat)
        matrix1 = csr_matrix(temp_mat)

        # Check Matrix
        dense_check = funm(temp_mat.todense(), lambda x: exp(x))
        self.CheckMat = csr_matrix(dense_check)
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        exp_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.fixed_solver_parameters.SetLoadBalance(permutation)
        nt.ExponentialSolvers.ComputeExponential(input_matrix, exp_matrix,
                                                 self.fixed_solver_parameters)
        exp_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_exponentialpade(self):
        '''Test our ability to compute the matrix exponential using the
        pade method.'''
        # Starting Matrix
        temp_mat = rand(
            self.matrix_dimension, self.matrix_dimension, density=1.0)
        temp_mat = 8 * (1.0 / self.matrix_dimension) * (temp_mat.T + temp_mat)
        matrix1 = csr_matrix(temp_mat)

        # Check Matrix
        dense_check = funm(temp_mat.todense(), lambda x: exp(x))
        self.CheckMat = csr_matrix(dense_check)
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        exp_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.iterative_solver_parameters.SetLoadBalance(permutation)
        nt.ExponentialSolvers.ComputeExponentialPade(input_matrix, exp_matrix,
                                                     self.iterative_solver_parameters)
        exp_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_logarithmfunction(self):
        '''Test our ability to compute the matrix logarithm.'''
        # Starting Matrix. Care taken to make sure eigenvalues are positive.
        temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                        density=1.0)
        temp_mat = (temp_mat.T + temp_mat)
        identity_matrix = identity(self.matrix_dimension)
        temp_mat = 4 * (1.0 / self.matrix_dimension) * \
            (temp_mat + identity_matrix * self.matrix_dimension)
        matrix1 = csr_matrix(temp_mat)

        # Check Matrix
        dense_check = funm(temp_mat.todense(), lambda x: log(x))
        self.CheckMat = csr_matrix(dense_check)
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        log_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.fixed_solver_parameters.SetLoadBalance(permutation)
        nt.ExponentialSolvers.ComputeLogarithm(input_matrix, log_matrix,
                                               self.fixed_solver_parameters)
        log_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_exponentialround(self):
        '''Test our ability to compute the matrix exponential using a round
        trip calculation.'''
        temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                        density=1.0)
        temp_mat = 0.25 * (0.5 / self.matrix_dimension) * \
            (temp_mat.T + temp_mat)
        self.CheckMat = csr_matrix(temp_mat)

        # Check Matrix
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(temp_mat))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        exp_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        round_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        nt.ExponentialSolvers.ComputeExponential(input_matrix, exp_matrix,
                                                 self.fixed_solver_parameters)
        nt.ExponentialSolvers.ComputeLogarithm(exp_matrix, round_matrix,
                                               self.fixed_solver_parameters)
        round_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_sinfunction(self):
        '''Test our ability to compute the matrix sine.'''
        # Starting Matrix
        temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                        density=1.0)
        temp_mat = (1.0 / self.matrix_dimension) * (temp_mat.T + temp_mat)
        matrix1 = csr_matrix(temp_mat)

        # Check Matrix
        dense_check = funm(temp_mat.todense(), lambda x: sin(x))
        self.CheckMat = csr_matrix(dense_check)
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        sin_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.fixed_solver_parameters.SetLoadBalance(permutation)
        nt.TrigonometrySolvers.Sine(input_matrix, sin_matrix,
                                    self.fixed_solver_parameters)
        sin_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_cosfunction(self):
        '''Test our ability to compute the matrix cosine.'''
        # Starting Matrix
        temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                        density=1.0)
        temp_mat = (1.0 / self.matrix_dimension) * (temp_mat.T + temp_mat)
        matrix1 = csr_matrix(temp_mat)

        # Check Matrix
        dense_check = funm(temp_mat.todense(),
                           lambda x: cos(x))
        self.CheckMat = csr_matrix(dense_check)
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        cos_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.fixed_solver_parameters.SetLoadBalance(permutation)
        nt.TrigonometrySolvers.Cosine(input_matrix, cos_matrix,
                                      self.fixed_solver_parameters)
        cos_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_hornerfunction(self):
        '''Test our ability to compute a matrix polynomial using horner's
        method.'''
        # Coefficients of the polynomial
        coef = [1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625]

        # Starting Matrix
        temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                        density=1.0)
        temp_mat = (1.0 / self.matrix_dimension) * (temp_mat.T + temp_mat)
        matrix1 = csr_matrix(temp_mat)

        # Check Matrix
        val, vec = eigh(temp_mat.todense())
        for i in range(0, len(val)):
            temp = val[i]
            val[i] = 0
            for j in range(0, len(coef)):
                val[i] = val[i] + coef[j] * (temp**j)
        temp_poly = dot(dot(vec, diag(val)), vec.T)
        self.CheckMat = csr_matrix(temp_poly)
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        poly_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())

        polynomial = nt.Polynomial(len(coef))
        for j in range(0, len(coef)):
            polynomial.SetCoefficient(j, coef[j])

        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.fixed_solver_parameters.SetLoadBalance(permutation)
        polynomial.HornerCompute(input_matrix, poly_matrix,
                                 self.fixed_solver_parameters)
        poly_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_patersonstockmeyerfunction(self):
        '''Test our ability to compute a matrix polynomial using the paterson
        stockmeyer method.'''
        # Coefficients of the polynomial
        coef = [1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625]

        # Starting Matrix
        temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                        density=1.0)
        temp_mat = (1.0 / self.matrix_dimension) * (temp_mat.T + temp_mat)
        matrix1 = csr_matrix(temp_mat)

        # Check Matrix
        val, vec = eigh(temp_mat.todense())
        for i in range(0, len(val)):
            temp = val[i]
            val[i] = 0
            for j in range(0, len(coef)):
                val[i] = val[i] + coef[j] * (temp**j)
        temp_poly = dot(dot(vec, diag(val)), vec.T)
        self.CheckMat = csr_matrix(temp_poly)
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        poly_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())

        polynomial = nt.Polynomial(len(coef))
        for j in range(0, len(coef)):
            polynomial.SetCoefficient(j, coef[j])

        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.fixed_solver_parameters.SetLoadBalance(permutation)
        polynomial.PatersonStockmeyerCompute(input_matrix, poly_matrix,
                                             self.fixed_solver_parameters)
        poly_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_chebyshevfunction(self):
        '''Test our ability to compute using Chebyshev polynomials.'''
        # Starting Matrix
        temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                        density=1.0)
        temp_mat = (1.0 / self.matrix_dimension) * (temp_mat.T + temp_mat)
        matrix1 = csr_matrix(temp_mat)

        # Function
        x = linspace(-1.0, 1.0, 200)
        y = [cos(i) + sin(i) for i in x]
        coef = chebfit(x, y, 10)

        # Check Matrix
        dense_check = funm(temp_mat.todense(), lambda x: chebval(x, coef))
        self.CheckMat = csr_matrix(dense_check)
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        poly_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())

        polynomial = nt.ChebyshevPolynomial(len(coef))
        for j in range(0, len(coef)):
            polynomial.SetCoefficient(j, coef[j])

        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.fixed_solver_parameters.SetLoadBalance(permutation)
        polynomial.Compute(input_matrix, poly_matrix,
                           self.fixed_solver_parameters)
        poly_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_recursivechebyshevfunction(self):
        '''Test our ability to compute using Chebyshev polynomials
        recursively.'''
        # Starting Matrix
        temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                        density=1.0)
        temp_mat = (1.0 / self.matrix_dimension) * (temp_mat.T + temp_mat)
        matrix1 = csr_matrix(temp_mat)

        # Function
        x = linspace(-1.0, 1.0, 200)
        # y = [scipy.special.erfc(i) for i in x]
        y = [exp(i) for i in x]
        coef = chebfit(x, y, 16 - 1)

        # Check Matrix
        dense_check = funm(temp_mat.todense(), lambda x: chebval(x, coef))
        self.CheckMat = csr_matrix(dense_check)
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        poly_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())

        polynomial = nt.ChebyshevPolynomial(len(coef))
        for j in range(0, len(coef)):
            polynomial.SetCoefficient(j, coef[j])

        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.fixed_solver_parameters.SetLoadBalance(permutation)
        polynomial.ComputeFactorized(
            input_matrix, poly_matrix, self.fixed_solver_parameters)
        poly_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_cgsolve(self):
        '''Test our ability to solve general matrix equations with CG.'''
        # Starting Matrix
        A = rand(self.matrix_dimension, self.matrix_dimension, density=1.0)
        B = rand(self.matrix_dimension, self.matrix_dimension, density=1.0)
        A = A.T.dot(A)

        # Check Matrix
        Ainv = inv(A)
        self.CheckMat = csr_matrix(Ainv.dot(B))
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(A))
            mmwrite(self.input_file2, csr_matrix(B))
        comm.barrier()

        # Result Matrix
        AMat = nt.DistributedSparseMatrix(self.input_file, False)
        XMat = nt.DistributedSparseMatrix(AMat.GetActualDimension())
        BMat = nt.DistributedSparseMatrix(self.input_file2, False)
        permutation = nt.Permutation(AMat.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.iterative_solver_parameters.SetLoadBalance(permutation)

        nt.LinearSolvers.CGSolver(AMat, XMat, BMat,
                                  self.iterative_solver_parameters)

        XMat.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_powermethod(self):
        '''Test our ability to compute eigenvalues with the power method.'''
        # Starting Matrix
        temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                        density=1.0)
        temp_mat = (temp_mat.T + temp_mat)

        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(temp_mat))
        comm.barrier()

        # Result Matrix
        max_value = 0.0
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        max_value = nt.EigenBounds.PowerBounds(input_matrix,
                                               self.iterative_solver_parameters)
        comm.barrier()

        vals, vec = eigsh(temp_mat, which="LM", k=1)
        relative_error = abs(max_value - vals[0])
        global_error = comm.bcast(relative_error, root=0)
        self.assertLessEqual(global_error, THRESHOLD)

    def test_hermitefunction(self):
        '''Test our ability to compute using Hermite polynomials.'''
        # Starting Matrix
        temp_mat = rand(self.matrix_dimension,
                        self.matrix_dimension,
                        density=1.0)
        temp_mat = (1.0 / self.matrix_dimension) * (temp_mat)
        matrix1 = csr_matrix(temp_mat)

        # Function
        x = linspace(-1.0, 1.0, 200)
        y = [cos(i) + sin(i) for i in x]
        coef = hermfit(x, y, 10)

        # Check Matrix
        dense_check = funm(temp_mat.todense(), lambda x: hermval(x, coef))
        self.CheckMat = csr_matrix(dense_check)
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        poly_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())

        polynomial = nt.HermitePolynomial(len(coef))
        for j in range(0, len(coef)):
            polynomial.SetCoefficient(j, coef[j])

        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.fixed_solver_parameters.SetLoadBalance(permutation)
        polynomial.Compute(input_matrix, poly_matrix,
                           self.fixed_solver_parameters)
        poly_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_cholesky(self):
        '''Test subroutine that computes the cholesky decomposition.'''
        # Starting Matrix
        temp_mat = scipy.sparse.rand(self.matrix_dimension,
                                     self.matrix_dimension,
                                     density=1.0)
        temp_mat = (1.0 / self.matrix_dimension) * (temp_mat)
        # Symmetric Positive Definite
        matrix1 = csr_matrix(temp_mat.T.dot(temp_mat))

        # Check Matrix
        dense_check = cholesky(matrix1.todense(), lower=True)
        self.CheckMat = csr_matrix(dense_check)
        if self.my_rank == 0:
            scipy.io.mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)

        cholesky_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        nt.LinearSolvers.CholeskyDecomposition(input_matrix, cholesky_matrix,
                                               self.fixed_solver_parameters)

        cholesky_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_pivotedcholesky(self):
        '''Test subroutine that computes the pivoted cholesky decomposition.'''
        # Starting Matrix
        # temp_mat = scipy.sparse.rand(self.matrix_dimension,
        #                              self.matrix_dimension,
        #                              density=1.0, random_state=2)
        # temp_mat = (1.0 / self.matrix_dimension) * (temp_mat)
        # rand2 = 10*rand(self.matrix_dimension, self.matrix_dimension,
        #              density=0.1)
        # temp_mat = temp_mat + rand2 + rand2.T
        # Symmetric Positive Definite
        # matrix1 = csr_matrix(temp_mat.T.dot(temp_mat))
        matrix1 = mmread(os.environ["CholTest"])
        rank = 2

        self.CheckMat = csr_matrix(matrix1)
        if self.my_rank == 0:
            scipy.io.mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        A = nt.DistributedSparseMatrix(self.input_file, False)
        L = nt.DistributedSparseMatrix(A.GetActualDimension())
        LT = nt.DistributedSparseMatrix(A.GetActualDimension())
        LLT = nt.DistributedSparseMatrix(A.GetActualDimension())
        memory_pool = nt.DistributedMatrixMemoryPool()

        nt.LinearSolvers.PivotedCholeskyDecomposition(A, L, rank,
                                                      self.fixed_solver_parameters)
        LT.Transpose(L)
        LLT.Gemm(L, LT, memory_pool)

        LLT.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_polarfunction(self):
        '''Test our ability to compute the matrix polar decomposition.'''
        # Starting Matrix
        temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                        density=1.0)
        matrix1 = csr_matrix(temp_mat)

        # Check Matrix
        dense_check_u, dense_check_h = polar(matrix1.todense())
        self.CheckMat = csr_matrix(dense_check_h)
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        u_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        h_matrix = nt.DistributedSparseMatrix(
            input_matrix.GetActualDimension())
        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.iterative_solver_parameters.SetLoadBalance(permutation)
        nt.SignSolvers.ComputePolarDecomposition(input_matrix, u_matrix, h_matrix,
                                                 self.iterative_solver_parameters)
        h_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

        comm.barrier()
        self.CheckMat = csr_matrix(dense_check_u)
        u_matrix.WriteToMatrixMarket(result_file)

        self.check_result()

    def test_eigendecomposition(self):
        '''Test our ability to compute the eigen decomposition.'''
        # Starting Matrix
        temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                        density=1.0, random_state=2)
        matrix1 = csr_matrix(temp_mat + temp_mat.T)

        # Check Matrix
        vals, vecs = eigh(matrix1.todense())
        self.CheckMat = csr_matrix(diag(vals))
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.iterative_solver_parameters.SetLoadBalance(permutation)
        vec_matrix = nt.DistributedSparseMatrix(self.matrix_dimension)
        val_matrix = nt.DistributedSparseMatrix(self.matrix_dimension)
        nt.EigenSolvers.EigenDecomposition(input_matrix, vec_matrix, val_matrix,
                                           self.matrix_dimension,
                                           self.iterative_solver_parameters)
        val_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_eigendecompositionhalf(self):
        '''Test our ability to compute the eigen decomposition.'''
        # Starting Matrix
        temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                        density=1.0)
        matrix1 = csr_matrix(temp_mat + temp_mat.T)

        # Check Matrix
        vals, vecs = eigh(matrix1.todense())
        valshalf = vals[:int(self.matrix_dimension/2)]
        self.CheckMat = csr_matrix(diag(valshalf))
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.iterative_solver_parameters.SetLoadBalance(permutation)
        vec_matrix = nt.DistributedSparseMatrix(self.matrix_dimension)
        val_matrix = nt.DistributedSparseMatrix(self.matrix_dimension)
        nt.EigenSolvers.EigenDecomposition(input_matrix, vec_matrix, val_matrix,
                                           int(self.matrix_dimension/2),
                                           self.iterative_solver_parameters)
        val_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_result()

    def test_svd(self):
        '''Test our ability to compute the eigen decomposition.'''
        # Starting Matrix
        temp_mat = rand(self.matrix_dimension, self.matrix_dimension,
                        density=1.0)
        matrix1 = csr_matrix(temp_mat + temp_mat.T)

        # Check Matrix
        u, s, vh = svd(matrix1.todense())
        self.CheckMat = csr_matrix(diag(s))
        if self.my_rank == 0:
            mmwrite(self.input_file, csr_matrix(matrix1))
        comm.barrier()

        # Result Matrix
        input_matrix = nt.DistributedSparseMatrix(self.input_file, False)
        permutation = nt.Permutation(input_matrix.GetLogicalDimension())
        permutation.SetRandomPermutation()
        self.iterative_solver_parameters.SetLoadBalance(permutation)
        left_matrix = nt.DistributedSparseMatrix(self.matrix_dimension)
        right_matrix = nt.DistributedSparseMatrix(self.matrix_dimension)
        val_matrix = nt.DistributedSparseMatrix(self.matrix_dimension)
        nt.EigenSolvers.SingularValueDecompostion(input_matrix, left_matrix,
                                                  right_matrix, val_matrix,
                                                  self.iterative_solver_parameters)
        val_matrix.WriteToMatrixMarket(result_file)
        comm.barrier()

        self.check_diag()


if __name__ == '__main__':
    unittest.main()
