!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> The driver program for the spec benchmark.
!! This program computes the inverse of a matrix using hotelling's method.
PROGRAM PremadeMatrixProgram
  USE DataTypesModule, ONLY : NTREAL
  USE InverseSolversModule, ONLY : Invert
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, WriteHeader
  USE PermutationModule, ONLY : Permutation_t, ConstructRandomPermutation
  USE ProcessGridModule, ONLY : ConstructProcessGrid, IsRoot, &
       & DestructProcessGrid
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructMatrixFromMatrixMarket, &
       & DestructMatrix, ConstructEmptyMatrix, FillMatrixIdentity
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, MatrixNorm, IncrementMatrix
  USE SolverParametersModule, ONLY : SolverParameters_t
  USE SquareRootSolversModule, ONLY : InverseSquareRoot
  USE MPI
  IMPLICIT NONE
  !! Variables for handling input parameters.
  CHARACTER(len=80) :: input_file
  INTEGER :: process_rows, process_columns, process_slices
  REAL(NTREAL) :: threshold
  REAL(NTREAL) :: converge
  TYPE(SolverParameters_t) :: solver_parameters
  TYPE(Permutation_t) :: permutation
  !! Matrices
  TYPE(Matrix_ps) :: Input
  TYPE(Matrix_ps) :: Identity
  TYPE(Matrix_ps) :: Result
  TYPE(Matrix_ps) :: Reference
  TYPE(Matrix_ps) :: Temp
  !! Temporary Variables
  CHARACTER(len=80) :: argument
  CHARACTER(len=80) :: argument_value
  INTEGER :: counter
  INTEGER :: provided, ierr
  INTEGER :: rank
  INTEGER :: II, loop_times
  REAL(NTREAL) :: error

  !! Setup MPI
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,rank, ierr)

  !! Process the input parameters.
  DO counter=1,command_argument_count(),2
     CALL get_command_argument(counter,argument)
     CALL get_command_argument(counter+1,argument_value)
     SELECT CASE(argument)
     CASE('--input')
        input_file = argument_value
     CASE('--process_rows')
        READ(argument_value,*) process_rows
     CASE('--process_columns')
        READ(argument_value,*) process_columns
     CASE('--process_slices')
        READ(argument_value,*) process_slices
     CASE('--threshold')
        READ(argument_value,*) threshold
     CASE('--converge')
        READ(argument_value,*) converge
     CASE('--loop_times')
        READ(argument_value,*) loop_times
     END SELECT
  END DO

  !! Setup the process grid.
  CALL ConstructProcessGrid(MPI_COMM_WORLD, process_rows, process_columns, &
       & process_slices)

  CALL WriteHeader("Command Line Parameters")
  CALL EnterSubLog
  CALL WriteElement(key="input", value=input_file)
  CALL WriteElement(key="process_rows", value=process_rows)
  CALL WriteElement(key="process_columns", value=process_columns)
  CALL WriteElement(key="process_slices", value=process_slices)
  CALL WriteElement(key="threshold", value=threshold)
  CALL WriteElement(key="converge", value=converge)
  CALL WriteElement(key="loop_times", value=loop_times)
  CALL ExitSubLog

  !! Read in the matrices from file.
  CALL ConstructMatrixFromMatrixMarket(Input, input_file)
  CALL ConstructEmptyMatrix(Identity, Input)
  CALL FillMatrixIdentity(Identity)

  !! Set Up The Solver Parameters.
  CALL ConstructRandomPermutation(permutation, Input%logical_matrix_dimension)
  solver_parameters = SolverParameters_t(&
       & converge_diff_in=converge, threshold_in=threshold, &
       & BalancePermutation_in=permutation, be_verbose_in=.TRUE.)

  !! Call the solver routine.
  !! Time this part.
  DO II = 1, loop_times
     CALL InverseSquareRoot(Input, Result, solver_parameters)
  END DO

  !! Check The Answer
  CALL MatrixMultiply(Input, Result, Temp, threshold_in=threshold)
  CALL MatrixMultiply(Result, Temp, Reference, threshold_in=threshold)
  CALL IncrementMatrix(Identity, Reference, alpha_in=-1.0_NTREAL)
  error = MatrixNorm(Reference)
  CALL WriteElement(key="error", value=error)
  IF (error .GT. converge*10) THEN
     CALL MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  END IF

  !! Cleanup
  CALL DestructMatrix(Input)
  CALL DestructMatrix(Reference)
  CALL DestructMatrix(Result)
  CALL DestructProcessGrid
  CALL MPI_Finalize(ierr)
END PROGRAM PremadeMatrixProgram
