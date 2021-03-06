!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Geometry Optimization
MODULE GeometryOptimizationModule
  USE DataTypesModule, ONLY : NTREAL
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteHeader, &
       & WriteElement, WriteListElement, WriteCitation
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, MatrixNorm, &
       & IncrementMatrix, MatrixTrace, ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, DestructMatrix, ConstructEmptyMatrix, &
       & FillMatrixIdentity, PrintMatrixInformation, CopyMatrix
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters
  USE SquareRootSolversModule, ONLY : SquareRoot, InverseSquareRoot
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: PurificationExtrapolate
  PUBLIC :: LowdinExtrapolate
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a new guess at the Density Matrix after updating the geometry.
  !> Based on the purification algorithm in \cite niklasson2010trace .
  SUBROUTINE PurificationExtrapolate(PreviousDensity, Overlap, nel, NewDensity,&
       & solver_parameters_in)
    !> Previous density to extrapolate from.
    TYPE(Matrix_ps), INTENT(IN) :: PreviousDensity
    !> The overlap matrix of the new geometry.
    TYPE(Matrix_ps), INTENT(IN) :: Overlap
    !> The number of electrons.
    INTEGER, INTENT(IN) :: nel
    !> The extrapolated density.
    TYPE(Matrix_ps), INTENT(INOUT) :: NewDensity
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: WorkingDensity
    TYPE(Matrix_ps) :: WorkingOverlap
    TYPE(Matrix_ps) :: AddBranch, SubtractBranch
    TYPE(Matrix_ps) :: TempMat1, TempMat2
    TYPE(Matrix_ps) :: Identity
    !! Local Variables
    REAL(NTREAL) :: subtract_trace
    REAL(NTREAL) :: add_trace
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: norm_value
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool1
    INTEGER :: outer_counter
    INTEGER :: total_iterations

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Density Matrix Extrapolator")
       CALL EnterSubLog
       CALL WriteElement(key="Method", value="Purification")
       CALL WriteCitation("niklasson2010trace")
       CALL PrintParameters(solver_parameters)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(NewDensity, PreviousDensity)
    CALL ConstructEmptyMatrix(WorkingDensity, PreviousDensity)
    CALL ConstructEmptyMatrix(WorkingOverlap, PreviousDensity)
    CALL ConstructEmptyMatrix(Identity, PreviousDensity)
    CALL FillMatrixIdentity(Identity)

    !! Compute the working hamiltonian.
    CALL CopyMatrix(PreviousDensity, WorkingDensity)
    CALL CopyMatrix(Overlap, WorkingOverlap)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(WorkingDensity, WorkingDensity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
       CALL PermuteMatrix(WorkingOverlap, WorkingOverlap, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF

    !! Finish Setup
    CALL CopyMatrix(WorkingDensity, NewDensity)
    CALL CopyMatrix(WorkingDensity, AddBranch)
    CALL CopyMatrix(WorkingDensity, SubtractBranch)

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       !! Figure Out Sigma Value. After which, XnS is stored in TempMat
       IF (outer_counter .GT. 1) THEN
          CALL MatrixMultiply(AddBranch, WorkingOverlap, TempMat1, &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
          CALL MatrixTrace(TempMat1, add_trace)
          CALL MatrixMultiply(SubtractBranch, WorkingOverlap, TempMat2, &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
          CALL MatrixTrace(TempMat2, subtract_trace)
          IF (ABS(nel - add_trace) .GT. ABS(nel - subtract_trace)) THEN
             !! Subtract Branch
             trace_value = subtract_trace
             CALL IncrementMatrix(AddBranch, NewDensity, -1.0_NTREAL)
             norm_value = MatrixNorm(NewDensity)
             CALL CopyMatrix(AddBranch, NewDensity)
             CALL CopyMatrix(TempMat2, TempMat1)
          ELSE
             !! Add Branch
             trace_value = add_trace
             CALL IncrementMatrix(SubtractBranch, NewDensity, -1.0_NTREAL)
             norm_value = MatrixNorm(NewDensity)
             CALL CopyMatrix(SubtractBranch, NewDensity)
          END IF
       ELSE
          CALL MatrixMultiply(NewDensity, WorkingOverlap, TempMat1, &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       END IF

       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Round", value=outer_counter-1)
          CALL EnterSubLog
          CALL WriteElement(key="Convergence", value=norm_value)
          CALL WriteElement(key="Trace", value=trace_value)
          CALL WriteElement(key="AddTrace", value=add_trace)
          CALL WriteElement(key="SubtractTrace", value=subtract_trace)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF

       !! Compute (I - XnS)Xn
       CALL IncrementMatrix(Identity, TempMat1, -1.0_NTREAL)
       CALL ScaleMatrix(TempMat1, -1.0_NTREAL)
       CALL MatrixMultiply(TempMat1, NewDensity, TempMat2, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

       !! Subtracted Version Xn - (I - XnS)Xn
       CALL CopyMatrix(NewDensity, SubtractBranch)
       CALL IncrementMatrix(TempMat2, SubtractBranch, -1.0_NTREAL)

       !! Added Version Xn + (I - XnS)Xn
       CALL CopyMatrix(NewDensity, AddBranch)
       CALL IncrementMatrix(TempMat2, AddBranch)

    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations", value=total_iterations)
       CALL PrintMatrixInformation(NewDensity)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(NewDensity, NewDensity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructMatrix(WorkingDensity)
    CALL DestructMatrix(WorkingOverlap)
    CALL DestructMatrix(Identity)
    CALL DestructMatrix(TempMat1)
    CALL DestructMatrix(TempMat2)
    CALL DestructMatrix(AddBranch)
    CALL DestructMatrix(SubtractBranch)
    CALL DestructMatrixMemoryPool(pool1)
    CALL DestructSolverParameters(solver_parameters)

  END SUBROUTINE PurificationExtrapolate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a new guess at the Density Matrix after updating the geometry.
  !> Based on the lowdin algorithm in \cite exner2002comparison .
  SUBROUTINE LowdinExtrapolate(PreviousDensity, OldOverlap, NewOverlap, &
       & NewDensity, solver_parameters_in)
    !> THe previous density to extrapolate from.
    TYPE(Matrix_ps), INTENT(IN)  :: PreviousDensity
    !> The old overlap matrix from the previous geometry.
    TYPE(Matrix_ps), INTENT(IN)  :: OldOverlap
    !> The new overlap matrix from the current geometry.
    TYPE(Matrix_ps), INTENT(IN)  :: NewOverlap
    !> The extrapolated density.
    TYPE(Matrix_ps), INTENT(INOUT) :: NewDensity
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: SQRMat
    TYPE(Matrix_ps) :: ISQMat
    TYPE(Matrix_ps) :: TempMat
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool1

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Density Matrix Extrapolator")
       CALL EnterSubLog
       CALL WriteElement(key="Method", value="Lowdin")
       CALL WriteCitation("exner2002comparison")
       CALL PrintParameters(solver_parameters)
    END IF

    CALL SquareRoot(OldOverlap, SQRMat, solver_parameters)
    CALL InverseSquareRoot(NewOverlap, ISQMat, solver_parameters)

    CALL MatrixMultiply(SQRMat, PreviousDensity, TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL MatrixMultiply(TempMat, SQRMat, NewDensity, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL MatrixMultiply(ISQMat, NewDensity, TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL MatrixMultiply(TempMat, ISQMat, NewDensity, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructMatrix(SQRMat)
    CALL DestructMatrix(ISQMat)
    CALL DestructMatrix(TempMat)
    CALL DestructSolverParameters(solver_parameters)

  END SUBROUTINE LowdinExtrapolate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE GeometryOptimizationModule
