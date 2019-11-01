!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> The driver is a playground for an openacc approach.
!! Run it like this ./bin/ACCConcept A.mtx B.mtx number_of_blocks
!! Mainly you'll want to look in the subroutine MultiplyLoop unless you want
!! to start modifying the blocking.
PROGRAM PremadeMatrixProgram
  USE DataTypesModule, ONLY : NTREAL
  USE MatrixMemoryPoolModule, ONLY : MatrixMemoryPool_lr, &
       & ConstructMatrixMemoryPool, DestructMatrixMemoryPool
  USE SMatrixModule, ONLY : Matrix_lsr, ConstructMatrixFromFile, &
       & DestructMatrix, ConstructEmptyMatrix, SplitMatrix, &
       & ComposeMatrix, TransposeMatrix
  USE SMatrixAlgebraModule, ONLY : MatrixMultiply, IncrementMatrix, MatrixNorm
  IMPLICIT NONE
  !! Variables for handling input parameters.
  TYPE(Matrix_lsr) :: mata, matb, matc_new, matc_old
  CHARACTER(len=80) :: filea, fileb, argv
  INTEGER :: number_of_blocks
  REAL(NTREAL) :: normval
  REAL :: stime, ftime

  !! Read in the matrices from file.
  CALL get_command_argument(1, filea)
  CALL ConstructMatrixFromFile(mata, filea)
  CALL get_command_argument(2, fileb)
  CALL ConstructMatrixFromFile(matb, fileb)
  CALL get_command_argument(3, argv)
  READ(argv,*) number_of_blocks

  !! The new multiplication approach.
  CALL BlockMultiply(mata, matb, matc_new, number_of_blocks)

  !! Check the Result
  !! Note that you shouldn't compare times with this call as this includes
  !! all of the memory allocation overhead which normally won't happen.
  CALL MatrixMultiply(mata, matb, matc_old, threshold_in=1e-6_NTREAL)
  CALL IncrementMatrix(matc_old, matc_new, alpha_in=-1.0_NTREAL)
  normval = MatrixNorm(matc_new)

  WRITE(*,*) "Error", normval
  IF (normval > 1e-10) THEN
     STOP 1
  END IF

  !! Cleanup
  CALL DestructMatrix(mata)
  CALL DestructMatrix(matb)
  CALL DestructMatrix(matc_new)
  CALL DestructMatrix(matc_old)
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Handles all of the setup for a realistic matrix multiplication setup.
  SUBROUTINE BlockMultiply(mata, matb, matc, num_blocks)
    !> The matrices to multiply.
    TYPE(Matrix_lsr), INTENT(IN) :: mata, matb
    !> The results
    TYPE(Matrix_lsr), INTENT(INOUT) :: matc
    !> The number of blocks to split into.
    INTEGER, INTENT(IN) :: num_blocks
    !! A blocked representation of the two matrices.
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: mata_split
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: matb_split
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: mata_split_t
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: matb_split_t
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: matc_split
    !! Temporary variables
    TYPE(MatrixMemoryPool_lr), DIMENSION(:,:), ALLOCATABLE, TARGET :: mpool
    INTEGER :: II, JJ

    !! Split the matrices into blocks.
    CALL CPU_TIME(stime)
    ALLOCATE(mata_split(num_blocks,1))
    ALLOCATE(matb_split(1,num_blocks))
    ALLOCATE(mata_split_t(num_blocks,1))
    ALLOCATE(matb_split_t(1,num_blocks))
    ALLOCATE(matc_split(num_blocks, num_blocks))

    CALL SplitMatrix(mata, num_blocks, 1, mata_split)
    CALL SplitMatrix(matb, 1, num_blocks, matb_split)
    !! In the parallel algorithm the local blocks are all transposed already.
    DO II = 1, num_blocks
       CALL TransposeMatrix(mata_split(II,1), mata_split_t(II,1))
       CALL TransposeMatrix(matb_split(1,II), matb_split_t(1,II))
    END DO

    CALL ConstructEmptyMatrix(matc, mata%rows, matb%columns, zero_in=.TRUE.)
    CALL SplitMatrix(matc, num_blocks, num_blocks, matc_split)

    CALL CPU_TIME(ftime)
    WRITE(*,*) "Split overhead time:", ftime - stime

    !! Allocate the memory pool.
    CALL CPU_TIME(stime)
    ALLOCATE(mpool(num_blocks, num_blocks))
    CALL BuildMemoryPool(matc_split, mpool)
    CALL CPU_TIME(ftime)
    WRITE(*,*) "Memory pool allocation time:", ftime - stime

    !! Multiply.
    CALL MultiplyLoop(mata_split_t, matb_split_t, matc_split, num_blocks, &
         & mpool)

    !! Merge the matrix C
    CALL CPU_TIME(stime)
    CALL ComposeMatrix(matc_split, num_blocks, num_blocks, matc)
    CALL CPU_TIME(ftime)
    WRITE(*,*) "Merge overhead time:", ftime - stime

    !! Cleanup
    DO II = 1, num_blocks
       CALL DestructMatrix(mata_split(II,1))
       CALL DestructMatrix(matb_split(1,II))
       DO JJ = 1, num_blocks
          CALL DestructMatrix(matc_split(II,JJ))
          CALL DestructMatrixMemoryPool(mpool(II,JJ))
       END DO
    END DO
    DEALLOCATE(mata_split)
    DEALLOCATE(matb_split)
    DEALLOCATE(matc_split)
    DEALLOCATE(mpool)
  END SUBROUTINE BlockMultiply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE BuildMemoryPool(matc, pool)
    !> The matrices being multiplied.
    TYPE(Matrix_lsr), DIMENSION(:,:), INTENT(IN) :: matc
    !> The pool to setup.
    TYPE(MatrixMemoryPool_lr), DIMENSION(:,:), INTENT(INOUT), TARGET :: pool
    INTEGER :: II, JJ

    DO II = 1, SIZE(matc,DIM=1)
       DO JJ = 1, SIZE(matc,DIM=2)
          CALL ConstructMatrixMemoryPool(pool(II,JJ), matc(II,JJ)%rows, &
               & matc(II,JJ)%columns, 1.0_NTREAL)
       END DO
    END DO
  END SUBROUTINE BuildMemoryPool
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This is the part of the code to modify.
  SUBROUTINE MultiplyLoop(mata, matb, matc, num_blocks, pool)
    !> The matrices to multiply (matc = mata*matb)
    TYPE(Matrix_lsr), DIMENSION(:,:), INTENT(IN) :: mata, matb
    TYPE(Matrix_lsr), DIMENSION(:,:), INTENT(INOUT) :: matc
    !> Number of blocks
    INTEGER, INTENT(IN) :: num_blocks
    !> The memory pool to avoid allocations.
    TYPE(MatrixMemoryPool_lr), DIMENSION(:,:), INTENT(INOUT), TARGET :: pool
    INTEGER :: II, JJ

    CALL CPU_TIME(stime)
    !$OMP PARALLEL PRIVATE(JJ)
    !$OMP DO
    DO II = 1, num_blocks
       DO JJ = 1, num_blocks
          CALL MatrixMultiply(mata(II,1), matb(1,JJ), matc(II,JJ), &
               & IsATransposed_in=.TRUE., IsBTransposed_in=.TRUE., &
               & threshold_in=1e-6_NTREAL, blocked_memory_pool_in=pool(II,JJ))
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    CALL CPU_TIME(ftime)
    WRITE(*,*) "Loop time:", ftime - stime
  END SUBROUTINE MultiplyLoop
END PROGRAM PremadeMatrixProgram
