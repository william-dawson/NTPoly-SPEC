# Generic modules
import sys

# NTPoly
import NTPolySwig as nt

# MPI Module
from mpi4py import MPI
comm = MPI.COMM_WORLD


##########################################################################
# This is the function we will map on to the matrix.
class TestOperation(nt.RealOperation):
    # Notice how you overload the call operator.
    def __call__(self):
        # This object contains a triplet called data for you to modify.
        if (self.data.index_row >= self.data.index_column):
            self.data.point_value *= 2
            return True
        return False

##########################################################################
if __name__ == "__main__":
    rank = comm.Get_rank()

    # Process The Input Parameters
    for i in range(1, len(sys.argv), 2):
        argument = sys.argv[i]
        argument_value = sys.argv[i + 1]
        if argument == '--input_matrix':
            input_matrix = argument_value
        elif argument == '--output_matrix':
            output_matrix = argument_value
        elif argument == '--process_slices':
            process_slices = int(argument_value)

    # Setup the process grid.
    nt.ConstructGlobalProcessGrid(process_slices)

    # Read in the matrices from file.
    Input = nt.Matrix_ps(input_matrix)
    Output = nt.Matrix_ps(Input.GetActualDimension())

    # Map
    op = TestOperation()
    nt.MatrixMapper.Map(Input, Output, op)

    # Print the density matrix to file.
    Output.WriteToMatrixMarket(output_matrix)
    nt.DestructGlobalProcessGrid()
