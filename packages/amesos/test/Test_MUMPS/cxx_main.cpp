#include "Amesos_ConfigDefs.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Util.h"
#include "Amesos_Mumps.h"
#include "Amesos_TestRowMatrix.h"
#include "Teuchos_ParameterList.hpp"
//  using namespace Trilinos_Util;   commented out to resolve bug #1886

//=============================================================================
bool CheckError(const Epetra_RowMatrix& A,
		const Epetra_MultiVector& x,
		const Epetra_MultiVector& b,
		const Epetra_MultiVector& x_exact)
{
  std::vector<double> Norm;
  int NumVectors = x.NumVectors();
  Norm.resize(NumVectors);
  Epetra_MultiVector Ax(x);
  A.Multiply(false,x,Ax);
  Ax.Update(1.0,b,-1.0);
  Ax.Norm2(&Norm[0]);
  bool TestPassed = false;
  double TotalNorm = 0.0;
  for (int i = 0 ; i < NumVectors ; ++i) {
    TotalNorm += Norm[i];
  }
  if (A.Comm().MyPID() == 0)
    std::cout << "||Ax - b||  = " << TotalNorm << std::endl;
  if (TotalNorm < 1e-5 )
    TestPassed = true;
  else
    TestPassed = false;

  Ax.Update (1.0,x,-1.0,x_exact,0.0);
  Ax.Norm2(&Norm[0]);
  for (int i = 0 ; i < NumVectors ; ++i) {
    TotalNorm += Norm[i];
  }
  if (A.Comm().MyPID() == 0)
    std::cout << "||Ax - b||  = " << TotalNorm << std::endl;
#ifdef HAVE_AMESOS_SMUMPS
  if (TotalNorm < 1e-2 )
#else
  if (TotalNorm < 1e-5 )
#endif
    TestPassed = true;
  else
    TestPassed = false;

  return(TestPassed);
}

//=============================================================================
int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int NumGlobalElements = 100;
  int NumVectors = 7;
  
  // =================== //
  // create a random map //
  // =================== //

  int* part = new int[NumGlobalElements];

  if (Comm.MyPID() == 0) {
    Epetra_Util Util;

    for( int i=0 ; i<NumGlobalElements ; ++i ) {
      unsigned int r = Util.RandomInt();
      part[i] = r%(Comm.NumProc());
    }
  }

  Comm.Broadcast(part,NumGlobalElements,0);

  // count the elements assigned to this proc
  int NumMyElements = 0;
  for (int i = 0 ; i < NumGlobalElements ; ++i) {
    if (part[i] == Comm.MyPID())
      NumMyElements++;
  }

  // get the loc2global list
  int* MyGlobalElements = new int[NumMyElements];
  int count = 0;
  for (int i = 0 ; i < NumGlobalElements ; ++i) {
    if (part[i] == Comm.MyPID() )
      MyGlobalElements[count++] = i;
  }

  Epetra_Map Map(NumGlobalElements,NumMyElements,MyGlobalElements,
                 0,Comm);

  delete [] part;
  
  // ===================== //
  // Create a dense matrix //
  // ===================== //
 
  Epetra_CrsMatrix Matrix(Copy,Map,NumGlobalElements);

  int* Indices = new int[NumGlobalElements];
  double* Values = new double[NumGlobalElements];

  for (int i = 0 ; i < NumGlobalElements ; ++i) 
    Indices[i] = i;

  for (int i = 0 ; i < NumMyElements ; ++i) {
    int iGlobal = MyGlobalElements[i];
    for (int jGlobal = 0 ; jGlobal < NumGlobalElements ; ++jGlobal) {
      if (iGlobal >= jGlobal) 
	Values[jGlobal] = 1.0 * (jGlobal + 1);
      else
	Values[jGlobal] = 1.0 * (iGlobal + 1);
    }
    Matrix.InsertGlobalValues(MyGlobalElements[i],
                              NumGlobalElements, Values, Indices);

  }

  Matrix.FillComplete();

  delete [] Indices;
  delete [] Values;
 
  // ======================== //
  // other data for this test //
  // ======================== //

  Amesos_TestRowMatrix A(&Matrix);
  Epetra_MultiVector x(Map,NumVectors);
  Epetra_MultiVector x_exact(Map,NumVectors);
  Epetra_MultiVector b(Map,NumVectors);
  x_exact.Random();
  A.Multiply(false,x_exact,b);

  // =========== //
  // AMESOS PART //
  // =========== //

  Epetra_LinearProblem Problem;
  Amesos_Mumps* Solver = new Amesos_Mumps(Problem);

  Problem.SetOperator(&A);
  Problem.SetLHS(&x);
  Problem.SetRHS(&b);

  Teuchos::ParameterList List;
  List.set("MaxProcs",2);
  AMESOS_CHK_ERR(Solver->SetParameters(List));

  AMESOS_CHK_ERR(Solver->SymbolicFactorization());
  AMESOS_CHK_ERR(Solver->NumericFactorization());
  AMESOS_CHK_ERR(Solver->Solve());

  bool TestPassed = true;

  TestPassed = TestPassed &&
    CheckError(A,x,b,x_exact);

  delete Solver;

  AMESOS_CHK_ERR( ! TestPassed ) ; 

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
