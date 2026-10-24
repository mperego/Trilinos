// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_ConfigDefs.hpp"

#if defined(HAVE_PIRO_TEMPUS) && defined(HAVE_PIRO_ROL)

#include "Piro_TempusSolver.hpp"
#include "Piro_ThyraProductME_Tempus_FinalObjective.hpp"
#include "Piro_Test_ThyraSupport.hpp"

#include "MassSpringDamperModel.hpp"

#include "ROL_ThyraVector.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_DefaultModelEvaluatorWithSolveFactory.hpp"
#include "Thyra_DefaultNominalBoundsOverrideModelEvaluator.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_YamlParameterListCoreHelpers.hpp"

#include "Tpetra_Core.hpp"

#include <cmath>

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using namespace Piro;
using namespace Piro::Test;

namespace {

typedef Thyra::ModelEvaluatorBase MEB;

// Linear solver of the models: GMRES without preconditioner, since the systems are 2x2.
// The adjoint integration uses the explicit adjoint model: the Stratimikos linear solvers
// do not support the transposed solves of the implicit adjoint.
RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactoryNew()
{
  const RCP<ParameterList> stratParams = rcp(new ParameterList("Stratimikos"));
  stratParams->set("Linear Solver Type", "Belos");
  stratParams->set("Preconditioner Type", "None");
  ParameterList& belosParams = stratParams->sublist("Linear Solver Types").sublist("Belos");
  belosParams.set("Solver Type", "Block GMRES");
  ParameterList& gmresParams = belosParams.sublist("Solver Types").sublist("Block GMRES");
  gmresParams.set("Convergence Tolerance", 1.0e-12);
  gmresParams.set("Maximum Iterations", 10);

  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
  linearSolverBuilder.setParameterList(stratParams);
  return createLinearSolveStrategy(linearSolverBuilder);
}

RCP<Thyra::ModelEvaluator<double> > withSolve(const RCP<Thyra::ModelEvaluator<double> >& model,
                                               const RCP<Thyra::LinearOpWithSolveFactoryBase<double> >& lowsFactory)
{
  return rcp(new Thyra::DefaultModelEvaluatorWithSolveFactory<double>(model, lowsFactory));
}

} // namespace

// The Tempus-driven objective (value and adjoint gradient) integrates from the consistent
// initial state given by the function of the Piro::TempusSolver: same results as with a
// model whose nominal initial state is the consistent one.
TEUCHOS_UNIT_TEST(Piro_TempusFinalObjective, ConsistentInitialState)
{
  const auto comm = Tpetra::getDefaultComm();
  const auto lowsFactory = lowsFactoryNew();
  const RCP<Thyra::ModelEvaluator<double> > model =
    withSolve(rcp(new MassSpringDamperModel(comm, false)), lowsFactory);
  const RCP<Thyra::ModelEvaluator<double> > adjointModel =
    withSolve(rcp(new MassSpringDamperModel(comm, true)), lowsFactory);

  // Backward Euler on [0, 0.2] with two time steps
  const RCP<ParameterList> appParams(new ParameterList("params"));
  Teuchos::updateParametersFromYamlFile("input_tempus_be_nox_solver.yaml", appParams.ptr());
  auto& tsc = appParams->sublist("Tempus").sublist("Demo Integrator").sublist("Time Step Control");
  const double timeStepSize = 0.1;
  tsc.set("Initial Time", 0.0);
  tsc.set("Final Time", 0.2);
  tsc.set("Minimum Time Step", timeStepSize/10);
  tsc.set("Initial Time Step", timeStepSize);
  tsc.set("Maximum Time Step", timeStepSize);

  // "Consistent" initial state, different from the nominal one (zero)
  const RCP<Thyra::VectorBase<double> > consistentState = model->getNominalValues().get_x()->clone_v();
  Thyra::set_ele(0, 0.3, consistentState.ptr());
  Thyra::set_ele(1, -0.2, consistentState.ptr());

  // Reference: the nominal initial state of the model is the consistent one
  const RCP<MEB::InArgs<double> > refNominalValues = rcp(new MEB::InArgs<double>(model->getNominalValues()));
  refNominalValues->set_x(consistentState);
  const RCP<Thyra::ModelEvaluator<double> > refModel =
    rcp(new Thyra::DefaultNominalBoundsOverrideModelEvaluator<double>(model, refNominalValues));

  // The objectives take the adjoint model directly: the solvers are only used for their
  // time interval, integrator and consistent initial state function
  const RCP<Thyra::ModelEvaluator<double> > noAdjointModel;
  const RCP<ObserverBase<double> > noObserver;
  const RCP<TempusSolver<double> > solver = tempusSolver(appParams, model, noAdjointModel, noObserver);
  const RCP<TempusSolver<double> > refSolver = tempusSolver(appParams, refModel, noAdjointModel, noObserver);

  // Parameters different from the nominal ones
  const RCP<Thyra::VectorBase<double> > p = model->getNominalValues().get_p(0)->clone_v();
  Thyra::set_ele(0, 1.2, p.ptr());
  Thyra::set_ele(1, 0.8, p.ptr());

  const double tol = 1.0e-8;
  int numCalls = 0;
  solver->setConsistentInitialStateFunction(
    [&](const RCP<Thyra::ModelEvaluator<double> >& icModel,
        const MEB::InArgs<double>& icInArgs,
        Thyra::VectorBase<double>& x0) {
      ++numCalls;
      TEST_ASSERT(icModel.get() == model.get());
      TEST_EQUALITY_CONST(icInArgs.get_t(), 0.0);
      TEST_COMPARE_FLOATING_ARRAYS(arrayFromVector(*icInArgs.get_p(0)), arrayFromVector(*p), tol);
      Thyra::copy(*consistentState, Teuchos::ptrFromRef(x0));
    });

  // The objective copies the function of the solver: create it after setting the function
  const int g_index = 0;
  const int Nt = 2;
  const RCP<Tempus::Integrator<double> > noIntegrator;
  ThyraProductME_TempusFinalObjective<double> objective(model, solver, noIntegrator, noIntegrator,
    adjointModel, g_index, *appParams, Nt, Teuchos::VERB_NONE);
  ThyraProductME_TempusFinalObjective<double> refObjective(refModel, refSolver, noIntegrator, noIntegrator,
    adjointModel, g_index, *appParams, Nt, Teuchos::VERB_NONE);
  ThyraProductME_TempusFinalObjective<double> nominalObjective(model, refSolver, noIntegrator, noIntegrator,
    adjointModel, g_index, *appParams, Nt, Teuchos::VERB_NONE);

  // The time integral depends on the initial state (trapezoidal rule, and the solution)
  objective.setResponseType(Tempus::AdjointResponseType::TimeIntegral);
  refObjective.setResponseType(Tempus::AdjointResponseType::TimeIntegral);
  nominalObjective.setResponseType(Tempus::AdjointResponseType::TimeIntegral);

  const ROL::ThyraVector<double> rol_p(p);
  double rol_tol = 1.0e-8;

  const double value = objective.value(rol_p, rol_tol);
  const double refValue = refObjective.value(rol_p, rol_tol);
  const double nominalValue = nominalObjective.value(rol_p, rol_tol);
  TEST_EQUALITY_CONST(numCalls, 1);
  TEST_FLOATING_EQUALITY(value, refValue, 1.0e-10);
  // Otherwise the test would not check that the consistent initial state is used
  TEST_ASSERT(std::abs(value - nominalValue) > 1.0e-6 * std::abs(refValue));

  const RCP<Thyra::VectorBase<double> > grad = p->clone_v();
  const RCP<Thyra::VectorBase<double> > refGrad = p->clone_v();
  ROL::ThyraVector<double> rol_grad(grad), rol_refGrad(refGrad);
  objective.gradient(rol_grad, rol_p, rol_tol);
  refObjective.gradient(rol_refGrad, rol_p, rol_tol);
  TEST_EQUALITY_CONST(numCalls, 2);
  TEST_COMPARE_FLOATING_ARRAYS(arrayFromVector(*grad), arrayFromVector(*refGrad), tol);
}

#endif /* HAVE_PIRO_TEMPUS && HAVE_PIRO_ROL */
