#include "pantherApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
pantherApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy material output, i.e., output properties on INITIAL as well as TIMESTEP_END
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

pantherApp::pantherApp(InputParameters parameters) : MooseApp(parameters)
{
  pantherApp::registerAll(_factory, _action_factory, _syntax);
}

pantherApp::~pantherApp() {}

void
pantherApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"pantherApp"});
  Registry::registerActionsTo(af, {"pantherApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
pantherApp::registerApps()
{
  registerApp(pantherApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
pantherApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  pantherApp::registerAll(f, af, s);
}
extern "C" void
pantherApp__registerApps()
{
  pantherApp::registerApps();
}
