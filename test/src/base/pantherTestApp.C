//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "pantherTestApp.h"
#include "pantherApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
pantherTestApp::validParams()
{
  InputParameters params = pantherApp::validParams();
  return params;
}

pantherTestApp::pantherTestApp(InputParameters parameters) : MooseApp(parameters)
{
  pantherTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

pantherTestApp::~pantherTestApp() {}

void
pantherTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  pantherApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"pantherTestApp"});
    Registry::registerActionsTo(af, {"pantherTestApp"});
  }
}

void
pantherTestApp::registerApps()
{
  registerApp(pantherApp);
  registerApp(pantherTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
pantherTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  pantherTestApp::registerAll(f, af, s);
}
extern "C" void
pantherTestApp__registerApps()
{
  pantherTestApp::registerApps();
}
