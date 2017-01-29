within Exergy.XClaRa.Components.Utilities.Blocks.Check;
model TestParameterizableTable1D
//___________________________________________________________________________//
// Component of the ClaRa library, version: 1.1.2                        //
//                                                                           //
// Licensed by the DYNCAP/DYNSTART research team under Modelica License 2.   //
// Copyright © 2013-2016, DYNCAP/DYNSTART research team.                     //
//___________________________________________________________________________//
// DYNCAP and DYNSTART are research projects supported by the German Federal //
// Ministry of Economic Affairs and Energy (FKZ 03ET2009/FKZ 03ET7060).      //
// The research team consists of the following project partners:             //
// Institute of Energy Systems (Hamburg University of Technology),           //
// Institute of Thermo-Fluid Dynamics (Hamburg University of Technology),    //
// TLK-Thermo GmbH (Braunschweig, Germany),                                  //
// XRG Simulation GmbH (Hamburg, Germany).                                   //
//___________________________________________________________________________//
//
// This model demonstrates the following usecase of a CombiTable1D:
// The parameter "table" of the CombiTable depends on the value of parameter "a" which is calculated by solving an initial equation system.
//
// This usecase is not supported by OpenModelica or SimulationX.
//
// As a workaround use ClaRa.Components.Utilities.Blocks.ParameterizableTable1D
// instead of Modelica.Blocks.Tables.CombiTable1D.
//
// For a realworld example see: ClaRa.Examples.ClaRaClosedLoop, which uses ClaRa.SubSystems.Boiler.SteamGenerator_L3
  parameter Real a(fixed=false);
  parameter Real b(fixed=false);
  // Modelica.Blocks.Tables.CombiTable1D table(table=[0,0;1,a]); // This line will fail with OpenModelica or SimulationX. Therefor use ParameterizableTable1D instead of CombiTable1D.
  ParameterizableTable1D table(table=[0,0;1,a]);
  Real a_times_table_input=a * table.u[1];
  Real table_result=table.y[1];
initial equation
    0 = sin(a) + b;
    0 = 3*a + 7*b + 1;
equation
    table.u[1] = 0.5;
    assert( abs(a_times_table_input - table_result) < 0.01, "Fehler table.y[1] != table.u[1] * a");
end TestParameterizableTable1D;
