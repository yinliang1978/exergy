within Exergy.XClaRa.Components.Utilities.Blocks.Check;
model TestMinMax
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

  import SI = Modelica.SIunits;

  Modelica.Blocks.Sources.Ramp ramp(
    duration=10,
    startTime=100,
    height=+420)
    annotation (Placement(transformation(extent={{-68,26},{-48,46}})));
  Modelica.Blocks.Sources.Sine sine(
    amplitude=25,
    freqHz=0.1,
    offset=-300)
    annotation (Placement(transformation(extent={{-62,-24},{-42,-4}})));
  Modelica.Blocks.Math.Add T
    annotation (Placement(transformation(extent={{-20,6},{0,26}})));
  Exergy.XClaRa.Components.Utilities.Blocks.TimeExtrema timeExtrema
    annotation (Placement(transformation(extent={{20,6},{40,26}})));
initial equation

equation
  connect(ramp.y, T.u1) annotation (Line(
      points={{-47,36},{-34,36},{-34,22},{-22,22}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sine.y, T.u2) annotation (Line(
      points={{-41,-14},{-32,-14},{-32,10},{-22,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(T.y, timeExtrema.u) annotation (Line(
      points={{1,16},{18,16}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (experiment(StopTime=300), experimentSetupOutput);
end TestMinMax;
