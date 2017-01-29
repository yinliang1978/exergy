within Exergy.XClaRa.Components.Utilities.Blocks.Check;
model testRandomGenerators2
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

  Modelica.Blocks.Sources.Ramp ramp1(
    duration=1000,
    height=1e6,
    offset=1e6,
    startTime=3001)
    annotation (Placement(transformation(extent={{-78,-52},{-58,-32}})));
  Modelica.Blocks.Sources.Ramp ramp3(
    height=1e5,
    duration=100,
    startTime=6000,
    offset=1e4)
    annotation (Placement(transformation(extent={{-36,-70},{-16,-50}})));
  Exergy.XClaRa.Components.Utilities.Blocks.Noise noise(
    mean_const=1e6,
    stdDev_const=1e4,
    varMeanValue=true,
    varStandardDeviation=true,
    startTime=1000,
    Tau_sample=300) annotation (Placement(transformation(extent={{4,
            -46},{24,-26}})));
equation

  connect(ramp1.y, noise.mean) annotation (Line(
      points={{-57,-42},{-28,-42},{-28,-30},{4,-30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ramp3.y, noise.sigma) annotation (Line(
      points={{-15,-60},{-6,-60},{-6,-40.2},{4,-40.2}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(graphics),
    experiment(StopTime=10000),
    __Dymola_experimentSetupOutput);
end testRandomGenerators2;
