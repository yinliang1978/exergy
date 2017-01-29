within Exergy.XClaRa.Components.BoundaryConditions.Check;
model TestIAPWSboundaries2
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
  extends ClaRa.Basics.Icons.PackageIcons.ExecutableExampleb60;
  inner ClaRa.SimCenter simCenter(redeclare replaceable
      TILMedia.VLEFluidTypes.TILMedia_InterpolatedWater fluid1)
    annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
  Modelica.Blocks.Sources.CombiTimeTable ramp(table=[0,5; 1,5; 2,-5; 3,-5; 3.1,
        10; 4,10], extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
    annotation (Placement(transformation(extent={{-92,14},{-72,34}})));
  Modelica.Blocks.Sources.CombiTimeTable ramp1(extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic, table=[0,
        2e5; 0.5,2e5; 2,5e5; 3,5e5; 3.1,3e5; 4,3e5])
    annotation (Placement(transformation(extent={{18,-66},{38,-46}})));
  BoundaryVLE_hxim_flow massFlowSource_T(variable_m_flow=true) annotation (Placement(transformation(extent={{-52,-20},{-32,0}})));
  BoundaryVLE_phxi pressureSink_pT(variable_p=true, Delta_p=1000) annotation (Placement(transformation(extent={{60,-20},{38,0}})));
equation
  connect(ramp1.y[1], pressureSink_pT.p) annotation (Line(
      points={{39,-56},{78,-56},{78,-4},{60,-4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowSource_T.m_flow, ramp.y[1]) annotation (Line(
      points={{-54,-4},{-64,-4},{-64,24},{-71,24}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(pressureSink_pT.steam_a, massFlowSource_T.steam_a) annotation (Line(
      points={{38,-10},{-32,-10}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (Diagram(graphics));
end TestIAPWSboundaries2;
