within Exergy.XClaRa.Components.Mills.HardCoalMills;
model RollerBowlMill_L1
  "A simple pulveriser without classifier based on Dolezal"
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

  extends ClaRa.Basics.Icons.RollerBowlMill;
  extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L1");

  parameter Modelica.SIunits.Time Tau_m=100 "time constant of pulveriser";
  parameter Modelica.SIunits.MassFlowRate m_flow_dust_0= 10
    "Initial coal dust flow"                                                         annotation(Dialog(group="Initialisation"));
  Modelica.Blocks.Continuous.TransferFunction transferFunction(
    initType=Modelica.Blocks.Types.Init.InitialOutput,
    y_start=m_flow_dust_0,
    a={Tau_m*10,Tau_m,1})
    annotation (Placement(transformation(extent={{-40,0},{-20,20}})));
  Modelica.Blocks.Interfaces.RealInput rawCoal "Connector of Real input signal"
    annotation (Placement(transformation(extent={{-128,-10},{-88,30}})));
  Modelica.Blocks.Interfaces.RealOutput coalDust
    "Connector of Real output signal"
    annotation (Placement(transformation(extent={{100,0},{120,20}})));
equation
  connect(transferFunction.y, coalDust) annotation (Line(
      points={{-19,10},{110,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(rawCoal, transferFunction.u) annotation (Line(
      points={{-108,10},{-42,10}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(graphics));
end RollerBowlMill_L1;
