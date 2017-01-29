within Exergy.XClaRa.Components.VolumesValvesFittings.Fittings.Check;
model Test_Junction
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
  import ClaRa;
  ClaRa.Components.VolumesValvesFittings.Fittings.FlueGasJunction_L2
    flueGasJunction(
    volume=0.1)                            annotation (Placement(transformation(
        extent={{-11,-8},{11,8}},
        rotation=180,
        origin={-15,14})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow source2(
    variable_m_flow=false,
    variable_T=false,
    m_flow_const=0.1,
    xi_const={0.7,0.1,0.1,0.05,0.01,0.02,0.01,0.01,0},
    T_const=298.15)   annotation (Placement(transformation(extent={{-74,4},{-54,24}})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow source1(
    variable_m_flow=false,
    variable_T=false,
    m_flow_const=0.1,
    xi_const={0.7,0.1,0.1,0.05,0.01,0.02,0.01,0.01,0},
    T_const=298.15)   annotation (Placement(transformation(extent={{-74,26},{-54,46}})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow sink(
    variable_m_flow=false,
    variable_T=false,
    m_flow_const=-0.2,
    xi_const={0.7,0.1,0.1,0.05,0.01,0.02,0.01,0.01,0})
                       annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={56,14})));
  inner ClaRa.SimCenter simCenter annotation (Placement(transformation(extent={{78,80},{98,100}})));
  ClaRa.Components.VolumesValvesFittings.Fittings.FlueGasJunction_L2
    flueGasJunction1(
    volume=0.1)                            annotation (Placement(transformation(
        extent={{-11,-8},{11,8}},
        rotation=0,
        origin={-15,-8})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow source01(
    variable_m_flow=false,
    variable_T=false,
    m_flow_const=0.2,
    xi_const={0.7,0.1,0.1,0.05,0.01,0.02,0.01,0.01,0})
                      annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-64,-8})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow sink01(
    variable_m_flow=false,
    variable_T=false,
    m_flow_const=-0.1,
    xi_const={0.7,0.1,0.1,0.05,0.01,0.02,0.01,0.01,0})
                       annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={56,-8})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow sink02(
    variable_m_flow=false,
    variable_T=false,
    m_flow_const=-0.1,
    xi_const={0.7,0.1,0.1,0.05,0.01,0.02,0.01,0.01,0})
                       annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={56,-30})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow source3(
    variable_m_flow=false,
    variable_T=false,
    m_flow_const=10,
    xi_const={0.7,0.1,0.1,0.05,0.01,0.02,0.01,0.01,0})
                     annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-62,-38})));
  ClaRa.Components.BoundaryConditions.BoundaryGas_pTxi gasSink_pT1(p_const=100000, xi_const={0.7,0.1,0.1,0.05,0.01,0.02,0.01,0.01,0})
                                                                                   annotation (Placement(transformation(extent={{68,-66},{48,-46}})));
  ClaRa.Components.VolumesValvesFittings.Fittings.FlueGasJunction_L2
    flueGasJunction2(
    volume=1)                              annotation (Placement(transformation(
        extent={{-11,-8},{11,8}},
        rotation=0,
        origin={29,-56})));
  ClaRa.Components.VolumesValvesFittings.Fittings.FlueGasJunction_L2
    flueGasJunction3(
    volume=1)                              annotation (Placement(transformation(
        extent={{11,8},{-11,-8}},
        rotation=0,
        origin={-33,-56})));
  ClaRa.Components.VolumesValvesFittings.Valves.ValveGas_L1
    flueGasValve_L1_2(
    redeclare model PressureLoss =
        ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.QuadraticZeta
        (A_cross=
           0.05, zeta=0.0002))
    annotation (Placement(transformation(extent={{-12,-90},{8,-80}})));
  ClaRa.Basics.ControlVolumes.GasVolumes.VolumeGas_L2 flueGasCell1(
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.GenericGeometry,
    m_flow_nom=0.5,
    redeclare model PressureLoss =
        ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.NoFriction_L2)
                                                                                                        annotation (Placement(transformation(extent={{-42,-96},{-22,-76}})));

  ClaRa.Components.VolumesValvesFittings.Valves.ValveGas_L1
    flueGasValve_L1_1(
    redeclare model PressureLoss =
        ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.QuadraticZeta
        (A_cross=
           0.01, zeta=0.0005))
    annotation (Placement(transformation(extent={{-12,-60},{8,-50}})));
equation
  connect(source2.gas_a, flueGasJunction.portB) annotation (Line(
      points={{-54,14},{-26,14}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(source1.gas_a, flueGasJunction.portC) annotation (Line(
      points={{-54,36},{-15,36},{-15,22}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(flueGasJunction.portA, sink.gas_a) annotation (Line(
      points={{-4,14},{46,14}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(sink01.gas_a, flueGasJunction1.portB) annotation (Line(
      points={{46,-8},{-4,-8}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(sink02.gas_a, flueGasJunction1.portC) annotation (Line(
      points={{46,-30},{-15,-30},{-15,-16}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(source01.gas_a, flueGasJunction1.portA) annotation (Line(
      points={{-54,-8},{-26,-8}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(source3.gas_a, flueGasJunction3.portC) annotation (Line(
      points={{-52,-38},{-33,-38},{-33,-48}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(flueGasJunction2.portB, gasSink_pT1.gas_a) annotation (Line(
      points={{40,-56},{48,-56}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(flueGasJunction3.portA, flueGasValve_L1_1.inlet) annotation (Line(
      points={{-22,-56},{-18,-56},{-18,-55.8333},{-12,-55.8333}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(flueGasValve_L1_1.outlet, flueGasJunction2.portA) annotation (Line(
      points={{8,-55.8333},{14,-55.8333},{14,-56},{18,-56}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(flueGasValve_L1_2.inlet, flueGasCell1.outlet) annotation (Line(
      points={{-12,-85.8333},{-12,-86},{-22,-86}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(flueGasCell1.inlet, flueGasJunction3.portB) annotation (Line(
      points={{-42,-86},{-52,-86},{-52,-56},{-44,-56}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(flueGasValve_L1_2.outlet, flueGasJunction2.portC) annotation (Line(
      points={{8,-85.8333},{29,-85.8333},{29,-64}},
      color={84,58,36},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                      graphics={  Text(
          extent={{-98,98},{100,58}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=10,
          textString="______________________________________________________________________________________________
PURPOSE:
>> Tester for gas junctions

______________________________________________________________________________________________
"),                    Text(
          extent={{-98,102},{16,82}},
          lineColor={0,128,0},
          fontSize=31,
          textString="TESTED -- 2015-01-27 //LN")}));
end Test_Junction;
