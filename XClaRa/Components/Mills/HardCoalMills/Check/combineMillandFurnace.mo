within Exergy.XClaRa.Components.Mills.HardCoalMills.Check;
model combineMillandFurnace
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
  extends ClaRa.Basics.Icons.PackageIcons.ExecutableExampleb50;
  Modelica.Blocks.Sources.Ramp PTarget1(
    startTime=10000,
    duration=1000,
    offset=1,
    height=-0.2)
    annotation (Placement(transformation(extent={{-94,-10},{-74,10}})));
  BoundaryConditions.BoundaryGas_Txim_flow fluelGasFlowSource1(
    m_flow_const=2.2*6.7362,
    variable_m_flow=true,
    T_const=393.15,
    variable_xi=false,
    xi_const={0,0,0.0005,0,0.8,0.1985,0,0.001,0})
                    annotation (Placement(transformation(extent={{-50,-53},{-30,-33}})));
  BoundaryConditions.BoundaryFuel_Txim_flow coalFlowSource(m_flow_const=2, variable_m_flow=true,
    LHV_calculationType="Verbandsformel",
    xi_const={0.732,0.05,0.05,0.05,0.025,0.025})                                                 annotation (Placement(transformation(extent={{-48,-16},{-28,4}})));
  Modelica.Blocks.Math.Gain gain2(k=1200e6/30e6)
    "INIT.boiler.Q_nom/combustionChamber.LHV_fixed"                                                              annotation (Placement(transformation(extent={{-66,-5},{-56,5}})));
  Modelica.Blocks.Sources.RealExpression m_Primary2(y=-1.1*coalFlowSource.fuel_a.m_flow
        *15) "combustionChamber.m_flow_air_req*1.1"
    annotation (Placement(transformation(extent={{16,-14},{-16,14}},
        rotation=180,
        origin={-84,-37})));
  Exergy.XClaRa.Components.Mills.HardCoalMills.VerticalMill_L3 mills1(
    millKoeff=Fundamentals.STV1(),
    N_mills=1,
    initChoice=ClaRa.Basics.Choices.Init.noInit,
    T_0=363.15) annotation (Placement(transformation(extent={{30,-33},
            {50,-13}})));
  Modelica.Blocks.Sources.Ramp ramp1(
    offset=1.50,
    startTime=1000,
    duration=10,
    height=0)
    annotation (Placement(transformation(extent={{-20,20},{0,40}})));
  Exergy.XClaRa.Components.Adapters.FuelFlueGas_join coalGas_join_burner4
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={8,-23})));

  inner ClaRa.SimCenter simCenter(redeclare replaceable
      TILMedia.VLEFluidTypes.TILMedia_InterpolatedWater fluid1,
      redeclare ClaRa.Basics.Media.Fuel.Coal_v2 fuelModel1)
    annotation (Placement(transformation(extent={{-100,80},{-60,100}})));
  Furnace.SimpleCombustionChamber combustionChamber(xi_NOx=0) annotation (Placement(transformation(extent={{60,-33},{80,-13}})));
  BoundaryConditions.BoundaryGas_pTxi flueGasPressureSink(                                  p_const=100000, xi_const={0,0,0,0,0.79,0.21,0,0,0})
                                                                                        annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={70,18})));
  BoundaryConditions.BoundarySlag_pT coalSink annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={70,-58})));
equation
  connect(PTarget1.y, gain2.u)
                              annotation (Line(
      points={{-73,0},{-67,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain2.y, coalFlowSource.m_flow) annotation (Line(
      points={{-55.5,0},{-55.5,0.5},{-48,0.5},{-48,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ramp1.y, mills1.classifierSpeed) annotation (Line(
      points={{1,30},{40,30},{40,-12.2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(coalGas_join_burner4.fuelFlueGas_outlet, mills1.inlet) annotation (Line(
      points={{18,-23},{30,-23}},
      color={175,175,175},
      smooth=Smooth.None));
  connect(coalFlowSource.fuel_a,coalGas_join_burner4.fuel_inlet)   annotation (
      Line(
      points={{-28,-6},{-20,-6},{-20,-17},{-2,-17}},
      color={0,0,0},
      pattern=LinePattern.Solid,
      smooth=Smooth.None));
  connect(fluelGasFlowSource1.gas_a, coalGas_join_burner4.flueGas_inlet)
    annotation (Line(
      points={{-30,-43},{-20,-43},{-20,-29},{-2,-29}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(flueGasPressureSink.gas_a, combustionChamber.flueGas_outlet)
    annotation (Line(
      points={{70,8},{70,-13}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(coalSink.slag_inlet, combustionChamber.slag_outlet) annotation (Line(
      points={{70.2,-48},{70.2,-32.8},{70,-32.8}},
      color={127,127,0},
      smooth=Smooth.None));
  connect(mills1.outlet, combustionChamber.inlet) annotation (Line(
      points={{50,-23},{60,-23}},
      color={175,175,175},
      smooth=Smooth.None));
  connect(m_Primary2.y, fluelGasFlowSource1.m_flow) annotation (Line(
      points={{-66.4,-37},{-50,-37}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                      graphics={Text(
          extent={{-100,-84},{100,-100}},
          lineColor={0,128,0},
          lineThickness=0.5,
          fillColor={102,198,0},
          fillPattern=FillPattern.Solid,
          textString="IDEA:
shows how to combine a mill model to a simple furnace model",
          horizontalAlignment=TextAlignment.Left)}),
    experiment(StopTime=20000),
    __Dymola_experimentSetupOutput);
end combineMillandFurnace;
