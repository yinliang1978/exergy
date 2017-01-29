within Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Check;
model Test_GasValves
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
  extends ClaRa.Basics.Icons.PackageIcons.ExecutableRegressiong100;

 model Regression
  extends ClaRa.Basics.Icons.RegressionSummary;

  Modelica.Blocks.Interfaces.RealInput Delta_p_Valve_1
      "Pressure loss of valve 1";
  Modelica.Blocks.Interfaces.RealInput Delta_p_Valve_2
      "Pressure loss of valve 2";
  Modelica.Blocks.Interfaces.RealInput Delta_p_Valve_3
      "Pressure loss of valve 3";
  Modelica.Blocks.Interfaces.RealInput Delta_p_Valve_4
      "Pressure loss of valve 4";
  Modelica.Blocks.Interfaces.RealInput Delta_p_Valve_5
      "Pressure loss of valve 5";
  Modelica.Blocks.Interfaces.RealInput Delta_p_Valve_6
      "Pressure loss of valve 6";
  Modelica.Blocks.Interfaces.RealInput Delta_p_coalDustValve_1
      "Pressure loss of coalDustValve 1";
  Modelica.Blocks.Interfaces.RealInput Delta_p_coalDustValve_2
      "Pressure loss of coalDustValve 2";
  Modelica.Blocks.Interfaces.RealInput Delta_p_coalDustValve_3
      "Pressure loss of coalDustValve 3";
  Modelica.Blocks.Interfaces.RealInput Delta_p_coalDustValve_4
      "Pressure loss of coalDustValve 4";
  Modelica.Blocks.Interfaces.RealInput Delta_p_coalDustValve_5
      "Pressure loss of coalDustValve 5";
  Modelica.Blocks.Interfaces.RealInput Delta_p_coalDustValve_6
      "Pressure loss of coalDustValve 6";

  Real y_Delta_p_Valve_1 = integrator_Delta_p_Valve_1.y;
  Real y_Delta_p_Valve_2 = integrator_Delta_p_Valve_2.y;
  Real y_Delta_p_Valve_3 = integrator_Delta_p_Valve_3.y;
  Real y_Delta_p_Valve_4 = integrator_Delta_p_Valve_4.y;
  Real y_Delta_p_Valve_5 = integrator_Delta_p_Valve_5.y;
  Real y_Delta_p_Valve_6 = integrator_Delta_p_Valve_6.y;
  Real y_Delta_p_coalDustValve_1 = integrator_Delta_p_coalDustValve_1.y;
  Real y_Delta_p_coalDustValve_2 = integrator_Delta_p_coalDustValve_2.y;
  Real y_Delta_p_coalDustValve_3 = integrator_Delta_p_coalDustValve_3.y;
  Real y_Delta_p_coalDustValve_4 = integrator_Delta_p_coalDustValve_4.y;
  Real y_Delta_p_coalDustValve_5 = integrator_Delta_p_coalDustValve_5.y;
  Real y_Delta_p_coalDustValve_6 = integrator_Delta_p_coalDustValve_6.y;

  protected
  Utilities.Blocks.Integrator integrator_Delta_p_Valve_1(u=
          Delta_p_Valve_1) annotation (Placement(transformation(
            extent={{-80,66},{-60,86}})));
  Utilities.Blocks.Integrator integrator_Delta_p_Valve_2(u=
          Delta_p_Valve_2) annotation (Placement(transformation(
            extent={{-80,-28},{-60,-8}})));
  Utilities.Blocks.Integrator integrator_Delta_p_Valve_3(u=
          Delta_p_Valve_3) annotation (Placement(transformation(
            extent={{-80,4},{-60,24}})));
  Utilities.Blocks.Integrator integrator_Delta_p_Valve_4(u=
          Delta_p_Valve_4) annotation (Placement(transformation(
            extent={{-80,-60},{-60,-40}})));
  Utilities.Blocks.Integrator integrator_Delta_p_Valve_5(u=
          Delta_p_Valve_5) annotation (Placement(transformation(
            extent={{-80,34},{-60,54}})));
  Utilities.Blocks.Integrator integrator_Delta_p_Valve_6(u=
          Delta_p_Valve_6) annotation (Placement(transformation(
            extent={{-78,-94},{-58,-74}})));
  Utilities.Blocks.Integrator integrator_Delta_p_coalDustValve_1(u=
          Delta_p_coalDustValve_1) annotation (Placement(
          transformation(extent={{14,-96},{34,-76}})));
  Utilities.Blocks.Integrator integrator_Delta_p_coalDustValve_2(u=
          Delta_p_coalDustValve_2) annotation (Placement(
          transformation(extent={{12,-62},{32,-42}})));
  Utilities.Blocks.Integrator integrator_Delta_p_coalDustValve_3(u=
          Delta_p_coalDustValve_3) annotation (Placement(
          transformation(extent={{12,-28},{32,-8}})));
  Utilities.Blocks.Integrator integrator_Delta_p_coalDustValve_4(u=
          Delta_p_coalDustValve_4) annotation (Placement(
          transformation(extent={{14,0},{34,20}})));
  Utilities.Blocks.Integrator integrator_Delta_p_coalDustValve_5(u=
          Delta_p_coalDustValve_5) annotation (Placement(
          transformation(extent={{14,34},{34,54}})));
  Utilities.Blocks.Integrator integrator_Delta_p_coalDustValve_6(u=
          Delta_p_coalDustValve_6) annotation (Placement(
          transformation(extent={{14,66},{34,86}})));

 end Regression;

  ValveGas_L1     valve1(
    openingInputIsActive=true,
    redeclare model PressureLoss =
        ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint
        (Delta_p_nom=2e5, m_flow_nom=10))
    annotation (Placement(transformation(extent={{-84,126},{-64,138}})));

  inner ClaRa.SimCenter simCenter(redeclare TILMedia.GasTypes.FlueGasTILMedia
                                        flueGasModel) annotation (
      Placement(transformation(extent={{162,232},{182,252}})));

  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T(
    m_flow_const=10,
    T_const=293.15,
    xi_const={0,0,0.0,0,0.77,0.23,0,0,0},
    gas_a(p(start=1000000)))              annotation (Placement(transformation(extent={{-114,122},{-94,142}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT(p_const=800000, xi_const={0,0,0.0,0,0.77,0.23,0,0,0})
                                                                 annotation (Placement(transformation(extent={{78,122},{58,142}})));
  Modelica.Blocks.Sources.Ramp ramp(
    duration=8,
    offset=1,
    startTime=1,
    height=-0.9)
    annotation (Placement(transformation(extent={{-114,154},{-94,174}})));
  ValveGas_L1     valve2(
    openingInputIsActive=true,
    useHomotopy=false,
    redeclare model PressureLoss = Fundamentals.QuadraticKV (Kvs=100))
    annotation (Placement(transformation(extent={{-64,98},{-44,110}})));
  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T1(
    m_flow_const=10,
    T_const=293.15,
    xi_const={0,0,0.0,0,0.77,0.23,0,0,0},
    gas_a(p(start=1000000)))              annotation (Placement(transformation(extent={{-114,94},{-94,114}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT1(p_const=800000, xi_const={0,0,0.0,0,0.77,0.23,0,0,0})
                                                                  annotation (Placement(transformation(extent={{78,94},{58,114}})));
  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T2(
    m_flow_const=10,
    T_const=293.15,
    xi_const={0,0,0.0,0,0.77,0.23,0,0,0},
    gas_a(p(start=1000000)))              annotation (Placement(transformation(extent={{-114,66},{-94,86}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT2(p_const=800000, xi_const={0,0,0.0,0,0.77,0.23,0,0,0})
                                                                  annotation (Placement(transformation(extent={{78,66},{58,86}})));
  BoundaryConditions.BoundaryFuel_Txim_flow coalFlowSource(m_flow_const=5, T_const=293.15) annotation (Placement(transformation(extent={{-146,-160},{-126,-140}})));
  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T3(
    m_flow_const=10,
    T_const=293.15,
    xi_const={0,0,0.0,0,0.77,0.23,0,0,0},
    gas_a(p(start=1000000)))              annotation (Placement(transformation(extent={{-146,-182},{-126,-162}})));
  Adapters.FuelFlueGas_join coalGas_join annotation (Placement(transformation(extent={{-114,-172},{-94,-152}})));
  Adapters.FuelFlueGas_split coalGas_split annotation (Placement(transformation(extent={{62,-172},{82,-152}})));
  BoundaryConditions.BoundaryFuel_pTxi coalSink(p_const=800000) annotation (Placement(transformation(extent={{116,-160},{96,-140}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT3(p_const=800000) annotation (Placement(transformation(extent={{116,-184},{96,-164}})));
  BoundaryConditions.BoundaryFuel_Txim_flow coalFlowSource1(m_flow_const=5, T_const=293.15) annotation (Placement(transformation(extent={{-146,-108},{-126,-88}})));
  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T4(
    m_flow_const=10,
    T_const=293.15,
    xi_const={0,0,0.0,0,0.77,0.23,0,0,0},
    gas_a(p(start=1000000)))              annotation (Placement(transformation(extent={{-146,-130},{-126,-110}})));
  Adapters.FuelFlueGas_join coalGas_join1 annotation (Placement(transformation(extent={{-114,-120},{-94,-100}})));
  Adapters.FuelFlueGas_split coalGas_split1 annotation (Placement(transformation(extent={{62,-120},{82,-100}})));
  BoundaryConditions.BoundaryFuel_pTxi coalSink1(p_const=800000) annotation (Placement(transformation(extent={{118,-108},{98,-88}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT4(p_const=800000) annotation (Placement(transformation(extent={{118,-132},{98,-112}})));
  ValveFuelFlueGas_L1 coalDustValve2(openingInputIsActive=true, redeclare model
      PressureLoss =
        Fundamentals.QuadraticKV (                                                                                       Kvs=100))
                                                                                        annotation (Placement(transformation(extent={{-64,-116},{-44,-104}})));
  BoundaryConditions.BoundaryFuel_Txim_flow coalFlowSource2(m_flow_const=5, T_const=293.15) annotation (Placement(transformation(extent={{-146,-60},{-126,-40}})));
  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T5(
    m_flow_const=10,
    T_const=293.15,
    xi_const={0,0,0.0,0,0.77,0.23,0,0,0},
    gas_a(p(start=1000000)))              annotation (Placement(transformation(extent={{-146,-82},{-126,-62}})));
  Adapters.FuelFlueGas_join coalGas_join2 annotation (Placement(transformation(extent={{-114,-72},{-94,-52}})));
  Adapters.FuelFlueGas_split coalGas_split2 annotation (Placement(transformation(extent={{64,-72},{84,-52}})));
  BoundaryConditions.BoundaryFuel_pTxi coalSink2(p_const=800000) annotation (Placement(transformation(extent={{118,-60},{98,-40}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT5(p_const=800000) annotation (Placement(transformation(extent={{118,-84},{98,-64}})));
  ValveFuelFlueGas_L1 coalDustValve1(openingInputIsActive=true, redeclare model
      PressureLoss =
        ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint
        (                                                                                                    m_flow_nom=10, Delta_p_nom=2e5)) annotation (Placement(transformation(extent={{-84,-68},{-64,-56}})));

  ValveGas_L1     valve3(
    openingInputIsActive=true,
    useHomotopy=false,
    redeclare model PressureLoss =
        ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.QuadraticNominalPoint
        (
        Delta_p_nom=2e5,
        rho_in_nom=1,
        m_flow_nom=10))
    annotation (Placement(transformation(extent={{-44,70},{-24,82}})));
  ValveGas_L1     valve4(
    openingInputIsActive=true,
    useHomotopy=false,
    redeclare model PressureLoss =
        ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.QuadraticZeta
        (A_cross=
           0.2, zeta=0.5))
    annotation (Placement(transformation(extent={{-24,42},{-4,54}})));
  ValveGas_L1     valve5(
    openingInputIsActive=true,
    useHomotopy=true,
    redeclare model PressureLoss = Fundamentals.Quadratic_EN60534 (Kvs=100))
    annotation (Placement(transformation(extent={{-4,14},{16,26}})));

  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T6(
    m_flow_const=10,
    T_const=293.15,
    xi_const={0,0,0.0,0,0.77,0.23,0,0,0},
    gas_a(p(start=1000000)))              annotation (Placement(transformation(extent={{-114,38},{-94,58}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT6(p_const=800000, xi_const={0,0,0.0,0,0.77,0.23,0,0,0})
                                                                  annotation (Placement(transformation(extent={{78,38},{58,58}})));
  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T7(
    m_flow_const=10,
    T_const=293.15,
    xi_const={0,0,0.0,0,0.77,0.23,0,0,0},
    gas_a(p(start=1000000)))              annotation (Placement(transformation(extent={{-114,10},{-94,30}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT7(p_const=800000, xi_const={0,0,0.0,0,0.77,0.23,0,0,0})
                                                                  annotation (Placement(transformation(extent={{78,10},{58,30}})));
  ValveGas_L1     valve6(
    openingInputIsActive=true,
    redeclare model PressureLoss =
        ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.Quadratic_FlowFunction,
    useHomotopy=true)
    annotation (Placement(transformation(extent={{16,-14},{36,-2}})));

  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T8(
    m_flow_const=10,
    T_const=293.15,
    xi_const={0,0,0.0,0,0.77,0.23,0,0,0},
    gas_a(p(start=1000000)))              annotation (Placement(transformation(extent={{-114,-18},{-94,2}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT8(p_const=800000, xi_const={0,0,0.0,0,0.77,0.23,0,0,0})
                                                                  annotation (Placement(transformation(extent={{78,-18},{58,2}})));
  ValveFuelFlueGas_L1 coalDustValve3(
    redeclare model PressureLoss =
        ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.QuadraticNominalPoint
        (
        Delta_p_nom=2e5,
        rho_in_nom=1,
        m_flow_nom=10),
    openingInputIsActive=true,
    checkValve=false)
                     annotation (Placement(transformation(extent={{-44,-168},{-24,-156}})));
  BoundaryConditions.BoundaryFuel_Txim_flow coalFlowSource3(m_flow_const=5, T_const=293.15) annotation (Placement(transformation(extent={{-146,-264},{-126,-244}})));
  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T9(
    m_flow_const=10,
    T_const=293.15,
    xi_const={0,0,0.0,0,0.77,0.23,0,0,0},
    gas_a(p(start=1000000)))              annotation (Placement(transformation(extent={{-146,-286},{-126,-266}})));
  Adapters.FuelFlueGas_join coalGas_join3 annotation (Placement(transformation(extent={{-114,-276},{-94,-256}})));
  Adapters.FuelFlueGas_split coalGas_split3 annotation (Placement(transformation(extent={{60,-276},{80,-256}})));
  BoundaryConditions.BoundaryFuel_pTxi coalSink3(p_const=800000) annotation (Placement(transformation(extent={{114,-264},{94,-244}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT9(p_const=800000) annotation (Placement(transformation(extent={{114,-288},{94,-268}})));
  BoundaryConditions.BoundaryFuel_Txim_flow coalFlowSource4(m_flow_const=5, T_const=293.15) annotation (Placement(transformation(extent={{-146,-212},{-126,-192}})));
  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T10(
    m_flow_const=10,
    T_const=293.15,
    xi_const={0,0,0.0,0,0.77,0.23,0,0,0},
    gas_a(p(start=1000000)))              annotation (Placement(transformation(extent={{-146,-234},{-126,-214}})));
  Adapters.FuelFlueGas_join coalGas_join4 annotation (Placement(transformation(extent={{-114,-224},{-94,-204}})));
  Adapters.FuelFlueGas_split coalGas_split4 annotation (Placement(transformation(extent={{60,-224},{80,-204}})));
  BoundaryConditions.BoundaryFuel_pTxi coalSink4(p_const=800000) annotation (Placement(transformation(extent={{116,-212},{96,-192}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT10(p_const=800000) annotation (Placement(transformation(extent={{116,-236},{96,-216}})));
  ValveFuelFlueGas_L1 coalDustValve4(openingInputIsActive=true, redeclare model
      PressureLoss =
        ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.QuadraticZeta
        (                                                                                                    zeta=0.5, A_cross=0.2)) annotation (Placement(transformation(extent={{-24,-220},{-4,-208}})));
  ValveFuelFlueGas_L1 coalDustValve5(
    openingInputIsActive=true,
    useHomotopy=true,
    redeclare model PressureLoss = Fundamentals.Quadratic_EN60534 (Kvs=100))
                      annotation (Placement(transformation(extent={{-4,-272},{16,-260}})));

  BoundaryConditions.BoundaryFuel_Txim_flow coalFlowSource5(m_flow_const=5, T_const=293.15) annotation (Placement(transformation(extent={{-146,-314},{-126,-294}})));
  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T11(
    m_flow_const=10,
    T_const=293.15,
    xi_const={0,0,0.0,0,0.77,0.23,0,0,0},
    gas_a(p(start=1000000)))              annotation (Placement(transformation(extent={{-146,-336},{-126,-316}})));
  Adapters.FuelFlueGas_join coalGas_join5 annotation (Placement(transformation(extent={{-114,-326},{-94,-306}})));
  Adapters.FuelFlueGas_split coalGas_split5 annotation (Placement(transformation(extent={{60,-326},{80,-306}})));
  BoundaryConditions.BoundaryFuel_pTxi coalSink5(p_const=800000) annotation (Placement(transformation(extent={{114,-314},{94,-294}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT11(p_const=800000) annotation (Placement(transformation(extent={{114,-338},{94,-318}})));
  ValveFuelFlueGas_L1 coalDustValve6(
    openingInputIsActive=true,
    redeclare model PressureLoss =
        ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.Quadratic_FlowFunction,
    useHomotopy=true) annotation (Placement(transformation(extent={{16,-322},{36,-310}})));

  Regression regression(
    Delta_p_Valve_1 = valve1.summary.outline.Delta_p,
    Delta_p_Valve_2 = valve2.summary.outline.Delta_p,
    Delta_p_Valve_3 = valve3.summary.outline.Delta_p,
    Delta_p_Valve_4 = valve4.summary.outline.Delta_p,
    Delta_p_Valve_5 = valve5.summary.outline.Delta_p,
    Delta_p_Valve_6 = valve6.summary.outline.Delta_p,
    Delta_p_coalDustValve_1 = coalDustValve1.summary.outline.Delta_p,
    Delta_p_coalDustValve_2 = coalDustValve2.summary.outline.Delta_p,
    Delta_p_coalDustValve_3 = coalDustValve3.summary.outline.Delta_p,
    Delta_p_coalDustValve_4 = coalDustValve4.summary.outline.Delta_p,
    Delta_p_coalDustValve_5 = coalDustValve5.summary.outline.Delta_p,
    Delta_p_coalDustValve_6 = coalDustValve6.summary.outline.Delta_p) annotation (Placement(transformation(extent={{162,212},{182,232}})));
equation
  connect(ramp.y, valve1.opening_in) annotation (Line(
      points={{-93,164},{-74,164},{-74,140}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ramp.y, valve2.opening_in) annotation (Line(
      points={{-93,164},{-54,164},{-54,112}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(coalFlowSource.fuel_a,coalGas_join.fuel_inlet)  annotation (Line(
      points={{-126,-150},{-120,-150},{-120,-156},{-114,-156}},
      color={0,0,0},
      pattern=LinePattern.Solid,
      smooth=Smooth.None));
  connect(gasFlowSource_T3.gas_a, coalGas_join.flueGas_inlet) annotation (Line(
      points={{-126,-172},{-120,-172},{-120,-168},{-114,-168}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(coalGas_split.fuel_outlet, coalSink.fuel_a) annotation (Line(
      points={{82,-156},{88,-156},{88,-150},{96,-150}},
      color={0,0,0},
      pattern=LinePattern.Solid,
      smooth=Smooth.None));
  connect(coalGas_split.flueGas_outlet, gasSink_pT3.gas_a) annotation (Line(
      points={{82,-168},{88,-168},{88,-174},{96,-174}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(coalFlowSource1.fuel_a,coalGas_join1.fuel_inlet)  annotation (Line(
      points={{-126,-98},{-120,-98},{-120,-104},{-114,-104}},
      color={0,0,0},
      pattern=LinePattern.Solid,
      smooth=Smooth.None));
  connect(gasFlowSource_T4.gas_a, coalGas_join1.flueGas_inlet) annotation (Line(
      points={{-126,-120},{-120,-120},{-120,-116},{-114,-116}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(coalGas_split1.fuel_outlet, coalSink1.fuel_a) annotation (Line(
      points={{82,-104},{90,-104},{90,-98},{98,-98}},
      color={0,0,0},
      pattern=LinePattern.Solid,
      smooth=Smooth.None));
  connect(coalGas_split1.flueGas_outlet, gasSink_pT4.gas_a) annotation (Line(
      points={{82,-116},{90,-116},{90,-122},{98,-122}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(ramp.y, coalDustValve2.opening_in)
                                            annotation (Line(
      points={{-93,164},{-54,164},{-54,-102}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(coalFlowSource2.fuel_a,coalGas_join2.fuel_inlet)  annotation (Line(
      points={{-126,-50},{-120,-50},{-120,-56},{-114,-56}},
      color={0,0,0},
      pattern=LinePattern.Solid,
      smooth=Smooth.None));
  connect(gasFlowSource_T5.gas_a, coalGas_join2.flueGas_inlet) annotation (Line(
      points={{-126,-72},{-120,-72},{-120,-68},{-114,-68}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(coalGas_split2.fuel_outlet, coalSink2.fuel_a) annotation (Line(
      points={{84,-56},{90,-56},{90,-50},{98,-50}},
      color={0,0,0},
      pattern=LinePattern.Solid,
      smooth=Smooth.None));
  connect(coalGas_split2.flueGas_outlet, gasSink_pT5.gas_a) annotation (Line(
      points={{84,-68},{90,-68},{90,-74},{98,-74}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(ramp.y, coalDustValve1.opening_in) annotation (Line(
      points={{-93,164},{-74,164},{-74,-54}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gasFlowSource_T2.gas_a, valve3.inlet) annotation (Line(
      points={{-94,76},{-70,76},{-70,75},{-44,75}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(valve3.outlet, gasSink_pT2.gas_a) annotation (Line(
      points={{-24,75},{18,75},{18,76},{58,76}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(gasFlowSource_T6.gas_a, valve4.inlet) annotation (Line(
      points={{-94,48},{-60,48},{-60,47},{-24,47}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(valve4.outlet, gasSink_pT6.gas_a) annotation (Line(
      points={{-4,47},{28,47},{28,48},{58,48}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(gasFlowSource_T.gas_a, valve1.inlet) annotation (Line(
      points={{-94,132},{-90,132},{-90,131},{-84,131}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(valve1.outlet, gasSink_pT.gas_a) annotation (Line(
      points={{-64,131},{-2,131},{-2,132},{58,132}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(gasSink_pT1.gas_a, valve2.outlet) annotation (Line(
      points={{58,104},{8,104},{8,103},{-44,103}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(valve2.inlet, gasFlowSource_T1.gas_a) annotation (Line(
      points={{-64,103},{-80,103},{-80,104},{-94,104}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(ramp.y, valve3.opening_in) annotation (Line(
      points={{-93,164},{-34,164},{-34,84}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ramp.y, valve4.opening_in) annotation (Line(
      points={{-93,164},{-14,164},{-14,56}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ramp.y, valve5.opening_in) annotation (Line(
      points={{-93,164},{6,164},{6,28}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gasFlowSource_T7.gas_a, valve5.inlet) annotation (Line(
      points={{-94,20},{-48,20},{-48,19},{-4,19}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(valve5.outlet, gasSink_pT7.gas_a) annotation (Line(
      points={{16,19},{38,19},{38,20},{58,20}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(gasFlowSource_T8.gas_a, valve6.inlet) annotation (Line(
      points={{-94,-8},{-38,-8},{-38,-9},{16,-9}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(valve6.outlet, gasSink_pT8.gas_a) annotation (Line(
      points={{36,-9},{48,-9},{48,-8},{58,-8}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(ramp.y, valve6.opening_in) annotation (Line(
      points={{-93,164},{26,164},{26,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(coalGas_join.fuelFlueGas_outlet, coalDustValve3.inlet) annotation (Line(
      points={{-94,-162},{-70,-162},{-70,-163},{-44,-163}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalDustValve3.outlet, coalGas_split.fuelFlueGas_inlet) annotation (Line(
      points={{-24,-163},{20,-163},{20,-162},{62,-162}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalGas_join2.fuelFlueGas_outlet, coalDustValve1.inlet) annotation (Line(
      points={{-94,-62},{-90,-62},{-90,-63},{-84,-63}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalDustValve1.outlet, coalGas_split2.fuelFlueGas_inlet) annotation (Line(
      points={{-64,-63},{0,-63},{0,-62},{64,-62}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalGas_split1.fuelFlueGas_inlet, coalDustValve2.outlet) annotation (Line(
      points={{62,-110},{10,-110},{10,-111},{-44,-111}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalDustValve2.inlet, coalGas_join1.fuelFlueGas_outlet) annotation (Line(
      points={{-64,-111},{-80,-111},{-80,-110},{-94,-110}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(ramp.y, coalDustValve3.opening_in) annotation (Line(
      points={{-93,164},{-34,164},{-34,-154}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(coalFlowSource3.fuel_a,coalGas_join3.fuel_inlet)
                                                          annotation (Line(
      points={{-126,-254},{-120,-254},{-120,-260},{-114,-260}},
      color={0,0,0},
      pattern=LinePattern.Solid,
      smooth=Smooth.None));
  connect(gasFlowSource_T9.gas_a, coalGas_join3.flueGas_inlet)
                                                              annotation (Line(
      points={{-126,-276},{-120,-276},{-120,-272},{-114,-272}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(coalGas_split3.fuel_outlet, coalSink3.fuel_a) annotation (Line(
      points={{80,-260},{86,-260},{86,-254},{94,-254}},
      color={0,0,0},
      pattern=LinePattern.Solid,
      smooth=Smooth.None));
  connect(coalGas_split3.flueGas_outlet, gasSink_pT9.gas_a) annotation (Line(
      points={{80,-272},{86,-272},{86,-278},{94,-278}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(coalFlowSource4.fuel_a,coalGas_join4.fuel_inlet)  annotation (Line(
      points={{-126,-202},{-120,-202},{-120,-208},{-114,-208}},
      color={0,0,0},
      pattern=LinePattern.Solid,
      smooth=Smooth.None));
  connect(gasFlowSource_T10.gas_a, coalGas_join4.flueGas_inlet)
                                                               annotation (Line(
      points={{-126,-224},{-120,-224},{-120,-220},{-114,-220}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(coalGas_split4.fuel_outlet, coalSink4.fuel_a) annotation (Line(
      points={{80,-208},{88,-208},{88,-202},{96,-202}},
      color={0,0,0},
      pattern=LinePattern.Solid,
      smooth=Smooth.None));
  connect(coalGas_split4.flueGas_outlet, gasSink_pT10.gas_a) annotation (Line(
      points={{80,-220},{88,-220},{88,-226},{96,-226}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(ramp.y, coalDustValve4.opening_in)
                                            annotation (Line(
      points={{-93,164},{-14,164},{-14,-206}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(coalGas_join3.fuelFlueGas_outlet, coalDustValve5.inlet) annotation (Line(
      points={{-94,-266},{-48,-266},{-48,-267},{-4,-267}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalDustValve5.outlet, coalGas_split3.fuelFlueGas_inlet) annotation (Line(
      points={{16,-267},{38,-267},{38,-266},{60,-266}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalGas_split4.fuelFlueGas_inlet, coalDustValve4.outlet) annotation (Line(
      points={{60,-214},{28,-214},{28,-215},{-4,-215}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalDustValve4.inlet, coalGas_join4.fuelFlueGas_outlet) annotation (Line(
      points={{-24,-215},{-60,-215},{-60,-214},{-94,-214}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(ramp.y, coalDustValve5.opening_in) annotation (Line(
      points={{-93,164},{6,164},{6,-258}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(coalFlowSource5.fuel_a,coalGas_join5.fuel_inlet)
                                                          annotation (Line(
      points={{-126,-304},{-120,-304},{-120,-310},{-114,-310}},
      color={0,0,0},
      pattern=LinePattern.Solid,
      smooth=Smooth.None));
  connect(gasFlowSource_T11.gas_a, coalGas_join5.flueGas_inlet)
                                                              annotation (Line(
      points={{-126,-326},{-120,-326},{-120,-322},{-114,-322}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(coalGas_split5.fuel_outlet, coalSink5.fuel_a) annotation (Line(
      points={{80,-310},{86,-310},{86,-304},{94,-304}},
      color={0,0,0},
      pattern=LinePattern.Solid,
      smooth=Smooth.None));
  connect(coalGas_split5.flueGas_outlet, gasSink_pT11.gas_a) annotation (Line(
      points={{80,-322},{86,-322},{86,-328},{94,-328}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(coalGas_join5.fuelFlueGas_outlet, coalDustValve6.inlet) annotation (Line(
      points={{-94,-316},{-38,-316},{-38,-317},{16,-317}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalDustValve6.outlet, coalGas_split5.fuelFlueGas_inlet) annotation (Line(
      points={{36,-317},{48,-317},{48,-316},{60,-316}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(ramp.y, coalDustValve6.opening_in) annotation (Line(
      points={{-93,164},{26,164},{26,-308}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,-420},{200,260}}),
                        graphics={Text(
          extent={{-206,214},{-54,190}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=8,
          textString="_______________________________________
PURPOSE:
>> Tester for gas valves
_______________________________________
"),                    Text(
          extent={{-218,260},{94,228}},
          lineColor={0,128,0},
          fontSize=20,
          textString="TESTED -- 2015-01-22 //LN"),
        Rectangle(
          extent={{-220,260},{200,-420}},
          lineColor={115,150,0},
          lineThickness=0.5)}),
    Icon(coordinateSystem(extent={{-100,-100},{100,100}}, preserveAspectRatio=false)),
    experiment(StopTime=10));
end Test_GasValves;
