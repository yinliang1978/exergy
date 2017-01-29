within Exergy.XClaRa.Components.HeatExchangers.Check;
model Test_HEXvle2gas_L3_1ph_BU
  "Example 1 at page Ca 15 in VDI Waermeatlas, 9th edition "
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
  //{0.001,0.001,0.001,0.001,0.79,0.2,0.004,0.001,0} composition for flueGasModel moist air mixture with a condensing component
  BoundaryConditions.BoundaryVLE_Txim_flow massFlowSource_T(
    m_flow_const=1,
    variable_m_flow=false,
    T_const=120 + 273.15) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={68,-14})));
  BoundaryConditions.BoundaryVLE_pTxi pressureSink_pT(T_const=303.15, p_const=11e5) annotation (Placement(transformation(extent={{100,0},{80,20}})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryGas_pTxi flueGasPressureSink(p_const=
        1e5) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-50,-30})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow flueGasFlowSource(
    variable_m_flow=false,
    variable_T=false,
    T_const(displayUnit="degC") = 293.15,
    m_flow_const=2*TILMedia.GasFunctions.density_pTxi(
                  flueGasFlowSource.medium,
                  1e5,
                  273.15 + 20),
    variable_xi=true)
    annotation (Placement(transformation(extent={{-60,0},{-40,20}})));
  inner ClaRa.SimCenter simCenter(
    useHomotopy=false,
    redeclare replaceable TILMedia.VLEFluidTypes.TILMedia_InterpolatedWater
                                                        fluid1,
    showExpertSummary=true) annotation (Placement(transformation(
          extent={{-100,-180},{-60,-160}})));
  HEXvle2gas_L3_1ph_BU_ntu hex_ntu(
    length=1,
    height=1,
    width=1,
    m_nom1=2,
    h_nom1=24000,
    diameter_i=12e-3,
    diameter_o=16e-3,
    m_nom2=1,
    N_passes=6,
    showExpertSummary=true,
    mass_struc=100,
    N_tubes=120,
    p_nom1=100000,
    p_start_shell=1000000,
    p_nom2=1000000,
    redeclare model PressureLossShell =
        ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.NoFriction_L2,
    redeclare model PressureLossTubes =
        ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L2,
    redeclare model HeatTransfer_Shell =
        Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.CharLine_L2
        (                                                                                                    alpha_nom=48.7, PL_alpha=[0,1; 1,1]),
    redeclare model HeatExchangerType =
        Basics.ControlVolumes.SolidVolumes.Fundamentals.HeatExchangerTypes.CrossCounterFlow
        (                                                                                                    N_rp=6),
    h_nom2=500e3,
    T_start_shell=20 + 273.15,
    initTypeShell=ClaRa.Basics.Choices.Init.noInit,
    h_start_tubes=375e3,
    initWall=ClaRa.Basics.Choices.Init.noInit,
    p_start_tubes=10e5,
    initTypeTubes=ClaRa.Basics.Choices.Init.noInit,
    T_w_i_start=343,
    T_w_a_start=343,
    redeclare model HeatTransferTubes =
        Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Constant_L2
        (                                                                                                    alpha_nom=1000),
    CF_geo=3) annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=180,
        origin={8,-8})));

  BoundaryConditions.BoundaryVLE_Txim_flow massFlowSource_T1(
    m_flow_const=1,
    variable_m_flow=false,
    T_const=120 + 273.15) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={68,-74})));
  BoundaryConditions.BoundaryVLE_pTxi pressureSink_pT1(
                                                      T_const=303.15, p_const=10e5) annotation (Placement(transformation(extent={{78,-60},{58,-40}})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryGas_pTxi flueGasPressureSink1(p_const=
        1e5) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-50,-90})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow flueGasFlowSource1(
    variable_m_flow=false,
    variable_T=false,
    T_const(displayUnit="degC") = 293.15,
    m_flow_const=2*TILMedia.GasFunctions.density_pTxi(
                  flueGasFlowSource1.medium,
                  1e5,
                  273.15 + 20),
    variable_xi=true) annotation (Placement(transformation(extent={{-60,
            -60},{-40,-40}})));
  HEXvle2gas_L3_1ph_BU_simple hex_simple(
    length=1,
    height=1,
    width=1,
    m_nom1=2,
    h_nom1=24000,
    diameter_i=12e-3,
    diameter_o=16e-3,
    m_nom2=1,
    h_nom2=1200e3,
    N_passes=6,
    T_w_i_start=353,
    T_w_a_start=333,
    showExpertSummary=true,
    N_tubes=120,
    p_nom1=100000,
    p_start_shell=1000000,
    p_nom2=1000000,
    initTypeShell=ClaRa.Basics.Choices.Init.noInit,
    redeclare model HeatTransfer_Shell =
        Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.CharLine_L2
        (                                                                                                    alpha_nom=48.7, PL_alpha=[0,1; 1,1]),
    p_start_tubes=10e5,
    initWall=ClaRa.Basics.Choices.Init.noInit,
    initTypeTubes=ClaRa.Basics.Choices.Init.steadyState,
    redeclare model PressureLossShell =
        Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.QuadraticNominalPoint_L2,
    CF_geo=3,
    T_start_shell=273.15 + 20,
    h_start_tubes=400e3,
    redeclare model HeatTransferTubes =
        Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Constant_L2
        (                                                                                                    alpha_nom=1000),
    redeclare model PressureLossTubes =
        Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.NoFriction_L2)
                                                         annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=180,
        origin={8,-68})));
  ClaRa.Visualisation.Quadruple quadruple1(largeFonts=false)
    annotation (Placement(transformation(extent={{22,0},{44,12}})));
  ClaRa.Visualisation.Quadruple quadruple2(largeFonts=false)
    annotation (Placement(transformation(extent={{20,-60},{44,-48}})));
  HEXvle2gas_L3_2ph_BU_simple hex_2ph(
    initTypeTubes=ClaRa.Basics.Choices.Init.steadyDensity,
    p_start_tubes=11e5,
    initWall=ClaRa.Basics.Choices.Init.noInit,
    flowOrientation=ClaRa.Basics.Choices.GeometryOrientation.vertical,
    T_start_shell=273.15 + 20,
    h_liq_start=400e3,
    h_vap_start=400e3,
    m_nom2=1,
    length=1,
    height=1,
    width=1,
    diameter_i=12e-3,
    diameter_o=16e-3,
    N_tubes=120,
    N_passes=6,
    CF_geo=3,
    redeclare model PressureLossTubes =
        Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.NoFriction_L3,
    redeclare model HeatTransfer_Shell =
        Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.CharLine_L2
        (                                                                                                    alpha_nom=48.7, PL_alpha=[0,1; 1,1]),
    redeclare model PressureLossShell =
        Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.QuadraticNominalPoint_L2,
    m_nom1=2,
    p_nom1=1e5,
    h_nom1=24e3,
    initTypeShell=ClaRa.Basics.Choices.Init.noInit,
    T_w_i_start=273.15 + 50,
    T_w_a_start=273.15 + 50,
    tubes(Tau_cond=0.3, Tau_evap=0.03),
    redeclare model HeatTransferTubes =
        Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Constant_L3
        (                                                                                                    alpha_nom={1000,1000}),
    z_out_tubes=0,
    level_rel_start=0.95)
                   annotation (Placement(transformation(extent={{0,-140},{20,-120}})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryGas_pTxi flueGasPressureSink2(p_const=
        1e5) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-50,-150})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow flueGasFlowSource2(
    variable_m_flow=false,
    variable_T=false,
    T_const(displayUnit="degC") = 293.15,
    m_flow_const=2*TILMedia.GasFunctions.density_pTxi(
                  flueGasFlowSource1.medium,
                  1e5,
                  273.15 + 20),
    variable_xi=true) annotation (Placement(transformation(extent={{-60,
            -120},{-40,-100}})));
  BoundaryConditions.BoundaryVLE_pTxi pressureSink_pT2(
                                                      T_const=303.15, p_const=10e5) annotation (Placement(transformation(extent={{100,-120},{80,-100}})));
  BoundaryConditions.BoundaryVLE_Txim_flow massFlowSource_T2(
    variable_m_flow=false,
    T_const=120 + 273.15,
    m_flow_const=1)       annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={70,-134})));
  ClaRa.Visualisation.Quadruple quadruple3(largeFonts=false)
    annotation (Placement(transformation(extent={{24,-122},{46,-110}})));
  VolumesValvesFittings.Valves.ValveVLE_L1 valveVLE_L1_1(redeclare model
      PressureLoss =
        VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (                                       m_flow_nom=10, Delta_p_nom=100))
                                                         annotation (Placement(transformation(extent={{56,-116},{76,-104}})));
  Sensors.GasTemperatureSensor gasTemperatureSensor_hex_ntu annotation (Placement(transformation(extent={{-10,-30},{-30,-8}})));
  Sensors.GasTemperatureSensor gasTemperatureSensor_hex_2ph annotation (Placement(transformation(extent={{-10,-150},{-30,-128}})));
  Sensors.GasTemperatureSensor gasTemperatureSensor_hex_simple annotation (Placement(transformation(extent={{-10,-90},{-30,-68}})));
  VolumesValvesFittings.Valves.ValveVLE_L1 valveVLE_L1_2(redeclare model
      PressureLoss =
        VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (                                       m_flow_nom=10))
                                                         annotation (Placement(transformation(extent={{56,4},{76,16}})));
  BoundaryConditions.GasCompositionByMassFractions
    gasCompositionByMassFractions(
    xi_ASH=0,
    xi_CO=0,
    xi_CO2=0,
    xi_SO2=0,
    xi_N2=0.79,
    xi_O2=0.21,
    xi_NO=0,
    xi_H2O=0,
    xi_NH3=0) annotation (Placement(transformation(extent={{-94,-6},{-74,14}})));
equation
  connect(hex_ntu.In2, massFlowSource_T.steam_a) annotation (Line(
      points={{18,-14},{58,-14}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(flueGasFlowSource.gas_a, hex_ntu.In1) annotation (Line(
      points={{-40,10},{8,10},{8,1.8}},
      color={118,106,98},
      thickness=0.5));
  connect(hex_simple.In2, massFlowSource_T1.steam_a) annotation (Line(
      points={{18,-74},{58,-74}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(pressureSink_pT1.steam_a, hex_simple.Out2) annotation (Line(
      points={{58,-50},{46,-50},{46,-62},{18,-62}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(flueGasFlowSource1.gas_a, hex_simple.In1) annotation (Line(
      points={{-40,-50},{8,-50},{8,-58.2}},
      color={118,106,98},
      thickness=0.5));
  connect(hex_ntu.eye, quadruple1.eye) annotation (Line(points={{19,-1.77636e-015},{20,-1.77636e-015},{20,6},{22,6}}, color={190,190,190}));
  connect(hex_simple.eye, quadruple2.eye) annotation (Line(points={{18,-60},{18,-60},{20,-60},{20,-56},{20,-54}},          color={190,190,190}));
  connect(flueGasFlowSource2.gas_a, hex_2ph.In1) annotation (Line(
      points={{-40,-110},{-40,-110},{10,-110},{10,-120.2}},
      color={118,106,98},
      thickness=0.5));
  connect(massFlowSource_T2.steam_a, hex_2ph.In2) annotation (Line(
      points={{60,-134},{20,-134}},
      color={0,131,169},
      thickness=0.5));
  connect(hex_2ph.eye, quadruple3.eye) annotation (Line(points={{20,-122},{22,-122},{22,-116},{24,-116}}, color={190,190,190}));
  connect(hex_2ph.Out2, valveVLE_L1_1.inlet) annotation (Line(
      points={{19.8,-124},{32,-124},{52,-124},{52,-110},{56,-110}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(valveVLE_L1_1.outlet, pressureSink_pT2.steam_a) annotation (Line(
      points={{76,-110},{80,-110}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(hex_2ph.Out1, gasTemperatureSensor_hex_2ph.inlet) annotation (Line(
      points={{10,-140},{10,-140},{10,-150},{-10,-150}},
      color={118,106,98},
      thickness=0.5));
  connect(gasTemperatureSensor_hex_2ph.outlet, flueGasPressureSink2.gas_a) annotation (Line(
      points={{-30,-150},{-40,-150}},
      color={118,106,98},
      thickness=0.5));
  connect(hex_simple.Out1, gasTemperatureSensor_hex_simple.inlet) annotation (Line(
      points={{8,-78},{8,-90},{-10,-90}},
      color={118,106,98},
      thickness=0.5));
  connect(gasTemperatureSensor_hex_simple.outlet, flueGasPressureSink1.gas_a) annotation (Line(
      points={{-30,-90},{-35,-90},{-40,-90}},
      color={118,106,98},
      thickness=0.5));
  connect(hex_ntu.Out1, gasTemperatureSensor_hex_ntu.inlet) annotation (Line(
      points={{8,-18},{8,-18},{8,-30},{-10,-30}},
      color={118,106,98},
      thickness=0.5));
  connect(gasTemperatureSensor_hex_ntu.outlet, flueGasPressureSink.gas_a) annotation (Line(
      points={{-30,-30},{-30,-30},{-40,-30}},
      color={118,106,98},
      thickness=0.5));
  connect(valveVLE_L1_2.outlet, pressureSink_pT.steam_a) annotation (Line(
      points={{76,10},{78,10},{80,10}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(hex_ntu.Out2, valveVLE_L1_2.inlet) annotation (Line(
      points={{18,-2},{30,-2},{52,-2},{52,10},{56,10}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(gasCompositionByMassFractions.X, flueGasFlowSource.xi) annotation (Line(points={{-72,4},{-64,4},{-60,4}}, color={0,0,127}));
  connect(gasCompositionByMassFractions.X, flueGasFlowSource1.xi) annotation (Line(points={{-72,4},{-70,4},{-70,2},{-68,2},{-68,-56},{-60,-56}}, color={0,0,127}));
  connect(gasCompositionByMassFractions.X, flueGasFlowSource2.xi) annotation (Line(points={{-72,4},{-68,4},{-68,-116},{-60,-116}}, color={0,0,127}));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-180},{100,100}}),
            graphics={                               Text(
          extent={{-100,100},{90,40}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=10,
          textString="
_________________________________________________________________________________________________________          
PURPOSE:
Compare model with VDI Waermeatlas, Ca15 Example 1
_________________________________________________________________________________________________________  
NOTE: 
> The fin geometry is modeled by tuning a correction value for the geometry, namely CF_geo such that the simulation fits
   the product of area of heat transfer and heat transfer coeffictient (US).
_________________________________________________________________________________________________________  
LOOK AT:
> The outlet temperatures of gas and water flow.
> Outlet temperature of water according to literature: T_h2o_out = 78 °C
> Outlet temperature of air according to literature: T_air_out =94 °C.
> The simulation meets these data quite well.
_________________________________________________________________________________________________________    
")}),
    Icon(coordinateSystem(initialScale=0.1),
         graphics),
    experiment(
      StopTime=5000,
      Tolerance=1e-005,
      __Dymola_Algorithm="Dassl"),
    __Dymola_experimentSetupOutput);
end Test_HEXvle2gas_L3_1ph_BU;
