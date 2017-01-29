within Exergy.XClaRa.Components.Furnace.Check;
model Test_burner_adiabatic_fuelDrying
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
  Adapters.FuelSlagFlueGas_split                                 coalSlagFlueGas_split_top
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-4,30})));
  BoundaryConditions.BoundarySlag_Tm_flow slagFlowSource_top(m_flow_const=0.0, T_const=658.15) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-4,74})));
  BoundaryConditions.BoundaryFuel_pTxi coalSink_top(T_const=658.15) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-30,62})));
  BoundaryConditions.BoundaryGas_pTxi flueGasPressureSink_top(p_const=100000, T_const=658.15) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={22,62})));
  inner ClaRa.SimCenter simCenter(
    redeclare TILMedia.VLEFluidTypes.TILMedia_InterpolatedWater
      fluid1,
    redeclare ClaRa.Basics.Media.Fuel.Coal_Reference fuelModel1,
    redeclare ClaRa.Basics.Media.Fuel.Slag_v2 slagModel,
    redeclare TILMedia.GasTypes.FlueGasTILMedia flueGasModel)
    annotation (Placement(transformation(extent={{98,62},{118,82}})));
  BoundaryConditions.BoundaryFuel_Txim_flow coalFlowSource(
    xi_const={0.661,0.0383,0.066,0.016,0.0057,0.135},
    m_flow_const=15,
    T_const=293.15,
    LHV_predefined=3.0e7,
    LHV_calculationType="Verbandsformel") annotation (Placement(transformation(extent={{-138,-30},{-118,-10}})));
  BoundaryConditions.BoundaryGas_Txim_flow flueGasFlowSource(
    variable_m_flow=false,
    variable_xi=false,
    m_flow_const=15*10.396,
    T_const=573.15,
    xi_const={0,0,0.0005,0,0.7681,0.2314,0,0,0})
                    annotation (Placement(transformation(extent={{-106,-42},{-86,-22}})));
  Adapters.FuelFlueGas_join coalGas_join_burner annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-54,-26})));
  Adapters.FuelSlagFlueGas_join
    coalSlagFlueGas_join
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-4,-54})));
  BoundaryConditions.BoundaryGas_Txim_flow flueGasFlowSourcee_bottom(m_flow_const=0, T_const=283.15) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={16,-94})));
  BoundaryConditions.BoundaryFuel_Txim_flow coalFlowSource_bottom(
    xi_const={0.86,0.035,0.025,0.014,0.007,0.0505},
    m_flow_const=0,
    T_const=293.15,
    LHV_predefined=3.0e7,
    LHV_calculationType="Verbandsformel")
                    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-24,-94})));
  BoundaryConditions.BoundarySlag_pT slagSink_bottom(T_const=373.15) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-4,-94})));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedTemperatureBottom(Q_flow=0)
                annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={66,-50})));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedTemperatureWall(Q_flow=0)
                annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={66,-26})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperatureTop(T=658.15)
                annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={56,8})));
  Burner.Burner_L2_Dynamic burner(
    T_slag=873,
    redeclare model HeatTransfer_Top =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Adiabat_L2,
    redeclare model ReactionZone =
        ClaRa.Components.Furnace.ChemicalReactions.CoalReactionZone (xi_NOx=0,
          xi_CO=0),
    redeclare model ParticleMigration =
        ClaRa.Components.Furnace.GeneralTransportPhenomena.ParticleMigration.MeanMigrationSpeed,
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.HollowBlock (
        z_in={0},
        z_out={5},
        width=10,
        length=burner.geo.z_out[1] - burner.geo.z_in[1],
        height=10),
    redeclare model Burning_time =
        ClaRa.Components.Furnace.GeneralTransportPhenomena.BurningTime.ConstantBurningTime)
    annotation (Placement(transformation(extent={{-18,-36},{42,-16}})));

  Adapters.FuelSlagFlueGas_split                                 coalSlagFlueGas_split_top1
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-2,-164})));
  BoundaryConditions.BoundarySlag_Tm_flow slagFlowSource_top1(
                                                             m_flow_const=0.0, T_const=658.15) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-2,-120})));
  BoundaryConditions.BoundaryFuel_pTxi coalSink_top1(
                                                    T_const=658.15) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-28,-132})));
  BoundaryConditions.BoundaryGas_pTxi flueGasPressureSink_top1(
                                                              p_const=100000, T_const=658.15) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={24,-132})));
  BoundaryConditions.BoundaryFuel_Txim_flow coalFlowSource1(
    m_flow_const=15,
    xi_const={0.661,0.0383,0.066,0.05,0.0057,0.135},
    T_const=298.15,
    LHV_predefined=3.0e7,
    LHV_calculationType="Verbandsformel") annotation (Placement(transformation(extent={{-136,-224},{-116,-204}})));
  BoundaryConditions.BoundaryGas_Txim_flow flueGasFlowSource1(
    variable_m_flow=false,
    variable_xi=false,
    m_flow_const=15*10.396,
    T_const=573.15,
    xi_const={0,0,0.0005,0,0.7681,0.2314,0,0,0})
                    annotation (Placement(transformation(extent={{-104,-236},{-84,-216}})));
  Adapters.FuelFlueGas_join coalGas_join_burner1 annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-52,-220})));
  Adapters.FuelSlagFlueGas_join
    coalSlagFlueGas_join1
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-2,-248})));
  BoundaryConditions.BoundaryGas_Txim_flow flueGasFlowSourcee_bottom1(
                                                                     m_flow_const=0, T_const=283.15) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={18,-288})));
  BoundaryConditions.BoundaryFuel_Txim_flow coalFlowSource_bottom1(
    xi_const={0.86,0.035,0.025,0.014,0.007,0.0505},
    m_flow_const=0,
    T_const=293.15,
    LHV_predefined=3.0e7,
    LHV_calculationType="Verbandsformel") annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-22,-288})));
  BoundaryConditions.BoundarySlag_pT slagSink_bottom1(
                                                     T_const=373.15) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-2,-288})));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedTemperatureBottom1(
                                                                             Q_flow=0)
                annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={68,-244})));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedTemperatureWall1(
                                                                           Q_flow=0)
                annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={68,-220})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperatureTop1(
                                                                             T=658.15)
                annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={58,-186})));
  Burner.Burner_L2_Dynamic_fuelDrying burner_drying(
    T_slag=873,
    redeclare model HeatTransfer_Top =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Adiabat_L2,
    redeclare model ReactionZone =
        ClaRa.Components.Furnace.ChemicalReactions.CoalReactionZone (                           xi_NOx=0, xi_CO=0),
    redeclare model ParticleMigration =
        ClaRa.Components.Furnace.GeneralTransportPhenomena.ParticleMigration.MeanMigrationSpeed,
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.HollowBlock (
        z_in={0},
        z_out={5},
        width=10,
        length=burner.geo.z_out[1] - burner.geo.z_in[1],
        height=10),
    redeclare model Burning_time =
        ClaRa.Components.Furnace.GeneralTransportPhenomena.BurningTime.ConstantBurningTime)
                                                                                                        annotation (Placement(transformation(extent={{-16,-230},{44,-210}})));
equation
  connect(coalFlowSource.fuel_a,coalGas_join_burner.fuel_inlet)  annotation (
      Line(
      points={{-118,-20},{-64,-20}},
      color={27,36,42},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));
  connect(flueGasFlowSource.gas_a, coalGas_join_burner.flueGas_inlet)
    annotation (Line(
      points={{-86,-32},{-64,-32}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalGas_join_burner.fuelFlueGas_outlet, burner.fuelFlueGas_inlet) annotation (Line(
      points={{-44,-26},{-18,-26}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalSlagFlueGas_split_top.fuelSlagFlueGas_inlet, burner.outlet) annotation (Line(
      points={{-4,20},{-4,-16}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalSlagFlueGas_join.fuelSlagFlueGas_outlet, burner.inlet)
    annotation (Line(
      points={{-4,-44},{-4,-36}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalSlagFlueGas_join.flueGas_inlet, flueGasFlowSourcee_bottom.gas_a)
    annotation (Line(
      points={{2,-64},{2,-74},{16,-74},{16,-84}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(slagSink_bottom.slag_inlet, coalSlagFlueGas_join.slag_outlet)
    annotation (Line(
      points={{-3.8,-84},{-3.8,-74},{-4,-74},{-4,-64}},
      color={234,171,0},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalSlagFlueGas_join.fuel_inlet,coalFlowSource_bottom.fuel_a)
    annotation (Line(
      points={{-10,-64},{-10,-74},{-24,-74},{-24,-84}},
      color={27,36,42},
      thickness=0.5,
      smooth=Smooth.None));
  connect(burner.heat_bottom, fixedTemperatureBottom.port) annotation (Line(
      points={{14,-36},{14,-50},{56,-50}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(burner.heat_wall, fixedTemperatureWall.port) annotation (Line(
      points={{42,-26},{56,-26}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(burner.heat_top, fixedTemperatureTop.port) annotation (Line(
      points={{14,-16},{14,8},{46,8}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalSink_top.fuel_a, coalSlagFlueGas_split_top.fuel_outlet) annotation (Line(
      points={{-30,52},{-30,48},{-10,48},{-10,40}},
      color={27,36,42},
      thickness=0.5,
      smooth=Smooth.None));
  connect(flueGasPressureSink_top.gas_a, coalSlagFlueGas_split_top.flueGas_outlet) annotation (Line(
      points={{22,52},{22,48},{2,48},{2,40}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalSlagFlueGas_split_top.slag_inlet, slagFlowSource_top.slag_outlet) annotation (Line(
      points={{-4,40},{-4,64}},
      color={234,171,0},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalFlowSource1.fuel_a,coalGas_join_burner1.fuel_inlet)  annotation (Line(
      points={{-116,-214},{-62,-214}},
      color={27,36,42},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));
  connect(flueGasFlowSource1.gas_a, coalGas_join_burner1.flueGas_inlet) annotation (Line(
      points={{-84,-226},{-62,-226}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalGas_join_burner1.fuelFlueGas_outlet, burner_drying.fuelFlueGas_inlet) annotation (Line(
      points={{-42,-220},{-16,-220}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalSlagFlueGas_split_top1.fuelSlagFlueGas_inlet, burner_drying.outlet) annotation (Line(
      points={{-2,-174},{-2,-210}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalSlagFlueGas_join1.fuelSlagFlueGas_outlet, burner_drying.inlet) annotation (Line(
      points={{-2,-238},{-2,-230}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalSlagFlueGas_join1.flueGas_inlet, flueGasFlowSourcee_bottom1.gas_a) annotation (Line(
      points={{4,-258},{4,-268},{18,-268},{18,-278}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(slagSink_bottom1.slag_inlet, coalSlagFlueGas_join1.slag_outlet) annotation (Line(
      points={{-1.8,-278},{-1.8,-268},{-2,-268},{-2,-258}},
      color={234,171,0},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalSlagFlueGas_join1.fuel_inlet,coalFlowSource_bottom1.fuel_a)  annotation (Line(
      points={{-8,-258},{-8,-268},{-22,-268},{-22,-278}},
      color={27,36,42},
      thickness=0.5,
      smooth=Smooth.None));
  connect(burner_drying.heat_bottom, fixedTemperatureBottom1.port) annotation (Line(
      points={{16,-230},{16,-244},{58,-244}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(burner_drying.heat_wall, fixedTemperatureWall1.port) annotation (Line(
      points={{44,-220},{58,-220}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(burner_drying.heat_top, fixedTemperatureTop1.port) annotation (Line(
      points={{16,-210},{16,-186},{48,-186}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalSink_top1.fuel_a, coalSlagFlueGas_split_top1.fuel_outlet) annotation (Line(
      points={{-28,-142},{-28,-146},{-8,-146},{-8,-154}},
      color={27,36,42},
      thickness=0.5,
      smooth=Smooth.None));
  connect(flueGasPressureSink_top1.gas_a, coalSlagFlueGas_split_top1.flueGas_outlet) annotation (Line(
      points={{24,-142},{24,-146},{4,-146},{4,-154}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(coalSlagFlueGas_split_top1.slag_inlet, slagFlowSource_top1.slag_outlet) annotation (Line(
      points={{-2,-154},{-2,-130}},
      color={234,171,0},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(extent={{-140,-300},{120,100}},
          preserveAspectRatio=false), graphics={
                                Text(
          extent={{-138,86},{-54,40}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=12,
          textString="________________________________________________________________
PURPOSE:
>>Tester for calculation of the adiabatic outlet temperature of a burner level

"),                             Text(
          extent={{-142,100},{116,78}},
          lineColor={0,128,0},
          fontSize=34,
          textString="TESTED -- 2016-03-09 //LN")}),                     Icon(
        coordinateSystem(extent={{-100,-100},{100,100}}, preserveAspectRatio=
            false)));
end Test_burner_adiabatic_fuelDrying;
