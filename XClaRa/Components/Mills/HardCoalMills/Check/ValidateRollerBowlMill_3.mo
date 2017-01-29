within Exergy.XClaRa.Components.Mills.HardCoalMills.Check;
model ValidateRollerBowlMill_3
  "A test scenario derived from the paper Niemczyk: 'Derivation and validation of a coal mill model for control'"
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
  import ModelicaServices.ExternalReferences.loadResource;
//   parameter Real K_10 = 1459.98;//1507.9;
//   parameter Real K_11 = 1.75388e7;//8.4736e6;
//   parameter Real K_2 = 0.127694;//0.1283;
//   parameter Real K_3 = 0.0153057;//0.0149602;
  Real Omega;
  Real Dev1;
  Real Dev2;

model Regression
  extends ClaRa.Basics.Icons.RegressionSummary;
  Modelica.Blocks.Interfaces.RealInput Delta_T_meas_int
      "Temperature deviation sim vs meas.";
  Modelica.Blocks.Interfaces.RealInput T_out "Outlet temperature";
  Modelica.Blocks.Interfaces.RealInput xi_h2o_out
      "H2O mass fraction at air outlet";
  Modelica.Blocks.Interfaces.RealInput m_flow_coal_out "Total coal outlet flow";

  Real y_Delta_T_meas_int = Delta_T_meas_int;

  Real y_xi_h2o_out_int = integrator2.y;

  Real y_T_out_int = integrator3.y;

  Real y_m_flow_coal_out_max = timeExtrema4.y_max;
  Real y_m_flow_coal_out_min = timeExtrema4.y_min;

  protected
  Utilities.Blocks.Integrator integrator2(u=xi_h2o_out);
  Utilities.Blocks.Integrator integrator3(u=T_out);
  Utilities.Blocks.TimeExtrema timeExtrema4(u=m_flow_coal_out);
end Regression;

  Exergy.XClaRa.Components.Mills.HardCoalMills.VerticalMill_L3 Mill2(
    M_c_0=3000,
    initChoice=ClaRa.Basics.Choices.Init.steadyDensity,
    T_0=96.15 + 273.15,
    M_pf_0=100,
    millKoeff=Fundamentals.STV4()) annotation (Placement(
        transformation(extent={{-8,-60},{12,-40}})));
  inner ClaRa.SimCenter simCenter(redeclare ClaRa.Basics.Media.Fuel.Coal_v2
                                      fuelModel1, T_amb=283.15)
    annotation (Placement(transformation(extent={{80,80},{100,100}})));
  Modelica.Blocks.Sources.CombiTimeTable W_c(
    tableOnFile=true,
    tableName="W_c",
    fileName=loadResource("modelica://ClaRa/TableBase/W_c.mat"))
    annotation (Placement(transformation(extent={{-100,-40},{-80,-20}})));
  Modelica.Blocks.Sources.CombiTimeTable
                               T_air(
    tableOnFile=true,
    fileName=loadResource("modelica://ClaRa/TableBase/T_air.mat"),
    tableName="T_air")
    annotation (Placement(transformation(extent={{-100,-100},{-80,-80}})));
  Modelica.Blocks.Sources.CombiTimeTable
                               omega(
    tableOnFile=true,
    tableName="omega",
    fileName=loadResource("modelica://ClaRa/TableBase/omega.mat"))
    annotation (Placement(transformation(extent={{-100,-8},{-80,12}})));
  Modelica.Blocks.Sources.CombiTimeTable
                               W_air(
    tableOnFile=true,
    fileName=loadResource("modelica://ClaRa/TableBase/W_air.mat"),
    tableName="W_air")
    annotation (Placement(transformation(extent={{-100,-72},{-80,-52}})));
  Modelica.Blocks.Sources.CombiTimeTable DeltaP_pa(
    tableOnFile=true,
    fileName=loadResource("modelica://ClaRa/TableBase/DeltaP_pa.mat"),
    tableName="DeltaP_pa")
    annotation (Placement(transformation(extent={{80,52},{100,72}})));
  Modelica.Blocks.Sources.CombiTimeTable DeltaP_mill_meas(
    tableOnFile=true,
    tableName="DeltaP_mill_meas",
    fileName=loadResource("modelica://ClaRa/TableBase/DeltaP_mill_meas.mat"))
    annotation (Placement(transformation(extent={{80,-52},{100,-32}})));
  Modelica.Blocks.Sources.CombiTimeTable
                               E_meas(
    tableOnFile=true,
    tableName="E_meas",
    fileName=loadResource("modelica://ClaRa/TableBase/E_meas.mat"))
    annotation (Placement(transformation(extent={{80,-20},{100,0}})));
  Modelica.Blocks.Sources.CombiTimeTable T_out_meas(
    tableOnFile=true,
    tableName="T_out_meas",
    fileName=loadResource("modelica://ClaRa/TableBase/T_out_meas.mat"),
    offset={-273.15})
    annotation (Placement(transformation(extent={{80,-84},{100,-64}})));
  Modelica.Blocks.Sources.CombiTimeTable
                               E_model(
    tableOnFile=true,
    tableName="E_model",
    fileName=loadResource("modelica://ClaRa/TableBase/E_model.mat"))
    annotation (Placement(transformation(extent={{52,-20},{72,0}})));
  Modelica.Blocks.Sources.CombiTimeTable DeltaP_mill_model(
    tableOnFile=true,
    tableName="DeltaP_mill_model",
    fileName=loadResource("modelica://ClaRa/TableBase/DeltaP_mill_model.mat"))
    annotation (Placement(transformation(extent={{52,-52},{72,-32}})));
  Modelica.Blocks.Sources.CombiTimeTable T_out_model(
    tableOnFile=true,
    offset={-273.15},
    tableName="T_out_model",
    fileName=loadResource("modelica://ClaRa/TableBase/T_out_model.mat"))
    annotation (Placement(transformation(extent={{52,-84},{72,-64}})));
  Modelica.Blocks.Sources.CombiTimeTable W_pf_meas(
    tableOnFile=true,
    tableName="W_pf_meas",
    fileName=loadResource("modelica://ClaRa/TableBase/W_pf_meas.mat"))
    annotation (Placement(transformation(extent={{80,18},{100,38}})));
  Modelica.Blocks.Sources.CombiTimeTable W_pf_model(
    tableOnFile=true,
    tableName="W_pf_model",
    fileName=loadResource("modelica://ClaRa/TableBase/W_pf_model.mat"))
    annotation (Placement(transformation(extent={{52,18},{72,38}})));
  BoundaryConditions.BoundaryFuel_Txim_flow coalFlowSource_XRG1(m_flow_const=10, variable_m_flow=true) annotation (Placement(transformation(extent={{-66,-48},{-46,-28}})));
  Exergy.XClaRa.Components.Adapters.FuelFlueGas_join coalGas_join_burner1
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-24,-50})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow
    fluelGasFlowSource_burner1(
    m_flow_const=11,
    variable_m_flow=true,
    variable_T=true,
    xi_const={0,0,0,0,0.79,0.21,0,0,0}) annotation (Placement(
        transformation(extent={{-64,-82},{-44,-62}})));
  BoundaryConditions.BoundaryFuel_pTxi coaSink_XRG1 annotation (Placement(transformation(
        extent={{9.5,-9.5},{-9.5,9.5}},
        rotation=270,
        origin={-10.5,-90.5})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryGas_pTxi idealGasPressureSink_XRG1(p_const=
        100000) annotation (Placement(transformation(
        extent={{9.5,-9.5},{-9.5,9.5}},
        rotation=270,
        origin={26.5,-91.5})));

  Exergy.XClaRa.Components.Adapters.FuelFlueGas_join coalGas_join_burner2
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={20,-70})));
  Regression regression(T_out = Mill2.summary.T_out,
    m_flow_coal_out = Mill2.summary.m_flow_tot_out,
    xi_h2o_out = Mill2.summary.xi_air_h2o_out,
    Delta_T_meas_int = Dev1) annotation (Placement(transformation(extent={{-100,20},{-80,40}})));
equation
  der(Dev1)=(Mill2.T_out-(T_out_meas.y[1]+273.15))^2;
  der(Dev2)=((Mill2.P_grind - E_meas.y[1])^2)*1e4;
  Omega=Dev1+Dev2;

  connect(omega.y[1], Mill2.classifierSpeed) annotation (Line(
      points={{-79,2},{2,2},{2,-39.2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(coalFlowSource_XRG1.fuel_a,coalGas_join_burner1.fuel_inlet)
    annotation (Line(
      points={{-46,-38},{-40,-38},{-40,-44},{-34,-44}},
      color={0,0,0},
      pattern=LinePattern.Solid,
      smooth=Smooth.None));
  connect(fluelGasFlowSource_burner1.gas_a,coalGas_join_burner1. flueGas_inlet)
    annotation (Line(
      points={{-44,-72},{-40,-72},{-40,-56},{-34,-56}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(coalGas_join_burner1.fuelFlueGas_outlet, Mill2.inlet) annotation (Line(
      points={{-14,-50},{-8,-50}},
      color={175,175,175},
      smooth=Smooth.None));
  connect(fluelGasFlowSource_burner1.m_flow, W_air.y[1]) annotation (Line(
      points={{-64,-66},{-72,-66},{-72,-62},{-79,-62}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(fluelGasFlowSource_burner1.T, T_air.y[1]) annotation (Line(
      points={{-64,-72},{-72,-72},{-72,-90},{-79,-90}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(coalFlowSource_XRG1.m_flow, W_c.y[1]) annotation (Line(
      points={{-66,-32},{-72,-32},{-72,-30},{-79,-30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(idealGasPressureSink_XRG1.gas_a,coalGas_join_burner2. flueGas_inlet)
    annotation (Line(
      points={{26.5,-82},{28,-82},{28,-80},{26,-80}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(coaSink_XRG1.fuel_a,coalGas_join_burner2.fuel_inlet)  annotation (
      Line(
      points={{-10.5,-81},{14,-81},{14,-80}},
      color={127,127,0},
      smooth=Smooth.None));
  connect(coalGas_join_burner2.fuelFlueGas_outlet, Mill2.outlet) annotation (Line(
      points={{20,-60},{20,-50},{12,-50}},
      color={175,175,175},
      smooth=Smooth.None));
  annotation (Diagram(graphics={Text(
          extent={{-96,80},{76,72}},
          lineColor={0,0,255},
          textString=
              "Note, that the Niemczyk model description does not handle coal drying while the model does!"), Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={115,150,0},
          lineThickness=0.5)}),
    experiment(StopTime=5960),
    __Dymola_experimentSetupOutput);
end ValidateRollerBowlMill_3;
