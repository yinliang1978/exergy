within Exergy.XClaRa.Components.FlueGasCleaning.Denitrification.Check;
model Test_Denitrification_NH3port
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
  inner ClaRa.SimCenter simCenter(redeclare TILMedia.GasTypes.FlueGasTILMedia
                                        flueGasModel) annotation (
      Placement(transformation(extent={{80,80},{100,100}})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow
    idealGasFlowSource_XRG1(
    m_flow_const=10,
    variable_m_flow=true,
    variable_T=true,
    medium=simCenter.flueGasModel,
    xi_const={0.01,0.01,0.73,0.01,0.065,0.036,0.01,0.13,0.0})
    annotation (Placement(transformation(extent={{-32,-46},{-12,-26}})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryGas_pTxi idealGasPressureSink(medium=
        simCenter.flueGasModel, p_const=100000) annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={52,-36})));
  Modelica.Blocks.Sources.Ramp massFlowRate1(
    offset=1e-3,
    height=9e-3,
    startTime=5,
    duration=5)  annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-70,-18})));
  Modelica.Blocks.Sources.Ramp Temperature1(
    duration=1,
    height=20,
    startTime=1,
    offset=273.15 + 250)
                        annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-70,-50})));
  Denitrification_L1_NH3port deNOx(redeclare model HeatTransfer =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Adiabat_L2)
    annotation (Placement(transformation(extent={{4,-46},{24,-26}})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow
    idealGasFlowSource_XRG2(
    variable_T=false,
    m_flow_const=0.0010,
    variable_m_flow=true,
    medium=simCenter.flueGasModel,
    T_const=523.15,
    xi_const={0/9,0/9,0/9,0/9,0/9,0/9,0/9,0/9,9/9}) annotation (
      Placement(transformation(extent={{-32,-14},{-12,6}})));
  Modelica.Blocks.Sources.Ramp massFlowRate2(
    offset=0.5e-4,
    height=5e-3,
    startTime=5,
    duration=30) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-70,14})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperatureTop1(T=293.15)
                annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-70,44})));
equation
  connect(massFlowRate1.y, idealGasFlowSource_XRG1.m_flow) annotation (Line(
      points={{-59,-18},{-50,-18},{-50,-30},{-32,-30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(Temperature1.y, idealGasFlowSource_XRG1.T)
                                             annotation (Line(
      points={{-59,-50},{-50,-50},{-50,-36},{-32,-36}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowRate2.y, idealGasFlowSource_XRG2.m_flow) annotation (Line(
      points={{-59,14},{-46,14},{-46,2},{-32,2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(idealGasFlowSource_XRG2.gas_a, deNOx.NH3_inlet) annotation (Line(
      points={{-12,-4},{14,-4},{14,-26}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(idealGasFlowSource_XRG1.gas_a, deNOx.inlet) annotation (Line(
      points={{-12,-36},{4,-36}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(deNOx.outlet, idealGasPressureSink.gas_a) annotation (Line(
      points={{24,-36},{42,-36}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(fixedTemperatureTop1.port, deNOx.heat) annotation (Line(
      points={{-60,44},{8.8,44},{8.8,-26.4}},
      color={191,0,0},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,
            -100},{100,100}}),
                      graphics={Text(
          extent={{-98,92},{-24,82}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=12,
          textString="________________________________________________________________
PURPOSE:
>>Tester for the Denitrificationl component"),
                                Text(
          extent={{-100,98},{26,88}},
          lineColor={0,128,0},
          fontSize=34,
          textString="TESTED -- 2014-10-08 //LN")}),
    experiment(StopTime=20),
    __Dymola_experimentSetupOutput,
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}})));
end Test_Denitrification_NH3port;
