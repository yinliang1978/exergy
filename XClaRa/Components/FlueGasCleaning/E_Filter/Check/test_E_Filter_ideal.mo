within Exergy.XClaRa.Components.FlueGasCleaning.E_Filter.Check;
model test_E_Filter_ideal "with flue gas model that includes argon"
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
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow
    idealGasFlowSource_XRG(
    m_flow_const=10,
    variable_m_flow=true,
    variable_T=true,
    xi_const={0.01,0,0.73,0,0.065,0.036,0,0.13,0.0}) annotation (
      Placement(transformation(extent={{-40,-28},{-20,-8}})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryGas_pTxi idealGasPressureSink(p_const=
        100000, xi_const={0.0,0,0.73,0,0.065,0.036,0,0.13,0.0})
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={70,-18})));
  inner ClaRa.SimCenter simCenter(redeclare TILMedia.GasTypes.FlueGasTILMedia
                                        flueGasModel) annotation (
      Placement(transformation(extent={{80,40},{100,60}})));
  Modelica.Blocks.Sources.Ramp massFlowRate(
    startTime=5,
    duration=1,
    height=-2,
    offset=1)    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-70,2})));
  Modelica.Blocks.Sources.Ramp Temperature(
    duration=1,
    startTime=1,
    height=50,
    offset=273.15 + 150)
                        annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-70,-30})));
  Exergy.XClaRa.Components.FlueGasCleaning.E_Filter.E_Filter_L2_simple
    e_Filter_dynamic(separationRate=0.9995, xi_start={0.0,0,0.73,0,
        0.065,0.036,0,0.13,0.0}) annotation (Placement(
        transformation(extent={{-4,-28},{16,-8}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperatureTop1(T=293.15)
                annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-70,30})));
equation

  connect(massFlowRate.y, idealGasFlowSource_XRG.m_flow) annotation (Line(
      points={{-59,2},{-50,2},{-50,-12},{-40,-12}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(Temperature.y, idealGasFlowSource_XRG.T)
                                             annotation (Line(
      points={{-59,-30},{-50,-30},{-50,-18},{-40,-18}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(idealGasFlowSource_XRG.gas_a, e_Filter_dynamic.inlet) annotation (
      Line(
      points={{-20,-18},{-4,-18}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(e_Filter_dynamic.outlet, idealGasPressureSink.gas_a) annotation (Line(
      points={{16,-18},{60,-18}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(fixedTemperatureTop1.port, e_Filter_dynamic.heat) annotation (Line(
      points={{-60,30},{6,30},{6,-8}},
      color={191,0,0},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,
            -50},{100,60}}),
                      graphics={Text(
          extent={{-98,54},{-24,44}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=12,
          textString="________________________________________________________________
PURPOSE:
>>Tester for the E_Filter_ideal component"),
                                Text(
          extent={{-100,60},{-20,50}},
          lineColor={0,128,0},
          fontSize=34,
          textString="TESTED -- 2014-10-08 //LN")}),
    experiment(StopTime=10),
    __Dymola_experimentSetupOutput,
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-50},{100,60}})));
end test_E_Filter_ideal;
