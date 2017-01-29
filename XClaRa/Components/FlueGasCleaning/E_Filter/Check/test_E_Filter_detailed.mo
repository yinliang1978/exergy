within Exergy.XClaRa.Components.FlueGasCleaning.E_Filter.Check;
model test_E_Filter_detailed
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
    variable_m_flow=false,
    m_flow_const=1,
    variable_T=false,
    T_const=473.15,
    xi_const={0.01,0,0.73,0,0.065,0.036,0,0.13,0.0}) annotation (
      Placement(transformation(extent={{-80,-40},{-60,-20}})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryGas_pTxi idealGasPressureSink(p_const=
        100000, xi_const={0.01,0,0.73,0,0.065,0.036,0,0.13,0.0})
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={30,-30})));
  inner ClaRa.SimCenter simCenter(redeclare TILMedia.GasTypes.FlueGasTILMedia
                                        flueGasModel) annotation (
      Placement(transformation(extent={{80,60},{100,80}})));
  E_Filter_L2_detailed e_Filter_dynamic(redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.GenericGeometry,                                                            use_dynamicMassbalance=true) annotation (Placement(transformation(extent={{-44,-40},{-24,-20}})));

  Modelica.Blocks.Sources.Ramp U_applied(
    duration=10,
    height=20e3,
    startTime=10,
    offset=1000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-70,0})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperatureTop1(T=293.15)
                annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-70,32})));
equation
  connect(U_applied.y, e_Filter_dynamic.U_applied) annotation (Line(
      points={{-59,0},{-41.4,0},{-41.4,-18.8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(idealGasFlowSource_XRG.gas_a, e_Filter_dynamic.inlet) annotation (
      Line(
      points={{-60,-30},{-44,-30}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(e_Filter_dynamic.outlet, idealGasPressureSink.gas_a) annotation (Line(
      points={{-24,-30},{20,-30}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(fixedTemperatureTop1.port, e_Filter_dynamic.heat) annotation (Line(
      points={{-60,32},{-34,32},{-34,-20}},
      color={191,0,0},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,
            -50},{100,80}}),
                      graphics={Text(
          extent={{-98,74},{-24,64}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=12,
          textString="________________________________________________________________
PURPOSE:
>>Tester for the E_Filter_detailed component"),
                                Text(
          extent={{-100,80},{-14,70}},
          lineColor={0,128,0},
          fontSize=34,
          textString="TESTED -- 2014-10-08 //LN")}),
    experiment(StopTime=30),
    __Dymola_experimentSetupOutput,
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-50},{100,80}})));
end test_E_Filter_detailed;
