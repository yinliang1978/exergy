within Exergy.XClaRa.Components.Sensors.Check;
model testGasComponentSensor
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
  GasComponentSensor SensorCO2(compositionDefinedBy=1, component=3)
    annotation (Placement(transformation(extent={{-38,-20},{-18,0}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT(                  variable_p=true,
    variable_xi=false,
    xi_const={0,0,0.4,0,0.4,0.2,0,0,0})                                             annotation (Placement(transformation(extent={{-68,-30},{-48,-10}})));
  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T(m_flow_const=-10) annotation (Placement(transformation(extent={{98,-30},{78,-10}})));
  inner ClaRa.SimCenter simCenter annotation (Placement(
        transformation(extent={{-100,-140},{-80,-120}})));
  Modelica.Blocks.Sources.Sine sine(
    freqHz=0.5,
    offset=100000,
    amplitude=20000,
    phase=0.017453292519943)
    annotation (Placement(transformation(extent={{-96,-8},{-76,12}})));
  ClaRa.Basics.ControlVolumes.GasVolumes.VolumeGas_L2 flueGasCell
    annotation (Placement(transformation(extent={{50,-30},{70,-10}})));
  GasComponentSensor SensorN2(compositionDefinedBy=1, component=5)
    annotation (Placement(transformation(extent={{-4,-20},{16,0}})));
  GasComponentSensor SensorO2(compositionDefinedBy=1, component=6)
    annotation (Placement(transformation(extent={{24,-20},{44,0}})));
  GasComponentSensor SensorCO1(component=3, compositionDefinedBy=2)
    annotation (Placement(transformation(extent={{-36,-86},{-16,-66}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT1(                  variable_p=true,
    variable_xi=false,
    xi_const={0,0,0.4,0,0.4,0.2,0,0,0})                                              annotation (Placement(transformation(extent={{-66,-96},{-46,-76}})));
  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T1(m_flow_const=-10) annotation (Placement(transformation(extent={{98,-96},{78,-76}})));
  Modelica.Blocks.Sources.Sine sine1(
    freqHz=0.5,
    offset=100000,
    amplitude=20000,
    phase=0.017453292519943)
    annotation (Placement(transformation(extent={{-96,-74},{-76,-54}})));
  ClaRa.Basics.ControlVolumes.GasVolumes.VolumeGas_L2 flueGasCell1
    annotation (Placement(transformation(extent={{50,-96},{70,-76}})));
  GasComponentSensor SensorN1(component=5, compositionDefinedBy=2)
    annotation (Placement(transformation(extent={{-12,-86},{8,-66}})));
  GasComponentSensor SensorO1(component=6, compositionDefinedBy=2)
    annotation (Placement(transformation(extent={{14,-86},{34,-66}})));
equation
  connect(sine.y, gasSink_pT.p) annotation (Line(
      points={{-75,2},{-72,2},{-72,-14},{-68,-14}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(flueGasCell.outlet, gasFlowSource_T.gas_a) annotation (Line(
      points={{70,-20},{78,-20}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(gasSink_pT.gas_a, SensorCO2.inlet)         annotation (Line(
      points={{-48,-20},{-38,-20}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(SensorO2.outlet, flueGasCell.inlet) annotation (Line(
      points={{44,-20},{50,-20}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(sine1.y, gasSink_pT1.p)
                                annotation (Line(
      points={{-75,-64},{-72,-64},{-72,-80},{-66,-80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(flueGasCell1.outlet, gasFlowSource_T1.gas_a)
                                                     annotation (Line(
      points={{70,-86},{78,-86}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(gasSink_pT1.gas_a, SensorCO1.inlet)        annotation (Line(
      points={{-46,-86},{-36,-86}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(SensorO1.outlet, flueGasCell1.inlet) annotation (Line(
      points={{34,-86},{50,-86}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(SensorO2.inlet, SensorN2.outlet) annotation (Line(
      points={{24,-20},{16,-20}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(SensorN2.inlet, SensorCO2.outlet) annotation (Line(
      points={{-4,-20},{-18,-20}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(SensorN1.inlet, SensorCO1.outlet) annotation (Line(
      points={{-12,-86},{-16,-86}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(SensorO1.inlet, SensorN1.outlet) annotation (Line(
      points={{14,-86},{8,-86}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (
    Diagram(coordinateSystem(extent={{-100,-140},{100,100}},
          preserveAspectRatio=false),
            graphics={            Text(
          extent={{-96,94},{102,54}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=9,
          textString="______________________________________________________________________________________________
PURPOSE:

______________________________________________________________________________________________
"),                    Text(
          extent={{-112,100},{88,80}},
          lineColor={0,128,0},
          fontSize=31,
          textString="TESTED -- YYYY-MM-DD //XX"),Text(
          extent={{-96,52},{68,38}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=8,
          textString="______________________________________________________________________________________________________________
Remarks: 
______________________________________________________________________________________________________________
"),                   Text(
          extent={{-96,64},{120,52}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=10,
          textString="______________________________________________________________________________________________
Scenario:  

______________________________________________________________________________________________
")}),
    experiment(StopTime=100),
    __Dymola_experimentSetupOutput,
    Icon(coordinateSystem(extent={{-100,-100},{100,100}}, preserveAspectRatio=
            false)));
end testGasComponentSensor;
