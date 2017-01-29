within Exergy.XClaRa.Components.Sensors.Check;
model testGasCompositionSensor
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
  GasCompositionSensor
                     SensorCO2(compositionDefinedBy=2, N=8)
    annotation (Placement(transformation(extent={{-10,-20},{10,0}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT(variable_p=true, variable_xi=false) annotation (Placement(transformation(extent={{-46,-30},{-26,-10}})));
  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T(m_flow_const=-10) annotation (Placement(transformation(extent={{80,-30},{60,-10}})));
  inner ClaRa.SimCenter simCenter annotation (Placement(
        transformation(extent={{-100,-120},{-80,-100}})));
  Modelica.Blocks.Sources.Sine sine(
    freqHz=0.5,
    offset=100000,
    amplitude=20000,
    phase=0.017453292519943)
    annotation (Placement(transformation(extent={{-82,-24},{-62,-4}})));
  ClaRa.Basics.ControlVolumes.GasVolumes.VolumeGas_L2 flueGasCell
    annotation (Placement(transformation(extent={{24,-30},{44,-10}})));
  GasCompositionSensor
                     SensorCO1(compositionDefinedBy=1, N=8)
    annotation (Placement(transformation(extent={{-8,-60},{12,-40}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT1(variable_p=true, variable_xi=false) annotation (Placement(transformation(extent={{-46,-70},{-26,-50}})));
  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T1(m_flow_const=-10) annotation (Placement(transformation(extent={{80,-70},{60,-50}})));
  Modelica.Blocks.Sources.Sine sine1(
    freqHz=0.5,
    offset=100000,
    amplitude=20000,
    phase=0.017453292519943)
    annotation (Placement(transformation(extent={{-82,-64},{-62,-44}})));
  ClaRa.Basics.ControlVolumes.GasVolumes.VolumeGas_L2 flueGasCell1
    annotation (Placement(transformation(extent={{24,-70},{44,-50}})));
equation
  connect(sine.y, gasSink_pT.p) annotation (Line(
      points={{-61,-14},{-46,-14}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(flueGasCell.outlet, gasFlowSource_T.gas_a) annotation (Line(
      points={{44,-20},{60,-20}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(gasSink_pT.gas_a, SensorCO2.inlet)         annotation (Line(
      points={{-26,-20},{-10,-20}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(sine1.y, gasSink_pT1.p)
                                annotation (Line(
      points={{-61,-54},{-46,-54}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gasSink_pT1.gas_a, SensorCO1.inlet)        annotation (Line(
      points={{-26,-60},{-8,-60}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(SensorCO1.outlet, flueGasCell1.inlet) annotation (Line(
      points={{12,-60},{24,-60}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(SensorCO2.outlet, flueGasCell.inlet) annotation (Line(
      points={{10,-20},{24,-20}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(gasFlowSource_T1.gas_a, flueGasCell1.outlet) annotation (Line(
      points={{60,-60},{44,-60}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (
    Diagram(coordinateSystem(extent={{-100,-120},{100,100}},
          preserveAspectRatio=false),
            graphics={            Text(
          extent={{-96,94},{102,54}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=10,
          textString="______________________________________________________________________________________________
PURPOSE:

______________________________________________________________________________________________
"),                    Text(
          extent={{-112,102},{88,82}},
          lineColor={0,128,0},
          fontSize=31,
          textString="TESTED -- YYYY-MM-DD //XX"),Text(
          extent={{-96,50},{68,36}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=9,
          textString="______________________________________________________________________________________________________________
Remarks: 
______________________________________________________________________________________________________________
"),                   Text(
          extent={{-96,66},{104,48}},
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
end testGasCompositionSensor;
