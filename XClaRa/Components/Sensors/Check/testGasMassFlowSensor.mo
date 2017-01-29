within Exergy.XClaRa.Components.Sensors.Check;
model testGasMassFlowSensor
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
  GasMassflowSensor gasMassFlowSensor
    annotation (Placement(transformation(extent={{-18,-20},{2,0}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT(variable_p=true,
    variable_xi=false,
    xi_const={0,0,0,0,0.8,0.2,0,0,0})                                               annotation (Placement(transformation(extent={{-52,-30},{-32,-10}})));
  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T(m_flow_const=-10, variable_xi=false) annotation (Placement(transformation(extent={{72,-30},{52,-10}})));
  inner ClaRa.SimCenter simCenter
    annotation (Placement(transformation(extent={{78,76},{98,96}})));
  Modelica.Blocks.Sources.Sine sine(
    freqHz=0.5,
    offset=100000,
    amplitude=20000,
    phase=0.017453292519943)
    annotation (Placement(transformation(extent={{-94,-10},{-74,10}})));
  ClaRa.Basics.ControlVolumes.GasVolumes.VolumeGas_L2 flueGasCell
    annotation (Placement(transformation(extent={{14,-30},{34,-10}})));
equation
  connect(sine.y, gasSink_pT.p) annotation (Line(
      points={{-73,0},{-66,0},{-66,-14},{-52,-14}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(flueGasCell.outlet, gasFlowSource_T.gas_a) annotation (Line(
      points={{34,-20},{52,-20}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(gasSink_pT.gas_a, gasMassFlowSensor.inlet) annotation (Line(
      points={{-32,-20},{-18,-20}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(gasMassFlowSensor.outlet, flueGasCell.inlet) annotation (Line(
      points={{2,-20},{14,-20}},
      color={84,58,36},
      smooth=Smooth.None));
  annotation (
    Diagram(graphics={            Text(
          extent={{-96,100},{102,60}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=10,
          textString="______________________________________________________________________________________________
PURPOSE:

______________________________________________________________________________________________
"),                    Text(
          extent={{-136,104},{64,84}},
          lineColor={0,128,0},
          fontSize=31,
          textString="TESTED -- YYYY-MM-DD //XX"),Text(
          extent={{-96,60},{68,46}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          textString="______________________________________________________________________________________________________________
Remarks: 
______________________________________________________________________________________________________________
",        fontSize=8),Text(
          extent={{-96,74},{104,56}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=10,
          textString="______________________________________________________________________________________________
Scenario:  

______________________________________________________________________________________________
")}),
    experiment(StopTime=100),
    __Dymola_experimentSetupOutput);
end testGasMassFlowSensor;
