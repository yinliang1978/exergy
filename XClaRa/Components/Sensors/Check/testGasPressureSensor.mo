within Exergy.XClaRa.Components.Sensors.Check;
model testGasPressureSensor
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
  GasPressureSensor gasPressureSensor
    annotation (Placement(transformation(extent={{-10,-16},{10,4}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT(                  variable_p=true,
    variable_xi=false,
    xi_const={0,0,0,0,0.8,0.2,0,0,0})                                               annotation (Placement(transformation(extent={{-32,-46},{-12,-26}})));
  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T(m_flow_const=-10) annotation (Placement(transformation(extent={{80,-46},{60,-26}})));
  inner ClaRa.SimCenter simCenter
    annotation (Placement(transformation(extent={{78,76},{98,96}})));
  Modelica.Blocks.Sources.Sine sine(
    freqHz=0.5,
    offset=100000,
    amplitude=20000,
    phase=0.017453292519943)
    annotation (Placement(transformation(extent={{-78,-22},{-58,-2}})));
  ClaRa.Basics.ControlVolumes.GasVolumes.VolumeGas_L2 flueGasCell
    annotation (Placement(transformation(extent={{18,-46},{38,-26}})));
equation
  connect(gasPressureSensor.port, gasSink_pT.gas_a) annotation (Line(
      points={{0,-16},{0,-36},{-12,-36}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(sine.y, gasSink_pT.p) annotation (Line(
      points={{-57,-12},{-44,-12},{-44,-30},{-32,-30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gasSink_pT.gas_a, flueGasCell.inlet) annotation (Line(
      points={{-12,-36},{18,-36}},
      color={84,58,36},
      smooth=Smooth.None));
  connect(flueGasCell.outlet, gasFlowSource_T.gas_a) annotation (Line(
      points={{38,-36},{60,-36}},
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
end testGasPressureSensor;
