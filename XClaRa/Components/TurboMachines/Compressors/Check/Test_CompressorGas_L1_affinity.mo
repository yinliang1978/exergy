within Exergy.XClaRa.Components.TurboMachines.Compressors.Check;
model Test_CompressorGas_L1_affinity
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
  inner ClaRa.SimCenter simCenter annotation (Placement(
        transformation(extent={{40,40},{60,60}})));
  CompressorGas_L1_affinity GasFanAdvanced(
    eta=0.85,
    J=0.1,
    rpm_fixed=3000,
    rpm_nom=3000,
    exp_hyd=0.5,
    Delta_p_eps=1,
    steadyStateTorque=true,
    useMechanicalPort=true,
    V_flow_max=1,
    Delta_p_max=1e5) annotation (Placement(transformation(
        extent={{-8,-8},{8,8}},
        rotation=0,
        origin={-22,-44})));
  BoundaryConditions.BoundaryGas_pTxi
                                gasSink_pT2(variable_p=true, p_const=120000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={14,-44})));
  Modelica.Blocks.Sources.Ramp step1(
    startTime=5,
    duration=5,
    height=1e5,
    offset=1.21e5)
                 annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={42,-50})));
  BoundaryConditions.BoundaryGas_pTxi
                                gasSink_pT3(variable_p=false, p_const=120000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-70,-44})));
  Modelica.Mechanics.Rotational.Sources.Speed speed1(
                                                    exact=true, useSupport=
        false)
    annotation (Placement(transformation(extent={{-68,2},{-48,22}})));
  Modelica.Mechanics.Rotational.Components.Inertia inertia1(J=10)
    annotation (Placement(transformation(extent={{-44,2},{-24,22}})));
  Modelica.Blocks.Sources.Ramp step2(
    duration=5,
    offset=2*Modelica.Constants.pi*3000/60,
    startTime=10,
    height=0.1*2*Modelica.Constants.pi*3000/60)
                 annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-86,12})));
equation
  connect(step1.y, gasSink_pT2.p)
                                annotation (Line(
      points={{31,-50},{24,-50}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(speed1.flange, inertia1.flange_a) annotation (Line(
      points={{-48,12},{-44,12}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(inertia1.flange_b, GasFanAdvanced.shaft) annotation (Line(
      points={{-24,12},{-22,12},{-22,-36}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(step2.y, speed1.w_ref) annotation (Line(
      points={{-75,12},{-70,12}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gasSink_pT3.gas_a, GasFanAdvanced.inlet) annotation (Line(
      points={{-60,-44},{-30,-44}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(GasFanAdvanced.outlet, gasSink_pT2.gas_a) annotation (Line(
      points={{-14,-44},{4,-44}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (
    Diagram(coordinateSystem(extent={{-100,-80},{60,60}}, preserveAspectRatio=false),
            graphics={          Text(
          extent={{-98,46},{-24,36}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=12,
          textString="________________________________________________________________
PURPOSE:
>> Tester for an affinity law based compressor model with different inputs


"),                             Text(
          extent={{-100,60},{-8,50}},
          lineColor={0,128,0},
          fontSize=34,
          textString="TESTED -- 2014-10-08 //LN")}),
    experiment(StopTime=15),
    __Dymola_experimentSetupOutput,
    Icon(coordinateSystem(extent={{-100,-100},{100,100}}, preserveAspectRatio=false)));
end Test_CompressorGas_L1_affinity;
