within Exergy.XClaRa.Components.Control.PredictorModels_3508.Check;
model testReheater
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
  Modelica.Blocks.Sources.Ramp p(
    height=-0.5,
    duration=200,
    startTime=2000,
    offset=1)
    annotation (Placement(transformation(extent={{-82,24},{-62,44}})));
  Exergy.XClaRa.Components.Control.PredictorModels_3508.TurbinesAndReheat_01_XRG
    turbinesAndReheat(m_flow_start_=1)
    annotation (Placement(transformation(extent={{18,4},{38,24}})));
  Modelica.Blocks.Sources.Ramp p1(
    height=-0.5,
    duration=200,
    startTime=2000,
    offset=1)
    annotation (Placement(transformation(extent={{-82,-4},{-62,16}})));
  Modelica.Blocks.Sources.Ramp m(
    height=-0.5,
    duration=200,
    startTime=2000,
    offset=1)
    annotation (Placement(transformation(extent={{-56,-26},{-36,-6}})));
protected
  ClaRa.Basics.Interfaces.SteamSignal steamSignal annotation (
      Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-6,16}), iconTransformation(extent={{-28,-6},{-8,14}})));
public
  inner ClaRa.SimCenter simCenter annotation (Placement(
        transformation(extent={{80,80},{100,100}})));
  TurbinesAndReheat_02_XRG turbinesAndReheat_02_XRG(m_flow_start_=1)
    annotation (Placement(transformation(extent={{20,-40},{40,-20}})));
  TurbinesAndReheat_00_XRG turbinesAndReheat_00_XRG(eta_Gen=1)
    annotation (Placement(transformation(extent={{20,-80},{40,-60}})));
equation
  connect(p.y, steamSignal.p_) annotation (Line(
      points={{-61,34},{-6,34},{-6,16}},
      color={0,0,127},
      smooth=Smooth.None), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}}));
  connect(p1.y, steamSignal.h_) annotation (Line(
      points={{-61,6},{-16,6},{-16,16},{-6,16}},
      color={0,0,127},
      smooth=Smooth.None), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}}));
  connect(m.y, steamSignal.m_flow_) annotation (Line(
      points={{-35,-16},{-6,-16},{-6,16}},
      color={0,0,127},
      smooth=Smooth.None), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}}));
  connect(steamSignal, turbinesAndReheat.inlet) annotation (Line(
      points={{-6,16},{8,16},{8,22},{18.2,22}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(steamSignal, turbinesAndReheat_02_XRG.inlet) annotation (Line(
      points={{-6,16},{8,16},{8,-22},{20.2,-22}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(turbinesAndReheat_00_XRG.inlet, steamSignal) annotation (Line(
      points={{20.2,-62},{8,-62},{8,16},{-6,16}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}}));
  annotation (Diagram(graphics={  Text(
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
")}));
end testReheater;
