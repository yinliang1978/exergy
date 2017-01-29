within Exergy.XClaRa.Components.Control.PredictorModels_3508;
model TurbinesAndReheat_02_XRG
  "A predictor for the generator power including the HP and IP/LP turbines aswell as the energy storage in the reheater"
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

  extends ClaRa.Basics.Icons.ComplexityLevel(complexity="02");
  parameter Modelica.SIunits.Pressure p_nom= 240e5
    "Nominal pressure at inlet of HP turbine"                                       annotation(Dialog(group="Nominal values"));
  parameter Modelica.SIunits.MassFlowRate m_flow_HP= 419
    "Nominal mass flow rate at inlet of HP turbine" annotation(Dialog(group="Nominal values"));
  parameter Modelica.SIunits.MassFlowRate m_flow_IP= 370
    "Nominal mass flow rate at inlet of IP turbine" annotation(Dialog(group="Nominal values"));
  parameter Modelica.SIunits.Power P_G_nom= 804.89e6
    "Nominal generator power of the block"                                       annotation(Dialog(group="Nominal values"));

 parameter Real CL_Ip_HP[:,2]= {{0.6000e7,    0.0267e7},
                                      {0.8000e7,    0.0419e7},
                                      {1.0000e7,    0.0601e7},
                                      {1.2000e7,    0.0815e7},
                                      {1.4000e7,    0.1060e7},
                                      {1.6000e7,    0.1339e7},
                                      {1.8000e7,    0.1652e7},
                                      {2.0000e7,    0.2001e7},
                                      {2.2000e7,    0.2386e7},
                                      {2.4000e7,    0.2807e7}}
    "Characteristic line IP pressure over HP pressure " annotation(Dialog(group="Part Load Definition"));
parameter Real CL_Deltah_HP[:,2]={{0.6000e7,    0.0822e7},
                                      {0.8000e7,    0.0782e7},
                                      {1.0000e7,    0.0746e7},
                                      {1.2000e7,    0.0713e7},
                                      {1.4000e7,    0.0682e7},
                                      {1.6000e7,    0.0652e7},
                                      {1.8000e7,    0.0624e7},
                                      {2.0000e7,    0.0598e7},
                                      {2.2000e7,    0.0572e7},
                                      {2.4000e7,    0.0548e7}}
    "Characteristic line enthalpy difference over HP inlet pressure "   annotation(Dialog(group="Part Load Definition"));

 parameter Real CL_Deltah_IP[:,2]= {{0.2669e6,    1.0678e6},
                                      {0.4191e6,    1.1287e6},
                                      {0.6014e6,    1.1770e6},
                                      {0.8147e6,    1.2172e6},
                                      {1.0601e6,    1.2517e6},
                                      {1.3389e6,    1.2819e6},
                                      {1.6522e6,    1.3086e6},
                                      {2.0009e6,    1.3325e6},
                                      {2.3858e6,    1.3541e6},
                                      {2.8073e6,    1.3736e6}}
    "Characteristic line enthalpy difference over IP inlet pressure "  annotation(Dialog(group="Part Load Definition"));

parameter Modelica.SIunits.Time Tau_HP= (0.2+0.5)/2
    "Time Constant for Energy Storage in HP turbine"                                     annotation(Dialog(group="Time Response Definition"));
parameter Modelica.SIunits.Time Tau_IP= (10+25)/2
    "Time Constant for Energy Storage in IP/LP turbine"                                     annotation(Dialog(group="Time Response Definition"));

  parameter Real m_flow_start_=1 "Initial mass flow rate in p.u." annotation(Dialog(group="Initialisation"));

  Modelica.Blocks.Continuous.FirstOrder energyStorage_HP_turbine(        T=Tau_HP,
    initType=Modelica.Blocks.Types.Init.InitialOutput,
    y_start=m_flow_start_)
    annotation (Placement(transformation(extent={{-62,24},{-42,44}})));
  Modelica.Blocks.Continuous.FirstOrder energyStroage_RH_IPLP_turbine(
               T=Tau_IP,
    initType=Modelica.Blocks.Types.Init.InitialOutput,
    y_start=m_flow_start_)
    annotation (Placement(transformation(extent={{-62,-36},{-42,-16}})));
  Modelica.Blocks.Math.Add add
    annotation (Placement(transformation(extent={{36,0},{56,20}})));
  ClaRa.Basics.Interfaces.SteamSignal inlet annotation (Placement(
        transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-98,0}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-98,80})));
  Modelica.Blocks.Interfaces.RealOutput P_gen_ "Generator power in p.u."
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));
  Modelica.Blocks.Tables.CombiTable1D convert2SpecificHPturbinePower(columns={2},
      table=CL_Deltah_HP)
    annotation (Placement(transformation(extent={{-24,60},{-4,80}})));
  Modelica.Blocks.Math.Gain gain(k=p_nom)
    annotation (Placement(transformation(extent={{-58,60},{-38,80}})));
  Modelica.Blocks.Math.Gain gain1(k=m_flow_HP)
    annotation (Placement(transformation(extent={{-24,24},{-4,44}})));
  Modelica.Blocks.Math.Product product
    annotation (Placement(transformation(extent={{14,36},{34,56}})));
  Modelica.Blocks.Math.Gain gain2(k=1/P_G_nom)
    annotation (Placement(transformation(extent={{64,0},{84,20}})));
  Modelica.Blocks.Math.Product product1
    annotation (Placement(transformation(extent={{12,-42},{32,-22}})));
  Modelica.Blocks.Math.Gain gain3(k=m_flow_IP)
    annotation (Placement(transformation(extent={{-26,-36},{-6,-16}})));
  Modelica.Blocks.Tables.CombiTable1D convert2IntermediatePressure(columns={2},
      table=CL_Ip_HP)
    annotation (Placement(transformation(extent={{-42,-84},{-22,-64}})));
  Modelica.Blocks.Math.Gain gain4(
                                 k=p_nom)
    annotation (Placement(transformation(extent={{-70,-84},{-50,-64}})));
  Modelica.Blocks.Tables.CombiTable1D convert2SpecificLPturbinePower(columns={2},
      table=CL_Deltah_IP)
    annotation (Placement(transformation(extent={{-12,-84},{8,-64}})));
equation

assert(Tau_HP>0 and Tau_IP>0, "Time constants must be greater than zero!");
  connect(gain.y, convert2SpecificHPturbinePower.u[1])
                                                  annotation (Line(
      points={{-37,70},{-26,70}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain2.y, P_gen_) annotation (Line(
      points={{85,10},{90,10},{90,0},{100,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(inlet.p_, gain.u) annotation (Line(
      points={{-98,0},{-84,0},{-84,70},{-60,70}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(inlet.m_flow_, energyStroage_RH_IPLP_turbine.u) annotation (Line(
      points={{-98,0},{-78,0},{-78,-26},{-64,-26}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(inlet.m_flow_, energyStorage_HP_turbine.u) annotation (Line(
      points={{-98,0},{-78,0},{-78,34},{-64,34}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(energyStorage_HP_turbine.y, gain1.u) annotation (Line(
      points={{-41,34},{-26,34}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain1.y, product.u2) annotation (Line(
      points={{-3,34},{4,34},{4,40},{12,40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(convert2SpecificHPturbinePower.y[1], product.u1) annotation (Line(
      points={{-3,70},{4,70},{4,52},{12,52}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.y, gain2.u) annotation (Line(
      points={{57,10},{62,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(product.y, add.u1) annotation (Line(
      points={{35,46},{34,46},{34,16}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(inlet.p_, gain4.u) annotation (Line(
      points={{-98,0},{-84,0},{-84,-74},{-72,-74}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(gain4.y, convert2IntermediatePressure.u[1]) annotation (Line(
      points={{-49,-74},{-44,-74}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(convert2IntermediatePressure.y, convert2SpecificLPturbinePower.u)
    annotation (Line(
      points={{-21,-74},{-14,-74}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(convert2SpecificLPturbinePower.y[1], product1.u2) annotation (Line(
      points={{9,-74},{10,-74},{10,-38}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain3.y, product1.u1) annotation (Line(
      points={{-5,-26},{10,-26}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(energyStroage_RH_IPLP_turbine.y, gain3.u) annotation (Line(
      points={{-41,-26},{-28,-26}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(product1.y, add.u2) annotation (Line(
      points={{33,-32},{34,-32},{34,4}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Icon(graphics={
        Polygon(
          points={{-94,10},{-74,18},{-74,-18},{-94,-10},{-94,10}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-100,80},{-94,80},{-94,10}},
          color={0,0,0},
          smooth=Smooth.None),
        Rectangle(
          extent={{-46,100},{6,62}},
          lineColor={0,0,0},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-74,18},{-74,80},{-34,80},{-20,92},{-20,68},{-6,80},{20,80},{
              20,18}},
          color={0,0,0},
          smooth=Smooth.None),
        Polygon(
          points={{20,20},{40,40},{40,-40},{20,-20},{20,20}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{80,40},{100,60},{100,-60},{80,-40},{80,40}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{80,40},{60,60},{60,-60},{80,-40},{80,40}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Line(
          points={{40,40},{40,80},{80,80},{80,40}},
          color={0,0,0},
          smooth=Smooth.None),
        Rectangle(
          extent={{-74,2},{20,-2}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{40,2},{60,-2}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-100,2},{-94,-2}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid)}), Diagram(graphics),
    Documentation(info="<html>
<p>A model calculating the time response of a HP, IP and LP turbine with reheating</p>
<p>Be carefully: This model extends the VDI/VDE guideline for Unit Control of Thermal Power Stations in the following:</p>
<p><ul>
<li>The generator output is generated by retarding the turbine mass flow according to two parallel first order elements</li>
<li>this is true only if the mass-specific power outlet is constant during load change which is not the case in general</li>
<li>to cope with variable specific power output curves a generic characteristic line was introduced</li>
<li>the default characteristic line refers to an ideal expansion and does not take reduced mass flows due to tapping of the LP and IP turbine into account!</li>
</ul></p>
<p><br/>This model extends the VDI/VDE guieline for Unit Control of Thermal Power Stations in the following:</p>
<p><ul>
<li>two characteristic lines for the specific power output in dependance of the turbine inlet pressure are used opeing two sparate paths for HP turbine and IP/LP turbine</li>
<li>in comparisson to the model <a href=\"ClaRa.Components.Control.PredictorModels_3508.TurbinesAndReheat_01_XRG\">TurbinesAndReheat_01_XRG</a> this allows the definition of a non-constant turbine power ratio (P_HP/P_tot)</li>
</ul></p>
<p><br/>Contact: Friedrich Gottelt, XRG Simulation</p>
</html>", revisions="<html>
<p><ul>
<li><b>v0.1 </b>2011-07-12: Initial implementation. Friedrich Gottelt, XRG Simulation GmbH</li>
</ul></p>
</html>"));
end TurbinesAndReheat_02_XRG;
