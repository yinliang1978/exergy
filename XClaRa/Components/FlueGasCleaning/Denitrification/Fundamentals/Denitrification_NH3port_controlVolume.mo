within Exergy.XClaRa.Components.FlueGasCleaning.Denitrification.Fundamentals;
model Denitrification_NH3port_controlVolume
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
  extends ClaRa.Basics.Icons.Box;

  outer ClaRa.SimCenter simCenter;
//## P A R A M E T E R S #######################################################################################
//_____________defintion of medium used in cell__________________________________________________________
   parameter TILMedia.GasTypes.BaseGas               medium = simCenter.flueGasModel
    "Medium to be used" annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"),
    Placement(transformation(extent={{80,80},{100,100}})));

// standard formation enthalpy (T = 298.15K /p = 1.0 bar) for  components involved in deNOx catalysis
final parameter Modelica.SIunits.MolarInternalEnergy Delta_f_H_NO = 91.271e3
    "Standrad formation enthalpy nitric oxide";
final parameter Modelica.SIunits.MolarInternalEnergy Delta_f_H_NH3 = -45.940e3
    "Standrad formation  enthalpy ammonia";
final parameter Modelica.SIunits.MolarInternalEnergy Delta_f_H_H2O = -241.826e3
    "Standrad formation  enthalpy water";

//## V A R I A B L E   P A R T##################################################################################

Modelica.SIunits.MassFlowRate idealMixture_m_flow "Idealized mass flow";

//Actual Molar flowrates at input ports
Modelica.SIunits.MolarFlowRate n_flow_NH3_in
    "Molar flow rate of ammonia at inlet";
Modelica.SIunits.MolarFlowRate n_flow_O2_in
    "Molar flow rate of oxygen at inlet";
Modelica.SIunits.MolarFlowRate n_flow_NOx_in
    "Molar flow rate of nitric oxides at inlet";

// Molar flowrates related to the catalysis
Modelica.SIunits.MolarFlowRate n_flow_NH3_Cat_in
    "Catalysis inlet ammonia mass flow";
Modelica.SIunits.MolarFlowRate n_flow_O2_Cat_in
    "Catalysis inlet oxygen mass flow";
Modelica.SIunits.MolarFlowRate n_flow_NOx_Cat_in
    "Catalysis inlet nitric oxide mass flow";
Modelica.SIunits.MolarFlowRate n_flow_N2_Cat_out
    "Catalysis outlet nitrogen mass flow";
Modelica.SIunits.MolarFlowRate n_flow_H2O_Cat_out
    "Catalysis outlet water mass flow";

//excess of NH3, NOx, O2
Modelica.SIunits.MolarFlowRate n_flow_NH3_exc "Molar flow of excess ammonia";
Modelica.SIunits.MolarFlowRate n_flow_O2_Cat_exc "Molar flow of excess oxygen";
Modelica.SIunits.MolarFlowRate n_flow_NOx_Cat_exc
    "Molar flow of excess nitric oxide";

// standard reaction enthalpy
Modelica.SIunits.MolarInternalEnergy Delta_R_H "Standard reaction enthalpy";
Modelica.SIunits.HeatFlowRate Qdot
    "Heat flow to environment (equals reaction heat)";
Modelica.SIunits.HeatFlowRate reactionHeat "Reaction heat";

Real NOx_separationRate "Efficiency of NOx separation";

Integer Case;

Real xi_idealMixture[medium.nc- 1] "Gas mixture ideal";
Real xi_CatalysisOut[medium.nc- 1] "Gas mixture at catalysis outlet";
Real h_idealMixture "Specific enthalpy of gas mixture ideal";
Real h_CatalysisOut "Specific enthalpy of gas mixture at catalysis outlet";

//____Connectors________________________________________________________________________________________________
  ClaRa.Basics.Interfaces.GasPortIn flueGas_a(Medium=medium)  annotation (Placement(
        transformation(extent={{-110,-10},{-90,10}}),iconTransformation(extent={{-110,
            -10},{-90,10}})));
  ClaRa.Basics.Interfaces.GasPortIn NH3_inlet(Medium=medium)  annotation (Placement(
        transformation(extent={{-10,90},{10,110}}), iconTransformation(extent={{-10,90},
            {10,110}})));
  ClaRa.Basics.Interfaces.GasPortOut flueGas_b(Medium=medium)  annotation (Placement(
        transformation(extent={{90,-10},{110,10}}),  iconTransformation(extent={{90,-10},
            {110,10}})));

//_____________________Media Objects_________________________________
  TILMedia.Gas_pT     gasMixture_a(gasType = medium,p = flueGas_a.p, T =  actualStream(flueGas_a.T_outflow),  xi = actualStream(flueGas_a.xi_outflow))
    annotation (Placement(transformation(extent={{-70,50},{-50,70}})));
  TILMedia.Gas_pT     gasMixture_NH3in(gasType = medium, p = NH3_inlet.p, T = actualStream(NH3_inlet.T_outflow),xi = actualStream(NH3_inlet.xi_outflow))
    annotation (Placement(transformation(extent={{-10,50},{10,70}})));

  TILMedia.Gas_ph     gasMixture_idealMixture(gasType = medium, p = flueGas_a.p, h = h_idealMixture, xi = xi_idealMixture)
    annotation (Placement(transformation(extent={{-40,10},{-20,30}})));
  TILMedia.Gas_ph     gasMixture_CatalysisOut(
    gasType = medium, p = flueGas_b.p,  h = h_CatalysisOut, xi = xi_CatalysisOut)
    annotation (Placement(transformation(extent={{66,-70},{86,-50}})));

equation
  // Think about correct usage of stream operator
  flueGas_a.T_outflow = 0;//inStream(flueGas_b.T_outflow);
  flueGas_a.xi_outflow = inStream(flueGas_b.xi_outflow);
  NH3_inlet.T_outflow = 0;//inStream(flueGas_b.T_outflow);
  NH3_inlet.xi_outflow = inStream(flueGas_b.xi_outflow);

  flueGas_a.p - flueGas_b.p = 0;
  NH3_inlet.p -  flueGas_b.p = 0;
  //gasMixture_idealMixture.p =  flueGas_a.p;  // different pressure levels can be considered lateron

  //+++++++++++++++  first control volume for mixing of NH3 and flueGas  +++++++++++++++
  // total mass balance
  0 = flueGas_a.m_flow + NH3_inlet.m_flow + idealMixture_m_flow;
  // component mass balance
  for i in 1:(medium.nc-1) loop
   gasMixture_a.xi[i]*flueGas_a.m_flow + gasMixture_NH3in.xi[i]*NH3_inlet.m_flow + xi_idealMixture[i]*idealMixture_m_flow = 0;
  end for;
  // energy balance
  0 = flueGas_a.m_flow*gasMixture_a.h + NH3_inlet.m_flow*gasMixture_NH3in.h + idealMixture_m_flow*h_idealMixture;

  //+++++++++++++ determination of stoichometric coefficients +++++++++++++++//

//   n_flow_NH3_in = idealMixture_m_flow*(1.0-sum({gasMixture_idealMixture.xi[i] for i in 1:medium.nc-1}))*1/
//     gasMixture_idealMixture.M_i[9];
  n_flow_NH3_in = idealMixture_m_flow*gasMixture_idealMixture.xi[9]*1/
    gasMixture_idealMixture.M_i[9];
  n_flow_O2_in = idealMixture_m_flow*gasMixture_idealMixture.xi[6]*1/
    gasMixture_idealMixture.M_i[6];
  n_flow_NOx_in = idealMixture_m_flow*gasMixture_idealMixture.xi[7]*1/
    gasMixture_idealMixture.M_i[7];

  // cases
  if (n_flow_NH3_in/n_flow_NOx_in >= 1.0) then // excess of NH3
    if
      (n_flow_O2_in/n_flow_NOx_in >= 1/4.0) then //excess of O2
      Case = 1;
      n_flow_NH3_Cat_in = n_flow_NOx_in;
      n_flow_O2_Cat_in = 1/4.0*n_flow_NOx_in;
      n_flow_NOx_Cat_in = n_flow_NOx_in;
    else // excess of NOx
      Case = 2;
      n_flow_O2_Cat_in = n_flow_O2_in;
      n_flow_NOx_Cat_in = 4.0*n_flow_O2_in;
      n_flow_NH3_Cat_in = 4.0*n_flow_O2_in;
    end if;
  else if (n_flow_O2_in/n_flow_NH3_in >= 1/4.0) then //deficit of NH3 and excess of O2
         Case = 3;
       n_flow_O2_Cat_in = 1/4.0*n_flow_NH3_in;
       n_flow_NOx_Cat_in = n_flow_NH3_in;
       n_flow_NH3_Cat_in = n_flow_NH3_in;
       else // deficit of O2
         Case = 4;
       n_flow_O2_Cat_in = n_flow_O2_in;
       n_flow_NH3_Cat_in = 4.0*n_flow_O2_in;
       n_flow_NOx_Cat_in = 4.0*n_flow_O2_in;
      end if;
  end if;

// excess of NH3, O2, NOx due to Stoichometry > < 1
n_flow_NH3_exc = n_flow_NH3_in-n_flow_NH3_Cat_in;
n_flow_O2_Cat_exc = n_flow_O2_in - n_flow_O2_Cat_in;
n_flow_NOx_Cat_exc = n_flow_NOx_in - n_flow_NOx_Cat_in;

// mass balance of catalysis
n_flow_N2_Cat_out = 4./6.*n_flow_H2O_Cat_out;
n_flow_NH3_Cat_in*gasMixture_idealMixture.M_i[9]+ n_flow_NOx_Cat_in*
    gasMixture_idealMixture.M_i[7] + n_flow_O2_Cat_in*gasMixture_idealMixture.M_i[6]  = n_flow_H2O_Cat_out* gasMixture_idealMixture.M_i[8]  + n_flow_N2_Cat_out* gasMixture_idealMixture.M_i[5];

// composition of flueGas after catalysis
// entire mass balance
flueGas_b.m_flow = n_flow_NH3_exc*gasMixture_idealMixture.M_i[9]+ n_flow_O2_Cat_exc*
    gasMixture_idealMixture.M_i[6]                                                                                + n_flow_NOx_Cat_exc*
    gasMixture_idealMixture.M_i[7]
                  + n_flow_N2_Cat_out*gasMixture_idealMixture.M_i[5]+ n_flow_H2O_Cat_out*
    gasMixture_idealMixture.M_i[8]
                  + idealMixture_m_flow*( sum({gasMixture_idealMixture.xi[i] for i in 1:5}) + gasMixture_idealMixture.xi[8]);

 // component mass balance
 for i in 1:(medium.nc-1) loop
    if i == 5 then
      xi_CatalysisOut[5]*flueGas_b.m_flow - gasMixture_idealMixture.xi[5]*idealMixture_m_flow -
        gasMixture_idealMixture.M_i[5]                                                                                                   *n_flow_N2_Cat_out = 0;
      else if i == 6 then
       xi_CatalysisOut[6]*flueGas_b.m_flow - gasMixture_idealMixture.xi[6]*idealMixture_m_flow +
          gasMixture_idealMixture.M_i[6]                                                                                                 *n_flow_O2_Cat_in = 0;
      else if i == 7 then
        xi_CatalysisOut[7]*flueGas_b.m_flow - gasMixture_idealMixture.xi[7]*idealMixture_m_flow +
            gasMixture_idealMixture.M_i[7]                                                                                               *n_flow_NOx_Cat_in = 0;
      else if i == 8 then
        xi_CatalysisOut[8]*flueGas_b.m_flow - gasMixture_idealMixture.xi[8]*idealMixture_m_flow -
              gasMixture_idealMixture.M_i[8]                                                                                             *n_flow_H2O_Cat_out = 0;
    else
     xi_CatalysisOut[i]*flueGas_b.m_flow - gasMixture_idealMixture.xi[i]*idealMixture_m_flow = 0;
     end if;
    end if;
   end if;
  end if;
 end for;

//reaction enthalpie of deNOx-Catalaysis
Delta_R_H = 4*Delta_f_H_NH3 + 4*Delta_f_H_NO  - 6*Delta_f_H_H2O;

reactionHeat = (n_flow_N2_Cat_out*gasMixture_idealMixture.M_i[5]
                +n_flow_H2O_Cat_out*gasMixture_idealMixture.M_i[8]) * Delta_R_H/(0.1*(4.0*gasMixture_idealMixture.M_i[5] + 6*gasMixture_idealMixture.M_i[8]));
 //energy balance with reactionheat
// 0 = -idealMixture_m_flow*gasMixture_idealMixture.h + flueGas_b.m_flow *gasMixture_CatalysisOut.h  +  reactionHeat + Qdot;
0 = -idealMixture_m_flow*gasMixture_idealMixture.h + flueGas_b.m_flow *h_CatalysisOut  - (Qdot-reactionHeat);
  Qdot = reactionHeat;
  NOx_separationRate = 1 + flueGas_b.m_flow*gasMixture_CatalysisOut.xi[7]/(flueGas_a.m_flow*gasMixture_a.xi[7]);

  //gasMixture_CatalysisOut.T = gasMixture_idealMixture.T;
 // gasMixture_CatalysisOut.p = flueGas_b.p;

flueGas_b.T_outflow =  gasMixture_CatalysisOut.T;
  flueGas_b.xi_outflow
                     = gasMixture_CatalysisOut.xi;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                   graphics), Diagram(coordinateSystem(extent={{-100,-100},{100,
            100}}, preserveAspectRatio=false),
                                      graphics={
        Rectangle(
          extent={{-76,72},{16,4}},
          lineColor={0,0,0},
          pattern=LinePattern.Dot,
          lineThickness=0.5),
        Line(
          points={{-50,60},{-46,60},{-40,60},{-32,60},{-32,52},{-32,46},{-32,34},
              {-32,32}},
          color={0,0,0},
          thickness=0.5,
          smooth=Smooth.Bezier),
        Line(
          points={{-10,60},{-14,60},{-20,60},{-28,60},{-28,52},{-28,46},{-28,34},
              {-28,32}},
          color={0,0,0},
          thickness=0.5,
          smooth=Smooth.Bezier),
        Polygon(
          points={{-32,30},{-34,34},{-30,34},{-32,30}},
          lineColor={0,0,0},
          lineThickness=0.5,
          smooth=Smooth.None,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-28,30},{-30,34},{-26,34},{-28,30}},
          lineColor={0,0,0},
          lineThickness=0.5,
          smooth=Smooth.None,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-30,0},{-32,4},{-28,4},{-30,0}},
          lineColor={0,0,0},
          lineThickness=0.5,
          smooth=Smooth.None,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-72,8},{-54,6}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="Control Volume"),
        Rectangle(
          extent={{-4,-24},{50,-36}},
          lineColor={0,0,0},
          pattern=LinePattern.Dot,
          lineThickness=0.5),
        Text(
          extent={{-18,-28},{60,-32}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="4NH3 +4NO +O2 -> 4N2 +6H2O"),
        Line(
          points={{-30,-4},{-30,-46},{-28,-52},{-24,-56},{-18,-58},{-12,-58},{-10,
              -58}},
          color={0,0,0},
          thickness=0.5,
          smooth=Smooth.Bezier),
        Polygon(
          points={{0,-2},{-2,2},{2,2},{0,-2}},
          lineColor={0,0,0},
          lineThickness=0.5,
          smooth=Smooth.None,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          origin={62,-58},
          rotation=90),
        Polygon(
          points={{0,-2},{-2,2},{2,2},{0,-2}},
          lineColor={0,0,0},
          lineThickness=0.5,
          smooth=Smooth.None,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          origin={-8,-30},
          rotation=90),
        Line(
          points={{-6,-30},{-22,-30},{-28,-28},{-30,-24},{-30,-18},{-30,-12}},
          color={0,0,0},
          thickness=0.5,
          smooth=Smooth.Bezier),
        Line(
          points={{9,14},{5,14},{-1,14},{-9,14},{-9,6},{-9,0},{-9,-12},{-9,-14}},
          color={0,0,0},
          thickness=0.5,
          smooth=Smooth.Bezier,
          origin={47,-46},
          rotation=90),
        Line(
          points={{-10,-58},{24,-58},{60,-58}},
          color={0,0,0},
          smooth=Smooth.None,
          thickness=0.5),
        Polygon(
          points={{0,-2},{-2,2},{2,2},{0,-2}},
          lineColor={0,0,0},
          lineThickness=0.5,
          smooth=Smooth.None,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          origin={62,-54},
          rotation=90),
        Rectangle(
          extent={{-36,-6},{56,-74}},
          lineColor={0,0,0},
          pattern=LinePattern.Dot,
          lineThickness=0.5),
        Text(
          extent={{-36,-70},{-18,-72}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="Control Volume")}));
end Denitrification_NH3port_controlVolume;
