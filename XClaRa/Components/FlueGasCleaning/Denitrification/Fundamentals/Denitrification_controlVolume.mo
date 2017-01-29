within Exergy.XClaRa.Components.FlueGasCleaning.Denitrification.Fundamentals;
model Denitrification_controlVolume
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
  inner parameter TILMedia.GasTypes.BaseGas               medium = simCenter.flueGasModel
    "Medium to be used in tubes" annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));

parameter Modelica.SIunits.Temperature T_NH3_O2_mixture = 273.15+200
    "Temperature of ammonia oxygen inlet";
parameter Real separationRate(min = 0, max = 1) = 0.995
    "Efficiency of NOx separation";

// standard formation enthalpy (T = 298.15K /p = 1.0 bar) for  components involved in deNOx catalysis
final parameter Modelica.SIunits.MolarInternalEnergy Delta_f_H_NO = 91.271e3
    "Standrad formation enthalpy nitric oxide";
final parameter Modelica.SIunits.MolarInternalEnergy Delta_f_H_NH3 = -45.940e3
    "Standrad formation  enthalpy ammonia";
final parameter Modelica.SIunits.MolarInternalEnergy Delta_f_H_H2O = -241.826e3
    "Standrad formation  enthalpy water";

//## V A R I A B L E   P A R T##################################################################################

// Quantaties for deNOx catalysis
Modelica.SIunits.MassFlowRate NH3_O2_m_flow "Mass flow of ammoinia oxygen mix";
Modelica.SIunits.MassFlowRate flueGasMixture_m_flow
    "Mass flow of flue gas mixture";
//Molar flow rates
Modelica.SIunits.MolarFlowRate n_flow_NOx_in
    "Molar flow rate of nitric oxides at inlet";
Modelica.SIunits.MolarFlowRate n_flow_O2_in
    "Molar flow rate of oxygen at inlet";
Modelica.SIunits.MolarFlowRate n_flow_NH3_req
    "Required molar flow rate of ammonia";
Modelica.SIunits.MolarFlowRate n_flow_O2_req
    "Required molar flow rate of oxygen";
// standard reaction enthalpy
Modelica.SIunits.MolarInternalEnergy Delta_R_H "Reaction enthalpy";
Modelica.SIunits.MassFraction NH3_O2_in_xi[medium.nc-1]
    "Inlet composition of ammonia oxygen mix";
Modelica.SIunits.MassFraction flueGasMixture_xi[medium.nc-1]
    "Inlet composition of flue gas";
Modelica.SIunits.HeatFlowRate Qdot "Heat flow to environment";

//____Connectors________________________________________________________________________________________________
  ClaRa.Basics.Interfaces.GasPortIn inlet(Medium=medium)  annotation (Placement(
        transformation(extent={{-110,-10},{-90,10}}),iconTransformation(extent={{-110,
            -10},{-90,10}})));
  ClaRa.Basics.Interfaces.GasPortOut outlet(Medium=medium)     annotation (Placement(
        transformation(extent={{90,-10},{110,10}}),  iconTransformation(extent={{90,-10},
            {110,10}})));

  //_____________________Media Objects_________________________________
  TILMedia.Gas_pT     flueGasInlet(p=inlet.p,
  T=inStream(inlet.T_outflow),
  xi=inStream(inlet.xi_outflow),
  gasType = medium)
    annotation (Placement(transformation(extent={{-70,50},{-50,70}})));

    TILMedia.Gas_pT     NH3_O2_in(
    p=inlet.p,
    T=T_NH3_O2_mixture,
    xi=NH3_O2_in_xi,
    gasType = medium)
    annotation (Placement(transformation(extent={{-10,50},{10,70}})));
  TILMedia.Gas_ph     flueGasMixture(xi=flueGasMixture_xi,
    gasType = medium)
    annotation (Placement(transformation(extent={{-40,10},{-20,30}})));
  TILMedia.Gas_pT     flueGasOutlet(
    p=inlet.p,
    T=outlet.T_outflow,
    xi=outlet.xi_outflow,
    gasType = medium)
    annotation (Placement(transformation(extent={{60,-70},{80,-50}})));

equation
  inlet.xi_outflow = zeros(medium.nc-1); // dummy values, flow reversal not allowed
  inlet.T_outflow = 273.15;

//   NH3_O2_in_xi = {0,0,0,0,0,if separationRate > 0 then NH3_O2_in.M_i[6]
//                                                                       *n_flow_O2_req/NH3_O2_m_flow else 0,0,0};
  for i in 1:(medium.nc-1) loop
      if i==6 then
        NH3_O2_in_xi[i] = if separationRate > 0 then
                            NH3_O2_in.M_i[6]*n_flow_O2_req/NH3_O2_m_flow
                          else 0;
        else
         NH3_O2_in_xi[i] = 0;
      end if;
  end for;

  // determination of required NH3 and O2 for deNOx catalysis
  n_flow_NOx_in = inlet.m_flow*flueGasInlet.xi[7]/flueGasInlet.M_i[7];
  n_flow_O2_in = inlet.m_flow*flueGasInlet.xi[6]/flueGasInlet.M_i[6];
  n_flow_NH3_req = separationRate*n_flow_NOx_in;
  n_flow_O2_req =
  if n_flow_O2_in >= 1/4.*n_flow_NH3_req then
    1e-12
  else
    1/4.*n_flow_NH3_req-n_flow_O2_in;

  NH3_O2_m_flow =NH3_O2_in.M_i[6]*n_flow_O2_req +NH3_O2_in.M_i[9]*n_flow_NH3_req;

  // balance equations for ideal mixture of inlet and NH3_O2_in
  flueGasMixture.p = inlet.p;
  flueGasMixture_m_flow + inlet.m_flow + NH3_O2_m_flow = 0; //mass balance
  flueGasMixture_m_flow*flueGasMixture.h + inlet.m_flow*flueGasInlet.h + NH3_O2_m_flow*NH3_O2_in.h = 0; // Enregy balance

  for i in 1:(medium.nc-1) loop //component mass balance
    flueGasMixture_m_flow*flueGasMixture.xi[i] + inlet.m_flow*flueGasInlet.xi[i] + NH3_O2_m_flow*NH3_O2_in.xi[i] = 0;
  end for;

//____________________________/ deNOx Catalysis \_______________________________________

// entire Mass balance of catalysis
outlet.m_flow - flueGasMixture_m_flow = 0;
//Component mass balance
for i in 1:(medium.nc-1) loop
    if i == 5 then
       flueGasOutlet.xi[5]*outlet.m_flow - flueGasMixture.xi[5]*flueGasMixture_m_flow + n_flow_NH3_req*
        flueGasMixture.M_i[5]                                                                                                = 0;
      else if i == 6 then
       flueGasOutlet.xi[6]*outlet.m_flow  -  flueGasMixture.xi[6]*flueGasMixture_m_flow - 1/4.*n_flow_NH3_req*
          flueGasMixture.M_i[6]                                                                                                    = 0;
      else if i == 7 then
       flueGasOutlet.xi[7]*outlet.m_flow - flueGasMixture.xi[7]*flueGasMixture_m_flow*(1-separationRate)  = 0;
      else if i == 8 then
       flueGasOutlet.xi[8]*outlet.m_flow  -  flueGasMixture.xi[8]*flueGasMixture_m_flow + 6/4.*n_flow_NH3_req*
              flueGasMixture.M_i[8]                                                                                                  = 0;
    else
    flueGasOutlet.xi[i]*outlet.m_flow  -  flueGasMixture.xi[i]*flueGasMixture_m_flow  = 0;
     end if;
    end if;
   end if;
  end if;
end for;

//inlet.m_flow + outlet.m_flow + NH3_O2_m_flow = 0;

//Energy balance
//reaction enthalpie of deNOx-Catalaysis
Delta_R_H = (4*Delta_f_H_NH3 + 4*Delta_f_H_NO  - 6*Delta_f_H_H2O);

//reactionHeat =(n_flow_N2_Cat_out*gasMixture_idealMixture.mm[5] + n_flow_H2O_Cat_out*gasMixture_idealMixture.mm[8])*Delta_R_H/(-(4.0*gasMixture_idealMixture.mm[5] + 6*gasMixture_idealMixture.mm[8]));

//energy balance

0 = -flueGasMixture_m_flow*(flueGasMixture.h
                        +flueGasMixture.xi[7]*Delta_f_H_NO/flueGasMixture.M_i[7]
                        +(1-sum(flueGasMixture.xi))*Delta_f_H_NH3/flueGasMixture.M_i[9]
                        +flueGasMixture.xi[8]*Delta_f_H_H2O/flueGasMixture.M_i[8])
      + outlet.m_flow*(flueGasOutlet.h
                       +flueGasOutlet.xi[7]*Delta_f_H_NO/flueGasOutlet.M_i[7]
                       +(1-sum(flueGasOutlet.xi))*Delta_f_H_NH3/flueGasOutlet.M_i[9]
                       +flueGasOutlet.xi[8]*Delta_f_H_H2O/flueGasOutlet.M_i[8]) - (Qdot);

flueGasOutlet.T = flueGasMixture.T;
//Qdot = 0;

 flueGasOutlet.p = outlet.p;

  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
            -100},{100,100}}),
                      graphics={
        Text(
          extent={{-39,4},{39,-4}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          origin={-1,-34},
          rotation=180,
          textString="4NH3 +4NO +O2"),
        Text(
          extent={{36,-30},{78,-38}},
          lineColor={0,0,0},
          textString="4N2 +6H2O"),
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
          points={{-28,30},{-30,34},{-26,34},{-28,30}},
          lineColor={0,0,0},
          lineThickness=0.5,
          smooth=Smooth.None,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-32,30},{-34,34},{-30,34},{-32,30}},
          lineColor={0,0,0},
          lineThickness=0.5,
          smooth=Smooth.None,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-10,-60},{24,-60},{60,-60}},
          color={0,0,0},
          smooth=Smooth.None,
          thickness=0.5),
        Line(
          points={{-30,-6},{-30,-48},{-28,-54},{-24,-58},{-18,-60},{-12,-60},{-10,
              -60}},
          color={0,0,0},
          thickness=0.5,
          smooth=Smooth.Bezier),
        Polygon(
          points={{-30,0},{-32,4},{-28,4},{-30,0}},
          lineColor={0,0,0},
          lineThickness=0.5,
          smooth=Smooth.None,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{0,-2},{-2,2},{2,2},{0,-2}},
          lineColor={0,0,0},
          lineThickness=0.5,
          smooth=Smooth.None,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          origin={58,-60},
          rotation=90),
        Line(
          points={{26,-34},{34,-34},{34,-34}},
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
          origin={34,-34},
          rotation=90),
        Rectangle(
          extent={{-28,-28},{76,-40}},
          lineColor={0,0,0},
          pattern=LinePattern.Solid,
          lineThickness=0.5)}), Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}), graphics));
end Denitrification_controlVolume;
