within Exergy.XClaRa.Components.Mills.HardCoalMills;
model VerticalMill_L3
  "Vertical roller mill such as ball-and-race mill and roller-bowl mills"
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

  extends ClaRa.Basics.Icons.RollerBowlMill;
  extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L3");
  ClaRa.Basics.Interfaces.Connected2SimCenter connected2SimCenter(
    powerIn=0,
    powerOut=0,
    powerAux=P_mills) if contributeToCycleSummary;

  outer ClaRa.SimCenter simCenter;

////////////////// PARAMETERS /////////////////////////
//________Materials and Media_______
  parameter ClaRa.Basics.Media.Fuel.PartialFuel coal=simCenter.fuelModel1
    "Medium to be used"                                                                       annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));
  parameter TILMedia.GasTypes.BaseGas  gas= simCenter.flueGasModel
    "Medium to be used in tubes"                                                                annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));

//________Mill definition___________
  parameter
    Exergy.XClaRa.Components.Mills.HardCoalMills.Fundamentals.RollerBowlMillDefinition
    millKoeff=Fundamentals.STV1() "Coefficients of the mill"
    annotation (choicesAllMatching, Dialog(group="Mill Definition"));
  parameter Integer N_mills= 1 "Number of equal mills in parallel" annotation(Dialog(group="Mill Definition"));

//________Initialisation____________
  parameter SI.Temperature T_0=simCenter.T_amb_start
    "Initial temperature in mill"
    annotation (Dialog(group="Initialisation"));
  parameter SI.Mass M_c_0=1000 "Initial Mass of raw coal on the table"
    annotation (Dialog(group="Initialisation"));
  parameter SI.Mass M_pf_0=100 "Initial Mass of pulverized coal on the table"
    annotation (Dialog(group="Initialisation"));
  parameter SI.Mass M_cair_0=100 "Initial Mass of pulverized coal in the air"
    annotation (Dialog(group="Initialisation"));
  parameter SI.MassFraction xi_coal_0[coal.nc-1]=coal.defaultComposition
    "Initial coal composition"                                                                      annotation(Dialog(group="Initialisation"));
  parameter ClaRa.Basics.Choices.Init initChoice=ClaRa.Basics.Choices.Init.noInit
    "Initialisation option"                                                                               annotation(Dialog(group="Initialisation"));

//________Summary and Visualisation_
  parameter Boolean contributeToCycleSummary = simCenter.contributeToCycleSummary
    "True if component shall contribute to automatic efficiency calculation"                annotation(Dialog(tab="Summary and Visualisation"));

//________Expert Settings___________
  parameter Boolean applyGrindingDelay = false
    "True if grinding process introduces a dead time"                                           annotation(Dialog(enable=applyGrindingDelay, tab="Expert Settings"));
  parameter SI.Time Tau_grind = 120 "Grinding dead time" annotation(Dialog(enable=applyGrindingDelay, tab="Expert Settings"));

///////////////// VARAIABLE DECLARATION ///////////////
//_________Masses_______________
protected
  SI.Mass M_c(start=M_c_0) "Mass of ungrinded coal on the table //mass_rct";
  SI.Mass M_pf(start=M_pf_0) "Mass of pulverized coal on the table //mass_pct";
  SI.Mass M_cair(start=M_cair_0)
    "Mass of pulverized coal carried by primary air //mass_pca";

//________Mass Flows___________
  SI.MassFlowRate m_flow_coal_in
    "Mass flow rate of raw coal entering the grinding table //m_flow_rct";
  SI.MassFlowRate m_flow_coal_pf "Mass flow rate of pulverized coal//m_flow_pc";
  SI.MassFlowRate W_c_ "Mass flow of raw coal to the mill //m_flow_rc_in";
  SI.MassFlowRate m_flow_coal_out
    "Mass flow rate of wet pulverized coal from the grinding zone //m_flow_wc_out";
  SI.MassFlowRate m_flow_coal_ret
    "Mass flow rate of coal returnig to the table //m_flow_pc_ret";
  SI.MassFlowRate m_flow_air "Primary air inlet mass flow rate //m_flow_air_in";
  SI.MassFlowRate m_flow_H2O_evap
    "Mass flow rate of evaporated coal H2O //m_flow_evap";
  SI.MassFlowRate m_flow_H2O_evap_max
    "Maximum evaporation flow until air saturation //m_flow_air_evap_max";
  SI.MassFlowRate m_flow_coal_evap_max
    "Maximum evaporation flow until coal dry out";

  SI.MassFlowRate m_flow_air_out "Primary air outlet mass flow rate";
  SI.MassFlowRate m_flow_rcg "Mass flow of raw coal entering the grinding zone";
  SI.MassFlowRate m_flow_dc_out "Mass flow of dried coal leaving the mill";

//________Mass Fractions_______
  SI.MassFraction xi_coal_in[coal.nc-1]
    "Mositure content of incoming raw coal //xi_rc_in";
  SI.MassFraction xi_coal_mix[coal.nc-1]
    "Average composition of wet coal after grinding //xi_wc_out";

  SI.MassFraction xi_dc_out[coal.nc-1]
    "Coal composition of dried coal at outlet";

  SI.MassFraction xi_air_in[gas.nc-1] "Composition of incoming air";
  SI.MassFraction xi_air_out[ gas.nc-1] "Composition of outgoing air";

//________Pressures____________
  SI.Pressure Delta_p_pa(displayUnit="Pa") "Primary air difference pressure";

//________Temperatures_________
  SI.Temperature T_out(start=T_0) "Classifier Temperature (outlet temperature)";
  SI.Temperature T_coal_in "Coal inlet temperature";
  SI.Temperature T_air_in "Primary air inlet temperature";

//________Coal specifics_______
  SI.EnthalpyMassSpecific LHV_dry(start=(33907*coal.defaultComposition[
        1] + 142324*(coal.defaultComposition[2] - coal.defaultComposition[
        3]/8.) + 10465*coal.defaultComposition[5] - 2512*((1 - sum(
        coal.defaultComposition)) + 9*coal.defaultComposition[2]))*
        1000) "Lower heating value after drying inside mill";
  SI.EnthalpyMassSpecific Delta_h_evap "Heat of vaporization";
  SI.HeatCapacityMassSpecific cp_w
    "Specific heat capacity of liquid water in the raw coal";
  SI.HeatCapacityMassSpecific cp_dc_in
    "Specific heat capacity of ideally dried coal at inlet condition";

//_______Mechanics_____________
  Real P_grind "Power consumed for grinding in p.u.";
  SI.Frequency rpm=classifierSpeed "Rotational speed of clasiifier";

  //______Auxilliary Variables_

  ClaRa.Basics.Functions.ClaRaDelay.ExternalTable pointer_W_c= ClaRa.Basics.Functions.ClaRaDelay.ExternalTable()
    "Pointer for delay memory allocation";
  TILMedia.VLEFluidObjectFunctions.VLEFluidPointer H2O_props=
      TILMedia.VLEFluidObjectFunctions.VLEFluidPointer(
      TILMedia.VLEFluidTypes.TILMedia_SplineWater.concatVLEFluidName,
      0,
      TILMedia.VLEFluidTypes.TILMedia_SplineWater.mixingRatio_propertyCalculation[1:end - 1]/sum(TILMedia.VLEFluidTypes.TILMedia_SplineWater.mixingRatio_propertyCalculation),
      TILMedia.VLEFluidTypes.TILMedia_SplineWater.nc_propertyCalculation,
      TILMedia.VLEFluidTypes.TILMedia_SplineWater.nc,
      TILMedia.Internals.redirectModelicaFormatMessage())
    "Pointer to external medium memory for evaporation enthalpy and heat capacity";

public
  Modelica.Blocks.Interfaces.RealOutput Delta_p_mill(unit="Pa")
    "Pressure Difference between inlet and outlet connector"
                                                    annotation(Placement(transformation(extent={{100,50},
            {140,90}})));
  Modelica.Blocks.Interfaces.RealInput classifierSpeed "Speed of classifier"
    annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=270,
        origin={0,108})));
  TILMedia.Gas_pT     gasIn(
    p=inlet.flueGas.p,
    T=inStream(inlet.flueGas.T_outflow),
    xi=inStream(inlet.flueGas.xi_outflow),
    gasType=gas)
    annotation (Placement(transformation(extent={{-88,-70},{-68,-50}})));

  Exergy.XClaRa.Components.Mills.HardCoalMills.Fundamentals.SummaryMill
    summary(
    T_coal_in=T_coal_in,
    T_air_in=T_air_in,
    rpm_classifier=classifierSpeed*60,
    P_grind=P_grind,
    m_flow_air_out=-outlet.flueGas.m_flow,
    mass_coal=N_mills*(M_cair + M_pf + M_c),
    m_flow_coal_in=inlet.fuel.m_flow,
    m_flow_air_in=inlet.flueGas.m_flow,
    m_flow_tot_in=inlet.fuel.m_flow + inlet.flueGas.m_flow,
    m_flow_coal_out=-outlet.fuel.m_flow,
    m_flow_tot_out=-outlet.fuel.m_flow - outlet.flueGas.m_flow,
    T_out=T_out,
    xi_coal_h2o_in=1 - sum(xi_coal_in),
    xi_coal_h2o_out=1 - sum(xi_dc_out),
    xi_air_h2o_in=xi_air_in[8],
    xi_air_h2o_out=xi_air_out[8],
    xi_air_h2o_sat=gasOut.xi_s) annotation (Placement(transformation(
          extent={{-40,-100},{-20,-80}})));

  ClaRa.Basics.Interfaces.FuelFlueGas_inlet inlet(flueGas(Medium=gas),
      fuelType=coal) "Combined gas-and-coal(raw, wet) inlet"
    annotation (Placement(transformation(extent={{-110,-8},{-90,12}}),
        iconTransformation(extent={{-110,-10},{-90,10}})));
  ClaRa.Basics.Interfaces.FuelFlueGas_outlet outlet(flueGas(Medium=
          gas), fuelType=coal) "Combined gas-and-coal(pulverised, dry) outlet"
                                                    annotation (
      Placement(transformation(extent={{90,-10},{110,10}})));

  TILMedia.Gas_pT     gasOut(
    p=outlet.flueGas.p,
    T=T_out,
    gasType=gas,
    xi=xi_air_out)
    annotation (Placement(transformation(extent={{70,-68},{90,-48}})));
  Modelica.Blocks.Interfaces.RealOutput P_mills(unit="W")
    "Mill power of all parallel mills"              annotation(Placement(transformation(extent={{100,22},
            {140,62}})));

initial equation
  xi_coal_mix = xi_coal_0;

  if initChoice == ClaRa.Basics.Choices.Init.noInit then
  //do nothing
  elseif initChoice == ClaRa.Basics.Choices.Init.steadyDensity then
    der(M_c)=0;
    der(M_cair)=0;
    der(M_pf)=0;
  elseif initChoice == ClaRa.Basics.Choices.Init.steadyTemperature then
    der(T_out)=0;
  elseif initChoice == ClaRa.Basics.Choices.Init.steadyState then
    der(M_c)=0;
    der(M_cair)=0;
    der(M_pf)=0;
    der(T_out)=0;
  else
    assert(false,"Unknown initialisation option");
  end if;

equation
////////////////////////////////////////////
/// Additional Media Data                ///
  Delta_h_evap = TILMedia.VLEFluidObjectFunctions.dewSpecificEnthalpy_Txi(T_coal_in, {1}, H2O_props) - TILMedia.VLEFluidObjectFunctions.bubbleSpecificEnthalpy_Txi(T_coal_in, {1}, H2O_props);
  cp_w = TILMedia.VLEFluidObjectFunctions.specificIsobaricHeatCapacity_pTxi(inlet.fuel.p, T_coal_in, {1}, H2O_props);

////////////////////////////////////////////
/// Grinding Process                     ///
  m_flow_coal_pf = millKoeff.K_5*m_flow_air*M_pf;
  m_flow_coal_out = millKoeff.K_4*M_cair*(1-min(rpm,millKoeff.K_6)/millKoeff.K_6);
  m_flow_coal_ret = millKoeff.K_9*M_cair;
  m_flow_rcg = millKoeff.K_1*M_c;

  //this is beyond Nimcyks model: a dead time taking particle transport
  //  from the entrance to the grinding table
  if applyGrindingDelay then
    m_flow_coal_in = ClaRa.Basics.Functions.ClaRaDelay.getDelayValuesAtTime(
       pointer_W_c,
       time,
       W_c_,
       time-Tau_grind);
  else
    m_flow_coal_in = W_c_;
  end if;

////////////////////////////////////////////
/// Coal Mass Balance                    ///

//_______Mass balances for the grinding table and the transport area:
  der(M_c) = m_flow_coal_in + m_flow_coal_ret - m_flow_rcg;
  der(M_pf) = millKoeff.K_1*M_c - m_flow_coal_pf;
  der(M_cair) = m_flow_coal_pf - m_flow_coal_out - m_flow_coal_ret;

//_______Species balance in grinding area
  der(xi_coal_mix) = (m_flow_coal_in * xi_coal_in - m_flow_coal_out * xi_coal_mix - xi_coal_mix*(der(M_c) + der(M_pf) + der(M_cair))) /(M_c + M_pf + M_cair);

//_______Drying of coal, after grinding process
  m_flow_dc_out = m_flow_coal_out - m_flow_H2O_evap;
  xi_dc_out = (m_flow_coal_out*xi_coal_mix - m_flow_H2O_evap*{0,0,0,0,0,0})/m_flow_dc_out;

////////////////////////////////////////////
/// Coal Drying                          ///
  m_flow_H2O_evap_max = m_flow_air*(gasOut.xi_s-gasIn.xi[8]);//Maximum H2O evaporation mass flow until air is saturated
  m_flow_coal_evap_max = m_flow_coal_in*(1-sum(xi_coal_in));//Maximum H2O evaporation mass flow until coal is dry

  if noEvent(m_flow_coal_evap_max <= m_flow_H2O_evap_max) then //Amount of coal H2O evaporation (if < m_flow_H2O_evap_max then ideal drying)
    m_flow_H2O_evap = m_flow_coal_out*(1-sum(xi_coal_mix));
  else
    m_flow_H2O_evap = m_flow_H2O_evap_max;
  end if;

  if inlet.fuel.LHV_calculationType=="Verbandsformel" then
    LHV_dry =  (33907*outlet.fuel.xi_outflow[1] + 142324*(outlet.fuel.xi_outflow[2] - outlet.fuel.xi_outflow[3]/8.) + 10465*outlet.fuel.xi_outflow[5] - 2512*((1 - sum(outlet.fuel.xi_outflow)) + 9*outlet.fuel.xi_outflow[2]))*1000;
  elseif inlet.fuel.LHV_calculationType=="predefined" then
    //LHV_dry = (inStream(inlet.fuel.LHV_outflow) + Delta_h_evap*(1-sum(xi_coal_mix)))/sum(xi_coal_mix); //Reduced LHV by amount of evaporated water during drying process
    LHV_dry = (inStream(inlet.fuel.LHV_outflow) + Delta_h_evap*(1-sum(xi_coal_mix)))*(1-(1-sum(xi_dc_out)))/(1-(1-sum(xi_coal_mix))) - Delta_h_evap*(1-sum(xi_dc_out));//Effenberger, lower heating value after drying
    //LHV_dry = inStream(inlet.fuel.LHV_outflow)* (1+summary.xi_coal_h2o_in/(1-summary.xi_coal_h2o_in))/(1+summary.xi_coal_h2o_out/(1-summary.xi_coal_h2o_out));
  else
    LHV_dry = -1;
    assert(inlet.fuel.LHV_calculationType == "predefined" or inlet.fuel.LHV_calculationType == "Verbandsformel", "Please check your LHV calculation settings inside boundaries.");
  end if;

////////////////////////////////////////////
/// Gas Mass Balance                     ///
  m_flow_air_out = m_flow_air + m_flow_H2O_evap; //no air mass storage
  xi_air_out = (m_flow_air*xi_air_in + m_flow_H2O_evap*cat(1,zeros(gas.condensingIndex-1), {1}, zeros(gas.nc-gas.condensingIndex-1)))/m_flow_air_out;

////////////////////////////////////////////
/// Gas Moisturing                       ///

////////////////////////////////////////////
/// Global Energy Balance                ///
// the energy balance as in equation (7) of [1] but with the derivative of the coal mass coming from the der(U) term
  der(T_out)=1/millKoeff.K_11*((gasIn.T-273.15)*gasIn.cp*m_flow_air
                              + (1-sum(xi_coal_in)) * m_flow_coal_in*cp_w*(T_coal_in-273.15)
                              + sum(xi_coal_in) * m_flow_coal_in*coal.cp*(T_coal_in-273.15)
                              - m_flow_air_out *gasOut.cp * (gasOut.T-273.15)
                              - m_flow_H2O_evap * Delta_h_evap
                              - sum(xi_dc_out) * m_flow_dc_out * coal.cp * (T_out-273.15)
                              - (1-sum(xi_dc_out)) * m_flow_dc_out * cp_w * (T_out-273.15)
                              + millKoeff.K_10*P_grind*100 - (der(M_c)+der(M_pf)+der(M_cair))*coal.cp*(T_out-273.15));

////////////////////////////////////////////
/// Effort for Grinding                  ///
  P_grind = 0.01*(millKoeff.K_2*M_pf+millKoeff.K_3*M_c)+millKoeff.E_e;
  P_mills = millKoeff.P_nom*P_grind*N_mills;

////////////////////////////////////////////
/// Hydraulics                           ///
  Delta_p_pa = millKoeff.K_12*abs(m_flow_air)*m_flow_air/max(gasIn.d,0.0001); //From mill volume flow measurement using an orifice
  Delta_p_mill = millKoeff.K_7 * Delta_p_pa + millKoeff.K_8*M_cair*100;

////////////////////////////////////////////
/// Connector Couplings                  ///
  T_air_in = inStream(inlet.flueGas.T_outflow);
  m_flow_air = inlet.flueGas.m_flow/N_mills;
  xi_air_in = inStream(inlet.flueGas.xi_outflow);
  T_coal_in = inStream(inlet.fuel.T_outflow);
  W_c_ = inlet.fuel.m_flow/N_mills;
  xi_coal_in = inStream(inlet.fuel.xi_outflow);

  inlet.fuel.p = inlet.flueGas.p;
  inlet.fuel.T_outflow = T_out; //DUMMY value - backflow is not supported!
  inlet.fuel.xi_outflow = inStream(outlet.fuel.xi_outflow); //DUMMY value - backflow is not supported!
  inlet.fuel.LHV_outflow = LHV_dry; //DUMMY value - backflow is not supported!
  inlet.fuel.cp_outflow = inStream(outlet.fuel.cp_outflow);

  inlet.flueGas.p = outlet.flueGas.p + Delta_p_pa;
  inlet.flueGas.T_outflow = T_out; //DUMMY value - backflow is not supported!
  inlet.flueGas.xi_outflow = inStream(outlet.flueGas.xi_outflow);//DUMMY value - backflow is not supported!

  outlet.fuel.m_flow = -m_flow_dc_out*N_mills;
  outlet.fuel.T_outflow = T_out;
  outlet.fuel.LHV_outflow = LHV_dry;
  outlet.fuel.LHV_calculationType = inlet.fuel.LHV_calculationType;
  outlet.fuel.xi_outflow = xi_dc_out;
  //outlet.fuel.cp_outflow = (inStream(inlet.fuel.cp_outflow) + cp_w * (1-sum(xi_coal_mix)))/sum(xi_coal_mix);
  outlet.fuel.cp_outflow = inStream(inlet.fuel.cp_outflow) - ((1-sum(xi_coal_in)) - (1-sum(xi_dc_out)))*(cp_w - cp_dc_in);
  cp_dc_in = (inStream(inlet.fuel.cp_outflow) - (1-sum(xi_coal_in))*cp_w)/(sum(xi_coal_in));

  outlet.flueGas.T_outflow = T_out;
  outlet.flueGas.xi_outflow = xi_air_out;
  outlet.flueGas.m_flow = -m_flow_air_out*N_mills;

annotation (Diagram(graphics), Icon(graphics));
end VerticalMill_L3;
