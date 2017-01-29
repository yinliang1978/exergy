within Exergy.XClaRa.SubSystems.Boiler;
model SteamGenerator_L3
  "A steam generation and reaheater model using lumped balance equations for mass and energy and two spray attemperators"
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

  extends Exergy.XClaRa.SubSystems.Boiler.CoalSupplyBoiler_base;
  ClaRa.Basics.Interfaces.Connected2SimCenter connected2SimCenter(
    powerIn=Q_flow_HP + Q_flow_IP,
    powerOut=0,
    powerAux=Q_flow_F_nom*QF_setl_ - Q_flow_HP - Q_flow_IP) if                                                                                                     contributeToCycleSummary;
  extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L3");

   parameter Modelica.SIunits.Pressure p_LS_nom= 300e5
    "Nominal life steam pressure"                                                  annotation(Dialog(group="Nominal values"));
   parameter Modelica.SIunits.Pressure p_RH_nom= 40e5 "Nominal reheat pressure"
                                                                               annotation(Dialog(group="Nominal values"));
   parameter Modelica.SIunits.SpecificEnthalpy h_LS_nom= 3000e3
    "Nominal life steam specific enthlapy"                                                  annotation(Dialog(group="Nominal values"));
   parameter Modelica.SIunits.SpecificEnthalpy h_RH_nom= 3500e3
    "Nominal reheat specific enthlapy"                                                  annotation(Dialog(group="Nominal values"));
  parameter Modelica.SIunits.Pressure Delta_p_nomHP = 40e5
    "Nominal main pressure loss"                                                         annotation(Dialog(group="Nominal values"));
  parameter Modelica.SIunits.Pressure Delta_p_nomIP = 4e5
    "Nominal reheat pressure loss"                                                        annotation(Dialog(group="Nominal values"));
  parameter Modelica.SIunits.MassFlowRate m_flow_nomLS = 419
    "Nominal life steam flow rate"                                                      annotation(Dialog(group="Nominal values"));
  parameter Modelica.SIunits.HeatFlowRate Q_flow_F_nom = 1340e6
    "Nominal firing power"                                                      annotation(Dialog(group="Nominal values"));
protected
  parameter Modelica.SIunits.Density rho_nom_HP= TILMedia.VLEFluidFunctions.density_phxi(medium, p_LS_nom, h_LS_nom)
    "Nominal density";
  parameter Modelica.SIunits.Density rho_nom_IP= TILMedia.VLEFluidFunctions.density_phxi(medium, p_RH_nom, h_RH_nom)
    "Nominal density";
public
  parameter Real CL_Delta_pHP_mLS_[:,:]=[0,0;0.1, 0.01; 0.2, 0.04; 0.3, 0.09; 0.4, 0.16; 0.5, 0.25; 0.6, 0.36; 0.7, 0.49; 0.8, 0.64; 0.9, 0.81; 1, 1]
    "Characteristic line of pressure drop as function of mass flow rate" annotation(Dialog(group="Part Load Definition"));

  parameter Real CL_Delta_pIP_mLS_[:,:]=[0,0;0.1, 0.01; 0.2, 0.04; 0.3, 0.09; 0.4, 0.16; 0.5, 0.25; 0.6, 0.36; 0.7, 0.49; 0.8, 0.64; 0.9, 0.81; 1, 1]
    "Characteristic line of reheat pressure drop as function of mass flow rate"
                                                                                annotation(Dialog(group="Part Load Definition"));
  parameter Real CL_yF_QF_[:,:] = [0,0.9;1,0.9]
    "Characteristic line of relative heat release in life steam as function of rel. firing power"
                                                                                    annotation(Dialog(group="Part Load Definition"));
  parameter Real CL_etaF_QF_[:,:] = [0,0.93;1,0.94]
    "Characteristic line of furnace efficiency as function of rel. firing power"
                                                                                 annotation(Dialog(group="Part Load Definition"));

  parameter Modelica.SIunits.Time Tau_dead = 120
    "Equivalent dead time of steam generation"                                        annotation(Dialog(group="Time Response Definition"));
  parameter Modelica.SIunits.Time Tau_bal = 200
    "Balancing time of steam generation"                                        annotation(Dialog(group="Time Response Definition"));
  parameter Boolean useHomotopy=simCenter.useHomotopy
    "True, if homotopy method is used during initialisation"  annotation(Dialog(tab="Initialisation", group="General"));
  parameter Modelica.SIunits.Pressure p_LS_start=300e5
    "Initial value of life steam pressure"                                                       annotation(Dialog(tab="Initialisation", group="High pressure part"));
  parameter Modelica.SIunits.SpecificEnthalpy h_LS_start=3000e3
    "Initial value of life steam specific enthalpy"                                                             annotation(Dialog(tab="Initialisation", group="High pressure part"));

//  parameter Real dotm_init_=1 "Initial mass flow rate in p.u." annotation(Dialog(group="Initialisation"));
  parameter ClaRa.Basics.Choices.Init initHP=ClaRa.Basics.Choices.Init.noInit
    "Type of initialisation of steam generation"  annotation(Dialog(tab="Initialisation", group="High pressure part"));

  parameter Modelica.SIunits.Pressure p_RH_start=40e5
    "Initial value of hot reheat pressure"                                                    annotation(Dialog(tab="Initialisation", group="Reheat part"));
  parameter Modelica.SIunits.SpecificEnthalpy h_RH_start=3500e3
    "Initial value of hot reheat specifc enthalpy"                                            annotation(Dialog(tab="Initialisation", group="Reheat part"));
  parameter ClaRa.Basics.Choices.Init initIP=ClaRa.Basics.Choices.Init.noInit
    "Type of initialisation of reheater" annotation(Dialog(tab="Initialisation", group="Reheat part"));

  parameter Modelica.SIunits.Volume volume_tot_HP=1000
    "Total volume of the live steam generator"                                                     annotation(Dialog(group="Geometry"));
  parameter Modelica.SIunits.Volume volume_tot_IP=1000
    "Total volume of the reheater"                                                    annotation(Dialog(group="Geometry"));

//___________Summary and Visualisation_____________________________________________//
  parameter Boolean showExpertSummary=false
    "|Summary and Visualisation||True, if expert summary should be applied";
  parameter Boolean showData=true
    "|Summary and Visualisation||True, if a data port containing p,T,h,s,m_flow shall be shown, else false";
  parameter Boolean contributeToCycleSummary = simCenter.contributeToCycleSummary
    "True if component shall contribute to automatic efficiency calculation"                  annotation(Dialog(tab="Summary and Visualisation"));

  outer ClaRa.SimCenter simCenter;

//___________Variables____________________________________________________________//
  Modelica.SIunits.SpecificEnthalpy h_IP(start=h_RH_start)
    "Specific enthalpy before IP injector";
  Modelica.SIunits.Pressure p_IP(start=p_RH_start)
    "Pressure at hot reheat outlet";
  Modelica.SIunits.Mass mass_IP "Mass in the reheater";
  Real drhodt_IP "Time derivative of reheater mean density";
  Modelica.SIunits.SpecificEnthalpy h_HP(start=h_LS_start)
    "Specific enthalpy before HP injector";
  Modelica.SIunits.Pressure p_HP(start=p_LS_start) "Live steam pressure";
  Modelica.SIunits.Mass mass_HP "Mass in the HP steam generator";
  Real drhodt_HP "Time dericative of the HP mean density";
  Modelica.SIunits.HeatFlowRate Q_flow_HP
    "Heat flow rate for HP steam generation";
  Modelica.SIunits.HeatFlowRate Q_flow_IP "Heat flow rate of the reheater";
protected
  Modelica.SIunits.SpecificEnthalpy h_inHP
    "Actual spec. enthalpy of the feedwater";
  Modelica.SIunits.SpecificEnthalpy h_outHP
    "Actual spec. enthalpy of the livesteam";
  Modelica.SIunits.SpecificEnthalpy h_sprayHP
    "Actual spec. enthalpy of the HP injection";
  Modelica.SIunits.SpecificEnthalpy h_inIP
    "Actual spec. enthalpy of the cold reheat";
  Modelica.SIunits.SpecificEnthalpy h_outIP
    "Actual spec. enthalpy of the hot reheat";
  Modelica.SIunits.SpecificEnthalpy h_sprayIP
    "Actual spec. enthalpy of the IP injection";
  Modelica.SIunits.MassFlowRate   m_flow_heatedHP
    "heated HP mass flow rate i.e. for energy and mass balance";
  Modelica.SIunits.MassFlowRate   m_flow_heatedIP
    "heated IP mass flow rate i.e. for energy and mass balance";

public
  Modelica.Blocks.Continuous.TransferFunction heatRelease(a={Tau_bal*Tau_dead,(Tau_bal + Tau_dead),1},
      initType=Modelica.Blocks.Types.Init.NoInit)
    "comprehends the coal supply, the heat release and the steam generation"
    annotation (Placement(transformation(extent={{-66,-56},{-46,-36}})));
  ClaRa.Components.Utilities.Blocks.ParameterizableTable1D convert2PressureDrop_HP(columns={2},
      table=CL_Delta_pHP_mLS_)
    annotation (Placement(transformation(extent={{-6,108},{14,128}})));
  TILMedia.VLEFluid_ph liveSteam(vleFluidType =    medium,       p=p_HP,
    h=(h_HP*(-m_flow_heatedHP) + HPInjection.m_flow*h_sprayHP)/(-livesteam.m_flow))
    annotation (Placement(transformation(extent={{-10,144},{10,164}})));
  TILMedia.VLEFluid_ph reheatedSteam(vleFluidType =    medium,   p=p_IP,
    h=(h_IP*(-m_flow_heatedIP) + IPInjection.m_flow*h_sprayIP)/(-reheat_out.m_flow))
    annotation (Placement(transformation(extent={{50,144},{70,164}})));

  Modelica.Blocks.Interfaces.RealOutput h_evap
    "evaporator outlet specific enthalpy" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={104,16}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={110,0})));
  ClaRa.Components.Utilities.Blocks.ParameterizableTable1D convert2PressureDrop_IP(columns={2},
      table=CL_Delta_pIP_mLS_)
    annotation (Placement(transformation(extent={{50,108},{70,128}})));
  ClaRa.Components.Utilities.Blocks.ParameterizableTable1D convert2HPFiring(columns={2}, table=
        CL_yF_QF_)
    annotation (Placement(transformation(extent={{-12,-56},{8,-36}})));
   ClaRa.Basics.Interfaces.FluidPortIn IPInjection(Medium=medium)
    "reheat spray injection"
     annotation (Placement(transformation(extent={{-110,124},{-90,144}}),
        iconTransformation(extent={{70,150},{90,170}})));
   ClaRa.Basics.Interfaces.FluidPortIn HPInjection(Medium=medium)
    "High pressure spray cooler"
     annotation (Placement(transformation(extent={{-110,90},{-90,110}}),
        iconTransformation(extent={{-28,150},{-8,170}})));
  ClaRa.Components.Utilities.Blocks.ParameterizableTable1D calculateEfficiency(columns={2}, table=
        CL_etaF_QF_)
    annotation (Placement(transformation(extent={{-12,-84},{8,-64}})));
  ClaRa.Basics.Interfaces.EyeOut eye_LS if showData annotation (Placement(
        transformation(
        extent={{-6,-6},{6,6}},
        rotation=90,
        origin={-40,188}), iconTransformation(
        extent={{-6,-6},{6,6}},
        rotation=90,
        origin={-40,226})));
protected
  ClaRa.Basics.Interfaces.EyeIn eye_intLS
    annotation (Placement(transformation(extent={{-43,159},{-41,161}})));
public
  ClaRa.Basics.Interfaces.EyeOut eye_RH if showData annotation (Placement(
        transformation(
        extent={{-6,-6},{6,6}},
        rotation=90,
        origin={88,188}), iconTransformation(
        extent={{-6,-6},{6,6}},
        rotation=90,
        origin={100,226})));
protected
  ClaRa.Basics.Interfaces.EyeIn eye_intRH
    annotation (Placement(transformation(extent={{85,171},{87,173}})));

protected
TILMedia.VLEFluid_ph fluidIn1(
    each vleFluidType=medium,
    h=actualStream(feedwater.h_outflow),
    p=feedwater.p)                                                           annotation (Placement(transformation(extent={{-50,-96},
            {-30,-76}},                                                                                                 rotation=0)));
TILMedia.VLEFluid_ph fluidOut1(
    each vleFluidType=medium,
    h=actualStream(livesteam.h_outflow),
    p=livesteam.p)                                                         annotation (Placement(transformation(extent={{-14,168},
            {14,188}},                                                                                                  rotation=0)));

TILMedia.VLEFluid_ph fluidIn2(
    each vleFluidType=medium,
    h=actualStream(reheat_in.h_outflow),
    p=reheat_in.p)                                                           annotation (Placement(transformation(extent={{74,-96},
            {94,-76}},                                                                                                  rotation=0)));
TILMedia.VLEFluid_ph fluidOut2(
    each vleFluidType=medium,
    h=actualStream(reheat_out.h_outflow),
    p=reheat_out.p)                                                         annotation (Placement(transformation(extent={{48,168},
            {74,188}},   rotation=0)));

public
    outer Exergy.Utilities.RefEnv refEnv "reference enviroment properties";
    Exergy.Utilities.ViewObjectNE viewObject(nEnergy={4,0,1,0});

    Exergy.Utilities.ViewPort viewOutput
    annotation (Placement(transformation(extent={{86,88},{106,108}})));

equation
  //  viewObject.h[1].E_flow = inlet.m_flow*noEvent(actualStream(inlet.h_outflow));
  viewObject.h[1].E_flow = feedwater.m_flow*fluidIn1.h;
  viewObject.h[1].Ex_flow = feedwater.m_flow*(fluidIn1.h-refEnv.T*fluidIn1.s);

  viewObject.h[2].E_flow = livesteam.m_flow*fluidOut1.h;
  viewObject.h[2].Ex_flow = livesteam.m_flow*(fluidOut1.h-refEnv.T*fluidOut1.s);

    //  viewObject.h[1].E_flow = inlet.m_flow*noEvent(actualStream(inlet.h_outflow));
  viewObject.h[3].E_flow = reheat_in.m_flow*fluidIn2.h;
  viewObject.h[3].Ex_flow = reheat_in.m_flow*(fluidIn2.h-refEnv.T*fluidIn2.s);

  viewObject.h[4].E_flow = reheat_out.m_flow*fluidOut2.h;
  viewObject.h[4].Ex_flow = reheat_out.m_flow*(fluidOut2.h-refEnv.T*fluidOut2.s);

  viewObject.w[1].E_flow = heatRelease.y*Q_flow_F_nom;
  viewObject.w[1].Ex_flow = heatRelease.y*Q_flow_F_nom;

 // viewObject.E[1].E =liveSteam.d*volume_tot_HP *(liveSteam.h);
 // viewObject.E[1].Ex =liveSteam.d*volume_tot_HP *(liveSteam.h);

 // viewObject.E[1].E =mass_HP *(liveSteam.h)  + mass_IP *(reheatedSteam.h);
 // viewObject.E[1].Ex = mass_HP *(liveSteam.h -  refEnv.T*liveSteam.s)
   //                   + mass_IP *(reheatedSteam.h -  refEnv.T*reheatedSteam.s);
//    viewObject.E[1].E = mass_HP *(liveSteam.h - liveSteam.p/liveSteam.d) + mass_IP *(reheatedSteam.h - reheatedSteam.p/reheatedSteam.d);
  // viewObject.E[1].Ex = mass_HP *(liveSteam.h - liveSteam.p/liveSteam.d-  refEnv.T*liveSteam.s)
 //                     + mass_IP *(reheatedSteam.h - reheatedSteam.p/reheatedSteam.d -  refEnv.T*reheatedSteam.s);

   //viewObject.E[1].Ex = shell.summary.fluid.mass * ( refEnv.T*(shell.summary.fluid.rho));

  connect(viewObject.viewOutput,viewOutput);

    h_inHP= if useHomotopy then homotopy(actualStream(feedwater.h_outflow), inStream(feedwater.h_outflow)) else actualStream(feedwater.h_outflow);
    h_outHP= if useHomotopy then homotopy(actualStream(livesteam.h_outflow), h_HP) else actualStream(livesteam.h_outflow);
    h_sprayHP= if useHomotopy then homotopy(actualStream(HPInjection.h_outflow), inStream(HPInjection.h_outflow)) else actualStream(HPInjection.h_outflow);
    h_inIP= if useHomotopy then homotopy(actualStream(reheat_in.h_outflow), inStream(reheat_in.h_outflow)) else actualStream(reheat_in.h_outflow);
    h_outIP= if useHomotopy then homotopy(actualStream(reheat_out.h_outflow), h_IP) else actualStream(reheat_out.h_outflow);
    h_sprayIP= if useHomotopy then homotopy(actualStream(IPInjection.h_outflow), inStream(IPInjection.h_outflow)) else actualStream(IPInjection.h_outflow);

// definition of connector variables
  livesteam.p=p_HP;
  livesteam.h_outflow=liveSteam.h;
  livesteam.xi_outflow=ones(medium.nc-1); //dummy

  HPInjection.p=convert2PressureDrop_HP.y[1]*Delta_p_nomHP*0.15 + p_HP;
  HPInjection.h_outflow = 2000e3;  //This is a generic value (to be refined in the future)
  HPInjection.xi_outflow=ones(medium.nc-1); //dummy

  feedwater.p = convert2PressureDrop_HP.y[1]*Delta_p_nomHP + p_HP;
  feedwater.h_outflow=1330e3; //This is a generic value (to be refined in the future)
  feedwater.xi_outflow=ones(medium.nc-1);//dummy

  reheat_out.p=p_IP;
  reheat_out.h_outflow=reheatedSteam.h;
  reheat_out.xi_outflow=ones(medium.nc-1);//dummy

  IPInjection.p=convert2PressureDrop_IP.y[1]*Delta_p_nomIP*0.5 + p_IP;
  IPInjection.h_outflow = 2000e3;  //This is a generic value (to be refined in the future)
  IPInjection.xi_outflow=ones(medium.nc-1); //dummy

  reheat_in.p= convert2PressureDrop_IP.y[1]*Delta_p_nomIP + p_IP;
  reheat_in.h_outflow=reheatedSteam.h; //This is a generic value (to be refined in the future);
  reheat_in.xi_outflow=ones(medium.nc-1);//dummy

// balance equations - definition of live steam and reheated steam conditions
  drhodt_HP =(feedwater.m_flow+m_flow_heatedHP)/volume_tot_HP;
  mass_HP=if useHomotopy then homotopy(liveSteam.d*volume_tot_HP, volume_tot_HP*rho_nom_HP) else liveSteam.d*volume_tot_HP;

  der(h_HP)=if useHomotopy then homotopy((feedwater.m_flow*h_inHP + m_flow_heatedHP*h_HP + Q_flow_HP + der(p_HP)*volume_tot_HP - drhodt_HP*volume_tot_HP*h_HP)/mass_HP, (m_flow_nomLS*h_inHP + m_flow_nomLS*h_HP + Q_flow_HP + der(p_HP)*volume_tot_HP - drhodt_HP*volume_tot_HP*h_HP)/mass_HP) else (feedwater.m_flow*h_inHP + m_flow_heatedHP*h_HP + Q_flow_HP + der(p_HP)*volume_tot_HP - drhodt_HP*volume_tot_HP*h_HP)/mass_HP;
  der(p_HP)*liveSteam.drhodp_hxi =(drhodt_HP-liveSteam.drhodh_pxi *der(h_HP));

  m_flow_heatedHP=  livesteam.m_flow + HPInjection.m_flow;
  Q_flow_HP=convert2HPFiring.y[1]*calculateEfficiency.y[1]*heatRelease.y*Q_flow_F_nom;

  drhodt_IP =(reheat_in.m_flow+m_flow_heatedIP)/volume_tot_IP;
  mass_IP=if useHomotopy then homotopy(reheatedSteam.d*volume_tot_IP, volume_tot_IP*rho_nom_IP) else reheatedSteam.d*volume_tot_IP;

  der(h_IP)=(reheat_in.m_flow*h_inIP + m_flow_heatedIP*h_IP + Q_flow_IP + der(p_IP)*volume_tot_IP - drhodt_IP*volume_tot_IP*h_IP)/mass_IP;
  der(p_IP)*reheatedSteam.drhodp_hxi = (drhodt_IP-reheatedSteam.drhodh_pxi *der(h_IP));

  m_flow_heatedIP=  reheat_out.m_flow + IPInjection.m_flow;
  Q_flow_IP=(1-convert2HPFiring.y[1])*calculateEfficiency.y[1]*heatRelease.y*Q_flow_F_nom;

  convert2PressureDrop_HP.u[1]=feedwater.m_flow/m_flow_nomLS;
  convert2PressureDrop_IP.u[1]=reheat_in.m_flow/m_flow_nomLS;
  h_evap=liveSteam.VLE.h_v+ 50e3;

//______define eye bus connectors________________________//
  eye_intLS.p=livesteam.p/1e5;
  eye_intLS.h=livesteam.h_outflow/1e3;
  eye_intLS.m_flow=-livesteam.m_flow;
  eye_intLS.T=liveSteam.T-273.15;
  eye_intLS.s=liveSteam.s/1e3;

  eye_intRH.p=reheat_out.p/1e5;
  eye_intRH.h=reheat_out.h_outflow/1e3;
  eye_intRH.m_flow=-reheat_out.m_flow;
  eye_intRH.T=reheatedSteam.T-273.15;
  eye_intRH.s=reheatedSteam.s/1e3;
//___end define eye bus connectors________________________//
initial equation
  heatRelease.y = QF_setl_;
  if initHP == ClaRa.Basics.Choices.Init.steadyEnthalpy then
    der(h_HP) = 0;
  else
    if initHP == ClaRa.Basics.Choices.Init.steadyPressure then
      der(p_HP) = 0;
    else
      if initHP == ClaRa.Basics.Choices.Init.steadyState then
        der(h_HP) = 0;
        der(p_HP) = 0;
      else
        if initHP == ClaRa.Basics.Choices.Init.noInit then
        //do nothing
        else
          assert(false, "Unsupported initialisation option");
        end if;
      end if;
    end if;
  end if;

  if initIP == ClaRa.Basics.Choices.Init.steadyEnthalpy then
    der(h_IP) = 0;
  else
    if initIP == ClaRa.Basics.Choices.Init.steadyPressure then
      der(p_IP) = 0;
    else
      if initIP == ClaRa.Basics.Choices.Init.steadyState then
        der(h_IP) = 0;
        der(p_IP) = 0;
      else
        if initIP == ClaRa.Basics.Choices.Init.noInit then
        //do nothing
        else
          assert(false, "Unsupported initialisation option");
        end if;
      end if;
    end if;
  end if;

equation
  connect(QF_setl_, heatRelease.u) annotation (Line(
      points={{-100,0},{-84,0},{-84,-46},{-68,-46}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(heatRelease.y, convert2HPFiring.u[1]) annotation (Line(
      points={{-45,-46},{-14,-46}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(calculateEfficiency.u[1], heatRelease.y) annotation (Line(
      points={{-14,-74},{-29,-74},{-29,-46},{-45,-46}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(eye_LS,eye_intLS)
                        annotation (Line(
      points={{-40,188},{-40,160},{-42,160}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  connect(eye_RH,eye_intRH)
                        annotation (Line(
      points={{88,188},{88,172},{86,172}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},
            {100,220}})),        Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,220}}), graphics),
    Documentation(info="<html>
<p>
This steam generator is based on lumped balance equations for mass and energy (for each of live steam generation and reheat one set) and on transfer functions and characteristic maps for heat release and pressure losses. The aim of this model is to capture the boiler's main dynamics without allowing deep insight into the process. Thus, the model merges the evaporation and superheat and place the spray injectors (one for each of the live steam generation and the reheater) at the very end of the boiler. This simplification does not allow detailed investigations of the control of steam temperature and feedwater mass flow (applying an enthalpy correction).<br>
</p>

<p>
<strong>Usage advices:</strong><br>
<ul>
<li>The default parameter set do not refer to any real plant</li>
<li>The injector spray is mixed right befor the livesteam connector and the reheat_out connector, respectively</li>
<li>The spray attemperators are ideal, i.e. have no storage behavior</li>
<li>The last superheater bundles are placed before the spray injector</li>
<li>The pressure loss between live steam and HPInjector is 15 % of the HP pressure loss (p_FW - p_LS)</li>
<li>The pressure loss between hot reheated steam and IPInjector is 50 % of the IP pressure loss (p_cRH - p_hRH)</li>
</ul>
</p>
</html>"));
end SteamGenerator_L3;
