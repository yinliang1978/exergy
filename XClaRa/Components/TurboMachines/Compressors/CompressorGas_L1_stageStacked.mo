within Exergy.XClaRa.Components.TurboMachines.Compressors;
model CompressorGas_L1_stageStacked
  "Advanced compressor or fan for ideal gas mixtures using the stage stacking method according to N. Gasparovic"
  import ClaRa;
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

  outer ClaRa.SimCenter simCenter;
extends ClaRa.Basics.Icons.Compressor;
parameter Boolean contributeToCycleSummary = simCenter.contributeToCycleSummary
    "True if component shall contribute to automatic efficiency calculation"                annotation(Dialog(tab="Summary and Visualisation"));
  ClaRa.Basics.Interfaces.Connected2SimCenter connected2SimCenter(
    powerIn=0,
    powerOut=-P_hyd,
    powerAux=-P_shaft + P_hyd) if                                                                                                  contributeToCycleSummary;

  model Outline
   extends ClaRa.Basics.Icons.RecordIcon;
   input SI.VolumeFlowRate V_flow "Volume flow rate";
   input SI.Power P_hyd "Hydraulic power";
   input Real Pi "Pressure ratio";
   input SI.PressureDifference Delta_p "Pressure difference";
   input SI.RPM rpm "Rotational speed";
   input Real Delta_alpha "Angle of VIGV";
   input Real eta_isen "Hydraulic efficiency";
   input Real eta_mech "Mechanic efficiency";
   input Real Y "Specific delivery work";
   input Real m_flow_rel "Relative mass flow";
   input Real Pi_rel "Relative pressure ratio";
  end Outline;

  model Summary
    extends ClaRa.Basics.Icons.RecordIcon;
    Outline outline;
    ClaRa.Basics.Records.FlangeGas  inlet;
    ClaRa.Basics.Records.FlangeGas  outlet;
  end Summary;

  import SI = ClaRa.Basics.Units;

  inner parameter TILMedia.GasTypes.BaseGas    medium = simCenter.flueGasModel;

  final parameter Boolean allow_reverseFlow = false;

  ClaRa.Basics.Interfaces.GasPortIn inlet(Medium=medium, m_flow(final start=
          m_flow_nom, min=if allow_reverseFlow then -Modelica.Constants.inf
           else 1e-5)) "inlet flow"
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
  ClaRa.Basics.Interfaces.GasPortOut outlet(Medium=medium, m_flow(max=if
          allow_reverseFlow then Modelica.Constants.inf else -1e-5))
    "outlet flow"
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));

  Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft if useMechanicalPort
     annotation (Placement(transformation(extent={{-10,90},{10,110}})));
protected
  ClaRa.Components.TurboMachines.Fundamentals.GetInputsRotary getInputsRotary
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={0,20})));
public
  TILMedia.Gas_pT flueGas_inlet( p = inlet.p, T = inStream(inlet.T_outflow), xi = inStream(inlet.xi_outflow), gasType = medium)
    annotation (Placement(transformation(extent={{-90,-12},{-70,8}})));

  TILMedia.Gas_pT flueGas_outlet( gasType = medium, T = T_out, p = outlet.p,  xi = flueGas_inlet.xi)
    annotation (Placement(transformation(extent={{70,-12},{90,8}})));

  Modelica.Blocks.Interfaces.RealInput Delta_alpha_input(value=Delta_alpha) if useExternalVIGVangle
    "VIGV angle input" annotation (Placement(transformation(extent={{-128,60},{-88,100}})));

  //__________________________/ Parameters \_____________________________
  parameter Integer N_stages = 12 "Number of Compressor Stages";
  parameter Integer N_VIGVstages = 1 "Number of VIGV Stages";
  parameter SI.RPM rpm_nom = 3000 "|Nominal Values|Nomial rotational speed";
  parameter Real eta_isen_stage_nom = 0.9
    "|Nominal Values|Nominal isentropic stage efficiency (axial: ca. 0.9, radial: ca. 0.82)";
  parameter Boolean useExternalVIGVangle= false
    "True, if an external source should be used to set VIGV angle";
  parameter SI.Angle Delta_alpha_fixed = 0
    "Fixed angle of VIGV (variable inlet guide vanes)"                                    annotation(Dialog(enable = not useExternalVIGVangle));
  parameter Real eta_mech = 0.99 "Mechanical efficiency";
  parameter Modelica.SIunits.Inertia J "Moment of Inertia" annotation(Dialog(group="Time Response Definitions", enable= not steadyStateTorque));
  parameter Boolean useMechanicalPort=false
    "|Fundamental Definitions|True, if a mechenical flange should be used";
  parameter Boolean steadyStateTorque=false
    "|Fundamental Definitions|True, if steady state mechanical momentum shall be used";
  parameter Boolean useBoundaryAssert=true
    "True, if simulation should stop when surge or choke boundary is hit"  annotation(Dialog(tab = "Advanced"));
  parameter SI.RPM rpm_fixed = 3000 "Constant rotational speed of pump" annotation (Dialog( group = "Fundamental Definitions", enable = not useMechanicalPort));
  parameter SI.Time Tau_aux=0.1 "Time constant of auxilliary kappa states"  annotation(Dialog(tab = "Advanced"));

  parameter SI.MassFlowRate  m_flow_nom = 100
    "|Nominal Values|Nominal mass flow";
  parameter Real Pi_nom = 7 "|Nominal Values|Nominal pressure ratio";
  parameter SI.Temperature T_in_nom = 293.15
    "|Nominal Values|Nominal inlet temperature";
  parameter SI.Pressure p_in_nom = 1.01325e5
    "|Nominal Values|Nominal inlet pressure";
  parameter SI.MassFraction xi_nom[medium.nc - 1]=
     {0,0,0,0,0.76,0.23,0,0,0} "|Nominal Values|Nominal gas composition";
  parameter Boolean useFixedEnthalpyCharacteristic=false
    "|Nominal Values|True, if a fixed nominal stage enthalpy characteristic should be used";
  parameter SI.Length diameter[N_stages] = ones(N_stages)*1.5
    "Individual or mean stage diameter"                                                           annotation(Dialog( group = "Nominal Values", enable = not useFixedEnthalpyCharacteristic));
  parameter Real psi_nom_fixed[N_stages]=ones(N_stages)*0.8
    "Fixed nominal enthalpy stage characteristic (axial: 0.2 to 0.9, radial: 0.9 to 1.4)"
                                                                                           annotation(Dialog( group = "Nominal Values", enable = useFixedEnthalpyCharacteristic));
  parameter String VIGVInfluence= "Lower"
    "Influence of VIGV on psi, phi and eta" annotation(Dialog(group="Parameters"), choices(choice="Higher", choice="Medium",  choice="Lower"));

  //________________________/ Variables \___________________________________
  Real Pi(final start=Pi_nom) "pressure ratio";
  SI.Power P_hyd "Hydraulic power";
  SI.VolumeFlowRate V_flow "Volume flow rate";
  Modelica.SIunits.AngularAcceleration a "Angular acceleration of the shaft";
  SI.Power P_shaft "Mechanical power at shaft";
  SI.RPM rpm "Rotational speed";
  Modelica.SIunits.Torque tau_fluid "Fluid torque";
  //SI.EnthalpyMassSpecific Delta_h;
  Real tau "Overall temperature ratio";
  Real tau_nom "Overall nominal temperature ratio";
  Real eta_isen "Overall isentropic efficiency";
  Real T_out;
  Real kappa "Overall heat capacity ratio";

  //SI.HeatCapacityMassSpecific cp_m;
  Real Delta_alpha "Angle of VIGV (variable inlet guide vanes)";
  Real Delta_alpha_int[i]
    "Angle of VIGV (variable inlet guide vanes) for internal calculation";
  Modelica.SIunits.SpecificEnergy Y_st[i]
    "Specific delivery work of each stage";
  Modelica.SIunits.SpecificEnergy Y "Specific delivery work";

  //______//Stage variables\\______________________________________________
  Real kappa_st[i] "Stage heat capacity ratio";
  Real kappa_in;
  Real kappa_out_st[i];

  Real psi_nom_st[i](each final start=1)
    "Nominal stage enthalpy characteristic";
  Real psi_rel_st[i]( each final start=1)
    "Relative stage enthalpy characteristic";
  Real psi_rel_st_vigv[i]( each final start=1)
    "Relative stage enthalpy characteristic for VIGV";
  Real phi_rel_st[i](each final start=1)
    "Relative stage performance characteristic";
  Real phi_surge_rel_st[i]( each final start=1)
    "Relative stage performance characteristic at surge point";
  // Real phi_max_rel_st[i]
  //   "Maximum relative stage performance characteristic";

  Real eta_isen_st[i]( each final start=eta_isen_stage_nom)
    "Nominal stage efficiency";
  Real eta_isen_rel_st[i](each final start = 1) "Relative stage efficiency";

  Real rpm_corr_st[i]( each final start = (rpm_nom/60)/T_in_nom^0.5)
    "Corrected rotational stage speed";
  Real rpm_corr_nom_st[i]( each final start = (rpm_nom/60)/T_in_nom^0.5)
    "Corrected nominal rotational stage speed";
  Real rpm_corr_rel_st[i](each final start = 1)
    "Relative corrected rotational speed";

  Real m_flow_corr_st[i]( each final start = m_flow_nom * T_in_nom^0.5 / p_in_nom)
    "Corrected stage mass flow";
  Real m_flow_corr_nom_st[i]( each final start = m_flow_nom * T_in_nom^0.5 / p_in_nom)
    "Corrected nominal stage mass flow";
  Real m_flow_corr_rel_st[i](each final start = 1) "Relative stage mass flow";

  Real epsilon_rel_st[i](each final start = 1)
    "Relative stage pressure characteristic";
  Real epsilon_rel_st_max[i](each final start = 1)
    "Relative stage pressure characteristic";
  Real epsilon_rel_st_vigv[i](each final start = 1)
    "Relative stage pressure characteristic for VIGV";

  Real C_1_st[i] "Constant";
  Real C_2_st[i] "Constant";

  Real tau_st[i](each final start = 1.001) "Stage temperature ratio";
  Real Pi_st[i]( each final start = Pi_nom^(1/N_stages)) "Stage pressure ratio";
  Real Pi_st_max[i]( each final start = Pi_nom^(1/N_stages))
    "Maximum stage pressure ratio at surge line";
  Real Pi_prod_st[i](  each final start = 1)
    "Overall pressure ratio until actual stage";                                          //not important, just for evaluation

  SI.Temperature T_out_st[i]( each final start = T_in_nom)
    "Calculated outlet stage temperature";
  SI.Pressure p_out_st[i](each final start = 1e5)
    "Calculated outlet stage pressure";

//Stufenvariablen nominal
  SI.HeatCapacityMassSpecific cp_in_nom
    "Nominal isobaric heat capacity at stage inlet";
  SI.HeatCapacityMassSpecific cv_in_nom
    "Nominal isochoric heat capacity at stage inlet";
  SI.HeatCapacityMassSpecific cp_out_nom_st[i]
    "Nominal isobaric heat capacity at stage outlet";
  SI.HeatCapacityMassSpecific cv_out_nom_st[i]
    "Nominal isochoric heat capacity at stage outlet";
  Real Pi_nom_st[i](each final start = Pi_nom^(1/N_stages))
    "Nominal stage pressure ratio";
  Real Pi_nom_st_vigv[i](each final start = Pi_nom^(1/N_stages))
    "Nominal stage pressure ratio for VIGV";
  Real Pi_prod_nom_st[i]( each final start = 1)
    "Overall nominal pressure ratio until actual stage";
                                                                            //(not important for design point calculation, just for evaluation)
  Real kappa_nom_st[i] "Nominal stage heat capacity ratio";
  Real kappa_nom_st_vigv[i] "Nominal stage heat capacity ratio for VIGV";
  SI.Temperature T_out_nom_st[i](each final start = T_in_nom)
    "Nominal stage outlet temperature";
  SI.Efficiency eta_isen_nom_st[i](each final start=eta_isen_stage_nom)
    "Nominal isentropic stage efficiency";
  Real tau_nom_st[i](each final start = 1.001)
    "Nominal stage temperature ratio";
  Real tau_nom_st_vigv[i](each final start = 1.001)
    "Nominal stage temperature ratio for VIGV";
  SI.EnthalpyMassSpecific h_out_nom_st[i] "Nominal stage outlet enthalpy";
  SI.EnthalpyMassSpecific h_in_nom_st "Nominal compressor inlet enthalpy";
  SI.EnthalpyMassSpecific Delta_h_nom_st[i] "Nominal stage enthalpy difference";
  SI.Pressure p_out_nom_st[i]( each final start = 1e5)
    "Nominal outlet stage pressure";

protected
 Real kappa_aux_st[i]
    "Auxiliary state for kappa (needed when composition changes)";
 Real kappa_aux
    "Auxiliary state for kappa over whole machine (needed when composition changes)";
 Real kappa_nom_aux_st[i]
    "Auxiliary state for kappa (needed when composition changes)";
 Real count_surge[i];
 Real count_choke[i];

 Real eta_isen_od_rpm;
 Real eta_isen_rel_od_rpm(final start = 1);
 Real kappa_out_nom_rp2;
 Real kappa_aux_nom_rp2;
 Real kappa_nom_rp2;
 Real kappa_out_rp2;
 Real kappa_aux_rp2;
 Real kappa_rp2;
 Real T_out_rp2;
 Real tau_rp2;
 Real tau_nom_rp2;
 Real Pi_rp2;
 Real m_flow_corr_rel_rp2;
 Real m_flow_corr_rp2;
 Real rpm_corr_rel_rp2(final start = 1);
 Real phi_rel_rp2(final start = 1);
 Real psi_rel_rp2(final start = 1);
 Real psi_nom_rp2;
 Real epsilon_rel_rp2(final start = 1);
 Real Delta_h_nom_rp2;
 SI.MassFlowRate m_flow_aux(start=m_flow_nom);

//Version 4 mit eta neu (passt sehr gut)
 final parameter Real vigv_coeff_eta_a = if VIGVInfluence=="Lower" then 7e-7 else if VIGVInfluence=="Higher" then 3e-6 else 2e-6;
 final parameter Real vigv_coeff_eta_b = if VIGVInfluence=="Lower" then -0.0003 else if VIGVInfluence=="Higher" then -0.0006 else -0.0004;
 final parameter Real vigv_coeff_eta_c = if VIGVInfluence=="Lower" then 0.0011 else if VIGVInfluence=="Higher" then 0.0014 else 0.0013;
 final parameter Real vigv_coeff_phi_a = if VIGVInfluence=="Lower" then -0.00004 else if VIGVInfluence=="Higher" then -0.0002 else -0.0001;
 final parameter Real vigv_coeff_phi_b = if VIGVInfluence=="Lower" then 0.0176 else if VIGVInfluence=="Higher" then 0.0294 else 0.0228;
 final parameter Real vigv_coeff_psi_a = if VIGVInfluence=="Lower" then -0.0001 else if VIGVInfluence=="Higher" then -0.0002 else -0.0001;
 final parameter Real vigv_coeff_psi_b = if VIGVInfluence=="Lower" then 0.0131 else if VIGVInfluence=="Higher" then 0.0218 else 0.0169;

 final parameter Integer i = N_stages;

     TILMedia.GasObjectFunctions.GasPointer GasPointerRp2OutNom=
       TILMedia.GasObjectFunctions.GasPointer(
       medium.concatGasName,
       8,
       medium.xi_default,
       medium.nc_propertyCalculation,
       medium.nc,
       0,
       0) "Pointer to external medium memory";

     TILMedia.GasObjectFunctions.GasPointer GasPointerRp2Out=
       TILMedia.GasObjectFunctions.GasPointer(
       medium.concatGasName,
       8,
       medium.xi_default,
       medium.nc_propertyCalculation,
       medium.nc,
       0,
       0) "Pointer to external medium memory";

     TILMedia.GasObjectFunctions.GasPointer GasPointerInletNom=
       TILMedia.GasObjectFunctions.GasPointer(
       medium.concatGasName,
       8,
       medium.xi_default,
       medium.nc_propertyCalculation,
       medium.nc,
       0,
       0) "Pointer to external medium memory";

     TILMedia.GasObjectFunctions.GasPointer GasPointerOutletNom[i]=
      {TILMedia.GasObjectFunctions.GasPointer(
       medium.concatGasName,
       8,
       medium.xi_default,
       medium.nc_propertyCalculation,
       medium.nc,
       0,
       0) for dummy in 1:i};

     TILMedia.GasObjectFunctions.GasPointer GasPointerInlet=
       TILMedia.GasObjectFunctions.GasPointer(
       medium.concatGasName,
       8,
       medium.xi_default,
       medium.nc_propertyCalculation,
       medium.nc,
       0,
       0) "Pointer to external medium memory";

     TILMedia.GasObjectFunctions.GasPointer GasPointerOutlet[i]=
      {TILMedia.GasObjectFunctions.GasPointer(
       medium.concatGasName,
       8,
       medium.xi_default,
       medium.nc_propertyCalculation,
       medium.nc,
       0,
       0) for dummy in 1:i};

    TILMedia.GasObjectFunctions.GasPointer GasPointerOutletVigv[i]=
     {TILMedia.GasObjectFunctions.GasPointer(
      medium.concatGasName,
      8,
      medium.xi_default,
      medium.nc_propertyCalculation,
      medium.nc,
      0,
      0) for dummy in 1:i};

public
  inner Summary summary(outline(
   V_flow = V_flow,
   P_hyd = P_hyd,
   Pi = Pi,
   Delta_p = outlet.p - inlet.p,
   rpm = rpm,
   Delta_alpha = Delta_alpha,
   eta_isen = eta_isen,
   eta_mech = eta_mech,
   Y = Y,
   m_flow_rel = m_flow_corr_rel_st[1],
   Pi_rel=Pi/Pi_nom),
    inlet(m_flow = inlet.m_flow,
          T = inStream(inlet.T_outflow),
          p = inlet.p,
          h = flueGas_inlet.h,
          xi = inStream(inlet.xi_outflow),
          H_flow = inlet.m_flow* flueGas_inlet.h),
    outlet(m_flow = -outlet.m_flow,
          T = outlet.T_outflow,
          p = outlet.p,
          h = flueGas_outlet.h,
          xi = outlet.xi_outflow,
          H_flow = -outlet.m_flow* flueGas_outlet.h)) annotation (Placement(transformation(extent={{-100,
            -114},{-80,-94}})));

public
  ClaRa.Basics.Interfaces.EyeOut eyeOut annotation (Placement(
        transformation(extent={{72,-78},{112,-42}}),
        iconTransformation(extent={{92,-70},{112,-50}})));
protected
  ClaRa.Basics.Interfaces.EyeIn eye_int annotation (Placement(
        transformation(extent={{48,-68},{32,-52}}),
        iconTransformation(extent={{90,-84},{84,-78}})));

initial equation

inlet.m_flow = m_flow_corr_st[1]*flueGas_inlet.p/ flueGas_inlet.T^0.5;

     if N_stages > 1 then
      for i in 1:N_stages loop
       kappa_st[i] =  TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(flueGas_inlet.p,flueGas_inlet.T,flueGas_inlet.xi,GasPointerInlet)/TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(flueGas_inlet.p,flueGas_inlet.T,flueGas_inlet.xi,GasPointerInlet);
       kappa_nom_st[i] =  kappa_st[i];
       kappa_nom_st_vigv[i] = kappa_st[i];
      end for;
       kappa = TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(flueGas_inlet.p,flueGas_inlet.T,flueGas_inlet.xi,GasPointerInlet)/TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(flueGas_inlet.p,flueGas_inlet.T,flueGas_inlet.xi,GasPointerInlet);
     else
      for i in 1:N_stages loop
       kappa_st[i] =  TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(flueGas_inlet.p,flueGas_inlet.T,flueGas_inlet.xi,GasPointerInlet)/TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(flueGas_inlet.p,flueGas_inlet.T,flueGas_inlet.xi,GasPointerInlet);
       kappa_nom_st[i] =  kappa_st[i];
       kappa_rp2 = TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(flueGas_inlet.p,flueGas_inlet.T,flueGas_inlet.xi,GasPointerInlet)/TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(flueGas_inlet.p,flueGas_inlet.T,flueGas_inlet.xi,GasPointerInlet);
       kappa_nom_rp2 = kappa_rp2;
       kappa_nom_st_vigv[i] = kappa_st[i];
      end for;
     end if;

equation
//____________________ Mechanics ___________________________
    if useMechanicalPort then
      der(getInputsRotary.rotatoryFlange.phi) = (2*Modelica.Constants.pi*rpm/60);
      J*a + tau_fluid + getInputsRotary.rotatoryFlange.tau = 0
      "Mechanical momentum balance";
    else
      rpm = rpm_fixed;
      getInputsRotary.rotatoryFlange.phi = 0.0;
    end if;

    if (steadyStateTorque) then
      a = 0;
    else
      a = 2*Modelica.Constants.pi/60*der(rpm);
    end if;
    tau_fluid = if noEvent(2*Modelica.Constants.pi*rpm/60<1e-8) then 0 else P_shaft/(2*Modelica.Constants.pi*rpm/60);

  //__________Angle of VIGV_______________
   if (not useExternalVIGVangle) then
     Delta_alpha=Delta_alpha_fixed;
   end if;

for i in 1:N_stages loop
  Delta_alpha_int[i] = if i <= N_VIGVstages then Delta_alpha else 0;
end for;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////DESIGN POINT////DESIGN POINT////DESIGN POINT////DESIGN POINT////DESIGN POINT////DESIGN POINT////DESIGN POINT////DESIGN POINT////DESIGN POINT////DESIGN POINT//////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
if N_stages == 1 then  ////FOR A SINGLE STAGE ONLY////FOR A SINGLE STAGE ONLY////FOR A SINGLE STAGE ONLY////FOR A SINGLE STAGE ONLY////FOR A SINGLE STAGE ONLY/////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    cp_in_nom = TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(p_in_nom,T_in_nom,xi_nom,GasPointerInletNom);
    cv_in_nom = TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(p_in_nom,T_in_nom,xi_nom,GasPointerInletNom);
    cp_out_nom_st[1] = TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(p_out_nom_st[1],T_out_nom_st[1],xi_nom,GasPointerOutletNom[1]);
    cv_out_nom_st[1] = TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(p_out_nom_st[1],T_out_nom_st[1],xi_nom,GasPointerOutletNom[1]);
    kappa_nom_aux_st[1] = (cp_in_nom/cv_in_nom + cp_out_nom_st[1]/cv_out_nom_st[1]) / 2;
    der(kappa_nom_st[1]) = 1/Tau_aux*(kappa_nom_aux_st[1]-kappa_nom_st[1]);
    //kappa_nom_st[1] =1.4;

    tau_nom_st[1] = 1 + (Pi_nom_st[1]^((kappa_nom_st[1]-1)/kappa_nom_st[1])-1) * 1/eta_isen_nom_st[1];
    T_out_nom_st[1] = tau_nom_st[1] * T_in_nom; //eq. 18
    p_out_nom_st[1] = p_in_nom * Pi_nom;
    h_in_nom_st = TILMedia.GasObjectFunctions.specificEnthalpy_pTxi(p_in_nom,T_in_nom,xi_nom,GasPointerInletNom);
    h_out_nom_st[1] = TILMedia.GasObjectFunctions.specificEnthalpy_pTxi(p_out_nom_st[1],T_out_nom_st[1],xi_nom,GasPointerOutletNom[1]);
    Delta_h_nom_st[1] = h_out_nom_st[1] - h_in_nom_st;

    //_______/Alternative design point for efficiency on new rpm-curve\________________________________________
    //See also Goettlich for documentation of this corrected calculation for a single stage only

    eta_isen_nom_st[1] = eta_isen_stage_nom;

    Pi_prod_nom_st[1] = Pi_nom; //Overall pressure ratio until actual stage
    Pi_nom_st[1] = Pi_nom;

    tau_nom_st[1] = tau_nom;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
else   ////FIRST STAGE////FIRST STAGE////FIRST STAGE////FIRST STAGE////FIRST STAGE////FIRST STAGE////FIRST STAGE////FIRST STAGE////FIRST STAGE////FIRST STAGE/////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   cp_in_nom = TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(p_in_nom,T_in_nom,xi_nom,GasPointerInletNom);
   cv_in_nom = TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(p_in_nom,T_in_nom,xi_nom,GasPointerInletNom);
   cp_out_nom_st[1] = TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(p_out_nom_st[1],T_out_nom_st[1],xi_nom,GasPointerOutletNom[1]);
   cv_out_nom_st[1] = TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(p_out_nom_st[1],T_out_nom_st[1],xi_nom,GasPointerOutletNom[1]);
   kappa_nom_aux_st[1] = (cp_in_nom/cv_in_nom + cp_out_nom_st[1]/cv_out_nom_st[1]) / 2;
   der(kappa_nom_st[1]) = 1/Tau_aux*(kappa_nom_aux_st[1]-kappa_nom_st[1]); //Auxiliary state for kappa (needed when composition changes)
   //kappa_nom_st[1] =1.4;

   tau_nom_st[1] = 1 + ((p_out_nom_st[1] / p_in_nom)^((kappa_nom_st[1]-1)/kappa_nom_st[1])-1) * 1/eta_isen_nom_st[1];
   T_out_nom_st[1] = tau_nom_st[1] * T_in_nom; //eq. 18
   p_out_nom_st[1] = p_in_nom * Pi_nom_st[1];
   h_in_nom_st = TILMedia.GasObjectFunctions.specificEnthalpy_pTxi(p_in_nom,T_in_nom,xi_nom,GasPointerInletNom);
   h_out_nom_st[1] = TILMedia.GasObjectFunctions.specificEnthalpy_pTxi(p_out_nom_st[1],T_out_nom_st[1],xi_nom,GasPointerOutletNom[1]);
   Delta_h_nom_st[1] = h_out_nom_st[1] - h_in_nom_st;

   eta_isen_nom_st[1] = eta_isen_stage_nom;

   Pi_prod_nom_st[1] = Pi_nom_st[1]; //Overall pressure ratio until actual stage

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   for i in 2:N_stages loop  ////STAGE LOOP////STAGE LOOP////STAGE LOOP////STAGE LOOP////STAGE LOOP////STAGE LOOP////STAGE LOOP////STAGE LOOP////STAGE LOOP////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

       cp_out_nom_st[i] = TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(p_out_nom_st[i],T_out_nom_st[i],xi_nom,GasPointerOutletNom[i]);
       cv_out_nom_st[i] = TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(p_out_nom_st[i],T_out_nom_st[i],xi_nom,GasPointerOutletNom[i]);
       kappa_nom_aux_st[i] = (cp_out_nom_st[i-1]/cv_out_nom_st[i-1] + cp_out_nom_st[i]/cv_out_nom_st[i]) / 2;
       der(kappa_nom_st[i]) = 1/Tau_aux*(kappa_nom_aux_st[i]-kappa_nom_st[i]); //Auxiliary state for kappa (needed when composition changes)
       //kappa_nom_st[i] =1.4;

       tau_nom_st[i] = 1 + (Pi_nom_st[i]^((kappa_nom_st[i]-1)/kappa_nom_st[i])-1) * 1/eta_isen_nom_st[i]; //eq. 25
       T_out_nom_st[i] = tau_nom_st[i] * T_out_nom_st[i-1]; //eq. 18
       Pi_nom_st[i] = p_out_nom_st[i] / p_out_nom_st[i-1];
       h_out_nom_st[i] = TILMedia.GasObjectFunctions.specificEnthalpy_pTxi(p_out_nom_st[i],T_out_nom_st[i],xi_nom,GasPointerOutletNom[i]);
       Delta_h_nom_st[i] = h_out_nom_st[i] - h_out_nom_st[i-1];

       eta_isen_nom_st[i] = eta_isen_stage_nom;

       Pi_prod_nom_st[i] = product(Pi_nom_st[j] for j in 1:i); //Overall pressure ratio until actual stage

       tau_nom_st[i]=2-1/tau_nom_st[i-1];

   end for;

   Pi_nom = product(Pi_nom_st);
   tau_nom = product(tau_nom_st);

end if;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////OFF DESIGN////OFF DESIGN////OFF DESIGN////OFF DESIGN////OFF DESIGN////OFF DESIGN////OFF DESIGN////OFF DESIGN////OFF DESIGN////OFF DESIGN////OFF DESIGN////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 if N_stages == 1 then ////FOR A SINGLE STAGE ONLY////FOR A SINGLE STAGE ONLY////FOR A SINGLE STAGE ONLY////FOR A SINGLE STAGE ONLY////FOR A SINGLE STAGE ONLY/////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      kappa_in = TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(inlet.p,inStream(inlet.T_outflow),inStream(inlet.xi_outflow),GasPointerInlet)/TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(inlet.p,inStream(inlet.T_outflow),inStream(inlet.
      xi_outflow),                                                                                                    GasPointerInlet);
      kappa_out_st[1] = TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(outlet.p,outlet.T_outflow,inStream(inlet.xi_outflow),GasPointerOutlet[1])/TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(outlet.p,outlet.T_outflow,inStream(inlet.
      xi_outflow),                                                                                                    GasPointerOutlet[1]);
      kappa_aux_st[1] = (kappa_in + kappa_out_st[1])/2.0;
      der(kappa_st[1]) = 1/Tau_aux*(kappa_aux_st[1]-kappa_st[1]);
      //kappa_st[1] =1.4;

      //_____________/Corrected RPM\_______________________________________________________
      rpm_corr_nom_st[1] = (rpm_nom/60)/T_in_nom^0.5;
      rpm_corr_st[1] = (rpm/60)/inStream(inlet.T_outflow)^0.5;
      rpm_corr_rel_st[1] = rpm_corr_st[1]/rpm_corr_nom_st[1];

      //_____________/Corrected mass flow\_________________________________________________
      m_flow_corr_nom_st[1] = m_flow_nom * T_in_nom^0.5 / p_in_nom * (1 + vigv_coeff_phi_a*Delta_alpha_int[1]^2 + vigv_coeff_phi_b*Delta_alpha_int[1]);//* (1 + 0.014 * Delta_alpha_int[1]); //eq. 33 (multiplied with additional inlet guide vane function)
      m_flow_corr_st[1] = m_flow_aux * inStream(inlet.T_outflow)^0.5 / inlet.p;
      m_flow_corr_rel_st[1] = m_flow_corr_st[1]/(m_flow_corr_nom_st[1]);

      //_____________/Stage performance characteristic\____________________________________
      phi_rel_st[1] = m_flow_corr_rel_st[1]/rpm_corr_rel_st[1]; //eq. 34 //1.0 wg behandlung als Bezugsdrehzahl

       //_____________/Additional values for stage enthalpy characteristic\_________________
       if useFixedEnthalpyCharacteristic == false then
         psi_nom_st[1] = Delta_h_nom_st[1]/(Modelica.Constants.pi^2 * diameter[1]^2 * (rpm_nom/60)^2/2) * (1 + vigv_coeff_psi_a*Delta_alpha_int[1]^2 + vigv_coeff_psi_b*Delta_alpha_int[1]);//(1 + 0.014 * Delta_alpha_int[1]);//eq. 19 (multiplied with additional inlet guide vane function)
       else
         psi_nom_st[1] = psi_nom_fixed[1] * (1 + vigv_coeff_psi_a*Delta_alpha_int[1]^2 + vigv_coeff_psi_b*Delta_alpha_int[1]);
       end if;

      C_1_st[1] = 2/psi_nom_st[1];
      C_2_st[1] = 1 - C_1_st[1];

      //_____________/Stage enthalpy characteristic\_______________________________________
      psi_rel_st[1] = C_1_st[1] + C_2_st[1] * phi_rel_st[1]; //eq. 15

      //_____________/Stage pressure characteristic\_______________________________________
      epsilon_rel_st[1] = (1 + C_1_st[1] - C_1_st[1] * phi_rel_st[1]) * phi_rel_st[1]; //eq. 39

      //_____________/Reference point adjustion for VIGV\________________________________________________
      //psi_rel_st_vigv[1] = (C_1_st[1] + C_2_st[1]) * (1 + 0.014 * Delta_alpha_int[1]);
      psi_rel_st_vigv[1] = 1 + vigv_coeff_psi_a*Delta_alpha_int[1]^2 + vigv_coeff_psi_b*Delta_alpha_int[1]; // 1 + 0.014 * Delta_alpha_int[1];
      epsilon_rel_st_vigv[1] =(-C_1_st[1]*(psi_rel_st_vigv[1]^2+1) + (C_1_st[1]^2+1)* psi_rel_st_vigv[1])/C_2_st[1]^2;

      der(kappa_nom_st_vigv[1]) =  1/Tau_aux*((TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(p_in_nom,T_in_nom,xi_nom,GasPointerInletNom)/TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(p_in_nom,T_in_nom,xi_nom,GasPointerInletNom)
                        + TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(p_in_nom*Pi_nom_st_vigv[1],T_in_nom*tau_nom_st_vigv[1],xi_nom,GasPointerOutletVigv[1])/TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(p_in_nom*Pi_nom_st_vigv[1],T_in_nom*tau_nom_st_vigv[1],xi_nom,GasPointerOutletVigv[1]))/2 -kappa_nom_st_vigv[1]);

      // tau_nom_st_vigv[1] = if N_VIGVstages > 0 then 1 + (Pi_nom_st_vigv[1]^((kappa_nom_st[1]-1)/kappa_nom_st[1])-1) * 1/(eta_isen_od_rpm *min(1 + vigv_coeff_eta_a*Delta_alpha_int[1]^3 + vigv_coeff_eta_b*Delta_alpha_int[1]^2 + vigv_coeff_eta_c*Delta_alpha_int[1],1)) else tau_nom_st[1]; //(1 - 0.000111 * Delta_alpha_int[1]^2)
      // Pi_nom_st_vigv[1] = if N_VIGVstages > 0 then (1 +  epsilon_rel_st_vigv[1] * (Pi_nom_st[1]^((kappa_nom_st[1]-1)/kappa_nom_st[1])-1))^(kappa_nom_st[1]/(kappa_nom_st[1]-1)) else Pi_nom_st[1];

       tau_nom_st_vigv[1] = if N_VIGVstages > 0 then 1 + (Pi_nom_st_vigv[1]^((kappa_nom_st_vigv[1]-1)/kappa_nom_st_vigv[1])-1) * 1/(eta_isen_od_rpm *min(1 + vigv_coeff_eta_a*Delta_alpha_int[1]^3 + vigv_coeff_eta_b*Delta_alpha_int[1]^2 + vigv_coeff_eta_c*Delta_alpha_int[1],1)) else tau_nom_st[1]; //(1 - 0.000111 * Delta_alpha_int[1]^2)
       Pi_nom_st_vigv[1] = if N_VIGVstages > 0 then (1 +  epsilon_rel_st_vigv[1] * (Pi_nom_st[1]^((kappa_nom_st[1]-1)/kappa_nom_st[1])-1))^(kappa_nom_st_vigv[1]/(kappa_nom_st_vigv[1]-1)) else Pi_nom_st[1];

      //_____________/Stage temperature ratio\_____________________________________________
      tau_st[1] = 1 + psi_rel_st[1] * rpm_corr_rel_st[1]^2 * (tau_nom_st_vigv[1]-1);  //eq. 22
      //tau_st[1] = 1 + (Pi_st[1]^((kappa_st[1]-1)/kappa_st[1])-1) * 1/eta_isen_st[1];
      T_out_st[1] = tau_st[1] * flueGas_inlet.T; //eq. 18

      //_____________/Stage pressure ratio\________________________________________________
      //Pi_st[1] = (1 + epsilon_rel_st[1] * rpm_corr_rel_st[1]^2 * (Pi_nom_st_vigv[1]^((kappa_st[1]-1)/kappa_st[1])-1))^(kappa_st[1]/(kappa_st[1]-1)); //eq. 31
      epsilon_rel_st[1] = 1/rpm_corr_rel_st[1]^2 * ((Pi_st[1]^((kappa_st[1]-1)/kappa_st[1])-1)/(Pi_nom_st_vigv[1]^((kappa_nom_st_vigv[1]-1)/kappa_nom_st_vigv[1])-1)); //eq. 31
      Pi_st[1] = outlet.p/inlet.p;

      //_____________/Isentropic stage efficiency\_________________________________________
      eta_isen_st[1] = (Pi_st[1]^((kappa_st[1]-1)/kappa_st[1])-1) / (tau_st[1] - 1); //eq. 23
      eta_isen_rel_st[1] = eta_isen_st[1]/(eta_isen_od_rpm * min(1 + vigv_coeff_eta_a*Delta_alpha_int[1]^3 + vigv_coeff_eta_b*Delta_alpha_int[1]^2 + vigv_coeff_eta_c*Delta_alpha_int[1],1));
      //eta_isen_st[1] = eta_isen_rel_st[1]*(eta_isen_od_rpm * min(1 - 0.0004*Delta_alpha_int[1]^2 + 0.0007*Delta_alpha_int[1],1));//(1 - 0.000111 * Delta_alpha_int[1]^2)
      //eta_isen_rel_st[1] = epsilon_rel_st[1]/psi_rel_st[1];

       p_out_st[1] = inlet.p * Pi_st[1];
       T_out = T_out_st[1];
       Pi = Pi_st[1];
       tau = tau_st[1];
       kappa = kappa_st[1];
       kappa_aux=kappa;
       eta_isen = eta_isen_st[1];
       Pi_prod_st[1] = Pi_st[1];

       Y_st[1]= (p_out_st[1] - inlet.p)/(0.5*(TILMedia.GasObjectFunctions.density_pTxi(inlet.p,inStream(inlet.T_outflow),inStream(inlet.xi_outflow),GasPointerInlet) + TILMedia.GasObjectFunctions.density_pTxi(outlet.p,outlet.T_outflow,inStream(inlet.
      xi_outflow),                                                                                                    GasPointerOutlet[1])));

   //____________/Surge line\____________________________________________________________
   if psi_nom_st[i] <= 0.9 then
     phi_surge_rel_st[1] = 1/3 +  (2 + psi_nom_st[1])/6; //eq. 66 and 40 (axial machine)
   else
     phi_surge_rel_st[1] = 2/7 + 5/7 * (2 + psi_nom_st[1])/4; //eq. from goettlich and 40 (radial machine)
   end if;

   epsilon_rel_st_max[1] = (1 + C_1_st[1] - C_1_st[1]*phi_surge_rel_st[1])*phi_surge_rel_st[1];
   Pi_st_max[1] = (1 + epsilon_rel_st_max[1]*rpm_corr_rel_st[1]^2*(Pi_nom_st_vigv[1]^((kappa_nom_st_vigv[1]-1)/kappa_nom_st_vigv[1])-1))^(kappa_st[1]/(kappa_st[1]-1));

   /////////////////////////////////////////////////////////////////////////////////////////////////////////
   //___________/Auxilliary reference point (rp2) for calculation of new reference efficiency at off design speed
   //This correction is only necessary if a single compressor stage is calculated.
   //For this correction a new reference point (rp2, on nominal speed line) is calculated for each off design speed line.
   //See also Goettlich for documentation of this corrected calculation for a single stage only.
   /////////////////////////////////////////////////////////////////////////////////////////////////////////

   rpm_corr_rel_rp2 = 1.0;
   m_flow_corr_rel_rp2 = rpm_corr_rel_st[1];

   m_flow_corr_rp2 = m_flow_corr_rel_rp2 * m_flow_nom * T_in_nom^0.5 / p_in_nom;// * (1 - 0.0005*Delta_alpha_int[1]^2 + 0.0231*Delta_alpha_int[1]);//(1 + 0.014 * Delta_alpha_int[1])

   phi_rel_rp2 = rpm_corr_rel_st[1]; // Goettlich eq. 6.12 (error in Goettlich betw. eq. 6.11 and 6.12: phi_rel_1 = 1.0)

   if useFixedEnthalpyCharacteristic == false then
     psi_nom_rp2 = Delta_h_nom_rp2/(Modelica.Constants.pi^2 * diameter[1]^2 * (rpm_nom/60)^2/2);// * (1 - 0.0003*Delta_alpha_int[1]^2 + 0.0176*Delta_alpha_int[1]);//(1 + 0.014 * Delta_alpha_int[1])
   else
     psi_nom_rp2 = psi_nom_fixed[1]; //* (1 - 0.0003*Delta_alpha_int[1]^2 + 0.0176*Delta_alpha_int[1]);//(1 + 0.014 * Delta_alpha_int[1])
   end if;

   epsilon_rel_rp2 = (1 + (2/psi_nom_rp2) - (2/psi_nom_rp2) * phi_rel_rp2) * phi_rel_rp2;
   psi_rel_rp2 = (2/psi_nom_rp2) + (1- 2/psi_nom_rp2) * phi_rel_rp2;

   //eta_isen_rel_od_rpm = epsilon_rel_rp2/psi_rel_rp2;
   //eta_isen_od_rpm = eta_isen_rel_od_rpm * eta_isen_stage_nom;

   tau_nom_rp2 = 1 + (Pi_nom^((kappa_nom_rp2-1)/kappa_nom_rp2)-1) * 1/eta_isen_stage_nom;//* min(1 - 0.0004*Delta_alpha_int[1]^2 + 0.0007*Delta_alpha_int[1],1));//*(1 - 0.000111 * Delta_alpha_int[1]^2)
   tau_rp2 = 1 + psi_rel_rp2 * rpm_corr_rel_rp2^2 * (tau_nom_rp2-1);
   T_out_rp2 = tau_rp2 * flueGas_inlet.T;
   Delta_h_nom_rp2 = TILMedia.GasObjectFunctions.specificEnthalpy_pTxi(p_in_nom*Pi_nom,T_in_nom*tau_nom_rp2,xi_nom,GasPointerRp2OutNom) - h_in_nom_st;

   kappa_out_nom_rp2 = TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(p_in_nom*Pi_nom,T_in_nom*tau_nom_rp2,xi_nom,GasPointerRp2OutNom)/TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(p_in_nom*Pi_nom,T_in_nom*tau_nom_rp2,xi_nom,GasPointerRp2OutNom);
   kappa_aux_nom_rp2 = (cp_in_nom/cv_in_nom + kappa_out_nom_rp2)/2.0;
   der(kappa_nom_rp2) = 1/Tau_aux*(kappa_aux_nom_rp2-kappa_nom_rp2);
   //kappa_nom_rp2 =1.4;

   kappa_out_rp2 = TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(inlet.p*Pi_rp2,T_out_rp2,inStream(inlet.xi_outflow),GasPointerRp2Out)/TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(inlet.p*Pi_rp2,T_out_rp2,inStream(inlet.
      xi_outflow),                                                                                                    GasPointerRp2Out);
   kappa_aux_rp2 = (kappa_in + kappa_out_rp2)/2.0;
   der(kappa_rp2) = 1/Tau_aux*(kappa_aux_rp2-kappa_rp2);
   //kappa_rp2 =1.4;

   Pi_rp2 = (1 + epsilon_rel_rp2 * rpm_corr_rel_rp2^2 * (Pi_nom^((kappa_rp2-1)/kappa_rp2)-1))^(kappa_rp2/(kappa_rp2-1));
   eta_isen_od_rpm = (Pi_rp2^((kappa_rp2-1)/kappa_rp2)-1) / (tau_rp2 - 1);
   eta_isen_rel_od_rpm = eta_isen_od_rpm/eta_isen_stage_nom;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  else ////FIRST STAGE////FIRST STAGE////FIRST STAGE////FIRST STAGE////FIRST STAGE////FIRST STAGE////FIRST STAGE////FIRST STAGE////FIRST STAGE////FIRST STAGE/////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    kappa_in = TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(inlet.p,inStream(inlet.T_outflow),inStream(inlet.xi_outflow),GasPointerInlet)/TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(inlet.p,inStream(inlet.T_outflow),inStream(inlet.
      xi_outflow),                                                                                                    GasPointerInlet);
    kappa_out_st[1] = TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(p_out_st[1],T_out_st[1],inStream(inlet.xi_outflow),GasPointerOutlet[1])/TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(p_out_st[1],T_out_st[1],inStream(inlet.
      xi_outflow),                                                                                                    GasPointerOutlet[1]);
    kappa_aux_st[1] = (kappa_in + kappa_out_st[1])/2.0;
    der(kappa_st[1]) = 1/Tau_aux*(kappa_aux_st[1]-kappa_st[1]);
    //kappa_st[1]=1.4;

    //_____________/Corrected RPM\_______________________________________________________
    rpm_corr_nom_st[1] = (rpm_nom/60)/T_in_nom^0.5;
    rpm_corr_st[1] = (rpm/60)/inStream(inlet.T_outflow)^0.5;
    rpm_corr_rel_st[1] = rpm_corr_st[1]/rpm_corr_nom_st[1];

    //_____________/Corrected mass flow\_________________________________________________

    m_flow_corr_nom_st[1] = m_flow_nom * T_in_nom^0.5 / p_in_nom * (1 + vigv_coeff_phi_a*Delta_alpha_int[1]^2 + vigv_coeff_phi_b*Delta_alpha_int[1]);//(1 + 0.014 * Delta_alpha_int[1]); //eq. 33 (multiplied with additional inlet guide vane function);
    m_flow_corr_st[1] = m_flow_aux * inStream(inlet.T_outflow)^0.5 / inlet.p;
    m_flow_corr_rel_st[1] = m_flow_corr_st[1]/(m_flow_corr_nom_st[1]);

    //_____________/Stage performance characteristic\____________________________________
    phi_rel_st[1] = (m_flow_corr_rel_st[1]/rpm_corr_rel_st[1]); //eq. 34

    //_____________/Additional values for stage enthalpy characteristic\_________________
    if useFixedEnthalpyCharacteristic == false then
        psi_nom_st[1] = Delta_h_nom_st[1]/(Modelica.Constants.pi^2 * diameter[1]^2 * (rpm_nom/60)^2/2) * (1 + vigv_coeff_psi_a*Delta_alpha_int[1]^2 + vigv_coeff_psi_b*Delta_alpha_int[1]);//(1 + 0.014 * Delta_alpha_int[1]);//eq. 19 (multiplied with additional inlet guide vane function)
    else
        psi_nom_st[1] = psi_nom_fixed[1] * (1 + vigv_coeff_psi_a*Delta_alpha_int[1]^2 + vigv_coeff_psi_b*Delta_alpha_int[1]);//(1 + 0.014 * Delta_alpha_int[1]);
    end if;

    C_1_st[1] = 2/psi_nom_st[1];
    C_2_st[1] = 1 - C_1_st[1];

    //_____________/Stage enthalpy characteristic\_______________________________________
    psi_rel_st[1] = (C_1_st[1] + C_2_st[1] * phi_rel_st[1]);//* (1 + 0.014 * Delta_alpha_int[1]); //eq. 15

    //_____________/Stage pressure characteristic\________________________________________________
    epsilon_rel_st[1] = (1 + C_1_st[1] - C_1_st[1] * phi_rel_st[1]) * phi_rel_st[1]; //eq. 39

    //_____________/Reference point adjustion for VIGV\________________________________________________
    //psi_rel_st_vigv[1] = (C_1_st[1] + C_2_st[1]) * (1 + 0.014 * Delta_alpha_int[1]);
    psi_rel_st_vigv[1] = 1 + vigv_coeff_psi_a*Delta_alpha_int[1]^2 + vigv_coeff_psi_b*Delta_alpha_int[1];//1 + 0.014 * Delta_alpha_int[1];
    epsilon_rel_st_vigv[1] =(-C_1_st[1]*(psi_rel_st_vigv[1]^2+1) + (C_1_st[1]^2+1)* psi_rel_st_vigv[1])/C_2_st[1]^2;

    der(kappa_nom_st_vigv[1]) =  1/Tau_aux*((TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(p_in_nom,T_in_nom,xi_nom,GasPointerInletNom)/TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(p_in_nom,T_in_nom,xi_nom,GasPointerInletNom)
                              + TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(p_in_nom*Pi_nom_st_vigv[1],T_in_nom*tau_nom_st_vigv[1],xi_nom,GasPointerOutletVigv[1])/TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(p_in_nom*Pi_nom_st_vigv[1],T_in_nom*tau_nom_st_vigv[1],xi_nom,GasPointerOutletVigv[1]))/2 -kappa_nom_st_vigv[1]);

    //tau_nom_st_vigv[1] = if N_VIGVstages > 0 then 1 + (Pi_nom_st_vigv[1]^((kappa_nom_st[1]-1)/kappa_nom_st[1])-1) * 1/(eta_isen_nom_st[1] * min(1 + vigv_coeff_eta_a*Delta_alpha_int[1]^3 + vigv_coeff_eta_b*Delta_alpha_int[1]^2 + vigv_coeff_eta_c*Delta_alpha_int[1],1)) else tau_nom_st[1];
    //Pi_nom_st_vigv[1] = if N_VIGVstages > 0 then (1 +  epsilon_rel_st_vigv[1] * (Pi_nom_st[1]^((kappa_nom_st[1]-1)/kappa_nom_st[1])-1))^(kappa_nom_st[1]/(kappa_nom_st[1]-1)) else Pi_nom_st[1];

    tau_nom_st_vigv[1] = if N_VIGVstages > 0 then 1 + (Pi_nom_st_vigv[1]^((kappa_nom_st_vigv[1]-1)/kappa_nom_st_vigv[1])-1) * 1/(eta_isen_nom_st[1] * min(1 + vigv_coeff_eta_a*Delta_alpha_int[1]^3 + vigv_coeff_eta_b*Delta_alpha_int[1]^2 + vigv_coeff_eta_c*Delta_alpha_int[1],1)) else tau_nom_st[1];
    Pi_nom_st_vigv[1] = if N_VIGVstages > 0 then (1 +  epsilon_rel_st_vigv[1] * (Pi_nom_st[1]^((kappa_nom_st[1]-1)/kappa_nom_st[1])-1))^(kappa_nom_st_vigv[1]/(kappa_nom_st_vigv[1]-1)) else Pi_nom_st[1];

    //_____________/Stage temperature ratio\_____________________________________________
    //psi_rel = 1/rpm_corr_rel^2*(tau-1)/(tau_nom-1); //eq. 22
    tau_st[1] = 1 + psi_rel_st[1] * rpm_corr_rel_st[1]^2 * (tau_nom_st_vigv[1]-1);
    //tau_st[1] = 1 + (Pi_st[1]^((kappa_st[1]-1)/kappa_st[1])-1) * 1/eta_isen_st[1];
    T_out_st[1] = tau_st[1] * flueGas_inlet.T; //eq. 18

    //_____________/Stage pressure ratio\________________________________________________
    //Pi_st[1] = (1 + epsilon_rel_st[1] * rpm_corr_rel_st[1]^2 * (Pi_nom_st[1]^((kappa_st[1]-1)/kappa_st[1])-1))^(kappa_st[1]/(kappa_st[1]-1)); //eq. 31
    epsilon_rel_st[1] = 1/rpm_corr_rel_st[1]^2 * ((Pi_st[1]^((kappa_st[1]-1)/kappa_st[1])-1)/(Pi_nom_st_vigv[1]^((kappa_nom_st_vigv[1]-1)/kappa_nom_st_vigv[1])-1)); //eq. 31
    p_out_st[1] = inlet.p * Pi_st[1];

    //_____________/Isentropic stage efficiency\_________________________________________
    eta_isen_st[1] =(Pi_st[1]^((kappa_st[1]-1)/kappa_st[1])-1) / (tau_st[1] - 1); //eq. 23
    eta_isen_rel_st[1] = eta_isen_st[1]/(eta_isen_nom_st[1]*min(1 + vigv_coeff_eta_a*Delta_alpha_int[1]^3 + vigv_coeff_eta_b*Delta_alpha_int[1]^2 + vigv_coeff_eta_c*Delta_alpha_int[1],1));
    //eta_isen_st[1] = eta_isen_rel_st[1]*(eta_isen_nom_st[1]*min(1 - 0.0004*Delta_alpha_int[1]^2 + 0.0007*Delta_alpha_int[1],1));
    //eta_isen_rel_st[1] = epsilon_rel_st[1]/psi_rel_st[1];

    Pi_prod_st[1] = Pi_st[1];

    Y_st[1]= (p_out_st[1] - inlet.p)/(0.5*(TILMedia.GasObjectFunctions.density_pTxi(inlet.p,inStream(inlet.T_outflow),inStream(inlet.xi_outflow),GasPointerInlet) + TILMedia.GasObjectFunctions.density_pTxi(p_out_st[1],T_out_st[1],inStream(inlet.
      xi_outflow),                                                                                                    GasPointerOutlet[1])));

    //____________/Surge line\____________________________________________________________
    if psi_nom_st[i] <= 0.9 then
      phi_surge_rel_st[1] = 1/3 +  (2 + psi_nom_st[1])/6; //eq. 66 and 40 (axial machine)
    else
      phi_surge_rel_st[1] = 2/7 + 5/7 * (2 + psi_nom_st[1])/4; //eq. from goettlich and 40 (radial machine)
    end if;

   epsilon_rel_st_max[1] = (1 + C_1_st[1] - C_1_st[1]*phi_surge_rel_st[1])*phi_surge_rel_st[1];
   Pi_st_max[1] = (1 + epsilon_rel_st_max[1]*rpm_corr_rel_st[1]^2*(Pi_nom_st_vigv[1]^((kappa_nom_st_vigv[1]-1)/kappa_nom_st_vigv[1])-1))^(kappa_st[1]/(kappa_st[1]-1));

    //___________/Dummies (only needed as corr. of eta_isen during calulation of a single stage)\___________
    eta_isen_od_rpm=0;
    eta_isen_rel_od_rpm=0;
    kappa_out_nom_rp2=0;
    kappa_aux_nom_rp2=0;
    kappa_nom_rp2=0;
    kappa_out_rp2=0;
    kappa_aux_rp2=0;
    kappa_rp2=0;
    T_out_rp2=0;
    tau_rp2=0;
    tau_nom_rp2=0;
    Pi_rp2=0;
    m_flow_corr_rel_rp2=0;
    rpm_corr_rel_rp2=0;
    phi_rel_rp2=0;
    psi_rel_rp2=0;
    psi_nom_rp2=0;
    epsilon_rel_rp2=0;
    m_flow_corr_rp2=0;
    Delta_h_nom_rp2=0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 for i in 2:N_stages loop ////STAGE LOOP////STAGE LOOP////STAGE LOOP////STAGE LOOP////STAGE LOOP////STAGE LOOP////STAGE LOOP////STAGE LOOP////STAGE LOOP///////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     kappa_out_st[i] = TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(p_out_st[i],T_out_st[i],inStream(inlet.xi_outflow),GasPointerOutlet[i])/TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(p_out_st[i],T_out_st[i],inStream(inlet.
        xi_outflow),                                                                                                    GasPointerOutlet[i]);
     kappa_aux_st[i] = (kappa_out_st[i-1] + kappa_out_st[i])/2.0;
     der(kappa_st[i]) = 1/Tau_aux*(kappa_aux_st[i]-kappa_st[i]);
     //kappa_st[i]=1.4;

     //_____________/Corrected RPM\_______________________________________________________
     rpm_corr_nom_st[i] = (rpm_nom/60)/T_out_nom_st[i-1]^0.5;
     rpm_corr_st[i] = (rpm/60)/T_out_st[i-1]^0.5;
     rpm_corr_rel_st[i] = rpm_corr_st[i]/(rpm_corr_nom_st[i]);

     //_____________/Corrected mass flow\_________________________________________________
     m_flow_corr_nom_st[i] = m_flow_nom * T_out_nom_st[i-1]^0.5 / p_out_nom_st[i-1] * (1 + vigv_coeff_phi_a*Delta_alpha_int[i]^2 + vigv_coeff_phi_b*Delta_alpha_int[i]);//(1 + 0.014 * Delta_alpha_int[i]);
     m_flow_corr_st[i] = m_flow_aux * T_out_st[i-1]^0.5 / p_out_st[i-1];
     m_flow_corr_rel_st[i] = m_flow_corr_st[i]/(m_flow_corr_nom_st[i]);

     //_____________/Stage performance characteristic\____________________________________
     phi_rel_st[i] = m_flow_corr_rel_st[i]/rpm_corr_rel_st[i]; //eq. 34

     //_____________/Additional values for stage enthalpy characteristic\_________________
     if useFixedEnthalpyCharacteristic == false then
       psi_nom_st[i] = Delta_h_nom_st[i]/(Modelica.Constants.pi^2 * diameter[i]^2 * (rpm_nom/60)^2/2) * (1 + vigv_coeff_psi_a*Delta_alpha_int[i]^2 + vigv_coeff_psi_b*Delta_alpha_int[i]);//(1 + 0.014 * Delta_alpha_int[i]);//eq. 19
     else
       psi_nom_st[i] = psi_nom_fixed[i] * (1 + vigv_coeff_psi_a*Delta_alpha_int[i]^2 + vigv_coeff_psi_b*Delta_alpha_int[i]);//(1 + 0.014 * Delta_alpha_int[i]);
     end if;

     C_1_st[i] = 2/psi_nom_st[i];
     C_2_st[i] = 1 - C_1_st[i];

     //_____________/Stage enthalpy characteristic\_______________________________________
     psi_rel_st[i] = C_1_st[i] + C_2_st[i] * phi_rel_st[i]; //eq. 15

     //_____________/Stage pressure characteristic\_______________________________________
     epsilon_rel_st[i] = (1 + C_1_st[i] - C_1_st[i] * phi_rel_st[i]) * phi_rel_st[i]; //eq. 39

    //_____________/Reference point adjustion for VIGV\________________________________________________
    psi_rel_st_vigv[i] = 1 + vigv_coeff_psi_a*Delta_alpha_int[i]^2 + vigv_coeff_psi_b*Delta_alpha_int[i];//1 + 0.014 * Delta_alpha_int[i];
    epsilon_rel_st_vigv[i] =(-C_1_st[i]*(psi_rel_st_vigv[i]^2+1) + (C_1_st[i]^2+1)* psi_rel_st_vigv[i])/C_2_st[i]^2;

    der(kappa_nom_st_vigv[i]) =  1/Tau_aux*((TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(p_out_nom_st[i-1],T_out_nom_st[i-1],xi_nom,GasPointerOutletNom[i-1])/TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(p_out_nom_st[i-1],T_out_nom_st[i-1],xi_nom,GasPointerOutletNom[i-1])
                              + TILMedia.GasObjectFunctions.specificIsobaricHeatCapacity_pTxi(p_out_nom_st[i-1]*Pi_nom_st_vigv[i],T_out_nom_st[i-1]*tau_nom_st_vigv[i],xi_nom,GasPointerOutletVigv[i])/TILMedia.GasObjectFunctions.specificIsochoricHeatCapacity_pTxi(p_out_nom_st[i-1]*Pi_nom_st_vigv[i],T_out_nom_st[i-1]*tau_nom_st_vigv[i],xi_nom,GasPointerOutletVigv[i]))/2 -kappa_nom_st_vigv[i]);

    //tau_nom_st_vigv[i] = if N_VIGVstages > 0 then 1 + (Pi_nom_st_vigv[i]^((kappa_nom_st[i]-1)/kappa_nom_st[i])-1) * 1/(eta_isen_nom_st[i] *min(1 + vigv_coeff_eta_a*Delta_alpha_int[i]^3 + vigv_coeff_eta_b*Delta_alpha_int[i]^2 + vigv_coeff_eta_c*Delta_alpha_int[i],1)) else tau_nom_st[i];
    //Pi_nom_st_vigv[i] = if N_VIGVstages > 0 then (1 +  epsilon_rel_st_vigv[i] * (Pi_nom_st[i]^((kappa_nom_st[i]-1)/kappa_nom_st[i])-1))^(kappa_nom_st[i]/(kappa_nom_st[i]-1)) else Pi_nom_st[i];

     tau_nom_st_vigv[i] = if N_VIGVstages > 0 then 1 + (Pi_nom_st_vigv[i]^((kappa_nom_st_vigv[i]-1)/kappa_nom_st_vigv[i])-1) * 1/(eta_isen_nom_st[i] *min(1 + vigv_coeff_eta_a*Delta_alpha_int[i]^3 + vigv_coeff_eta_b*Delta_alpha_int[i]^2 + vigv_coeff_eta_c*Delta_alpha_int[i],1)) else tau_nom_st[i];
     Pi_nom_st_vigv[i] = if N_VIGVstages > 0 then (1 +  epsilon_rel_st_vigv[i] * (Pi_nom_st[i]^((kappa_nom_st[i]-1)/kappa_nom_st[i])-1))^(kappa_nom_st_vigv[i]/(kappa_nom_st_vigv[i]-1)) else Pi_nom_st[i];

     //_____________/Stage temperature ratio\_____________________________________________
     tau_st[i] = 1 + psi_rel_st[i] * rpm_corr_rel_st[i]^2 * (tau_nom_st_vigv[i]-1); //eq. 22
     //tau_st[i] = 1 + (Pi_st[i]^((kappa_st[i]-1)/kappa_st[i])-1) * 1/eta_isen_st[i];
     T_out_st[i] = tau_st[i] * T_out_st[i-1]; //eq. 18

     //_____________/Stage pressure ratio\________________________________________________
     //Pi_st[i] = (1 + epsilon_rel_st[i] * rpm_corr_rel_st[i]^2 * (Pi_nom_st_vigv[i]^((kappa_st[i]-1)/kappa_st[i])-1))^(kappa_st[i]/(kappa_st[i]-1)); //eq. 31
     epsilon_rel_st[i] = 1/rpm_corr_rel_st[i]^2 * ((Pi_st[i]^((kappa_st[i]-1)/kappa_st[i])-1)/(Pi_nom_st_vigv[i]^((kappa_nom_st_vigv[i]-1)/kappa_nom_st_vigv[i])-1)); //eq. 31
     p_out_st[i] = p_out_st[i-1] * Pi_st[i];

     //_____________/Isentropic stage efficiency\_________________________________________
     eta_isen_st[i] = (Pi_st[i]^((kappa_st[i]-1)/kappa_st[i])-1) / (tau_st[i] - 1);//eq. 23
     eta_isen_rel_st[i] = eta_isen_st[i]/(eta_isen_nom_st[i]*min(1 + vigv_coeff_eta_a*Delta_alpha_int[i]^3 + vigv_coeff_eta_b*Delta_alpha_int[i]^2 + vigv_coeff_eta_c*Delta_alpha_int[i],1));  //Nicht wirklich benoetigt
     //eta_isen_st[i] = eta_isen_rel_st[i]*(eta_isen_nom_st[i]*min(1 - 0.0004*Delta_alpha_int[i]^2 + 0.0007*Delta_alpha_int[i],1));
     //eta_isen_rel_st[i] = epsilon_rel_st[i]/psi_rel_st[i];

     Pi_prod_st[i] = product(Pi_st[j] for j in 1:i); //Overall pressure ratio until actual stage

     //____________/Surge line\____________________________________________________________
     if psi_nom_st[i] <= 0.9 then
       phi_surge_rel_st[i] = 1/3 +  (2 + psi_nom_st[i])/6; //eq. 66 and 40 (axial machine)
     else
       phi_surge_rel_st[i] = 2/7 + 5/7 * (2 + psi_nom_st[i])/4; //eq. from goettlich and 40 (radial machine)
     end if;

     epsilon_rel_st_max[i] = (1 + C_1_st[i] - C_1_st[i]*phi_surge_rel_st[i])*phi_surge_rel_st[i];
     Pi_st_max[i] = (1 + epsilon_rel_st_max[i]*rpm_corr_rel_st[i]^2*(Pi_nom_st_vigv[i]^((kappa_nom_st_vigv[i]-1)/kappa_nom_st_vigv[i])-1))^(kappa_st[i]/(kappa_st[i]-1));

     Y_st[i]= (p_out_st[i] - p_out_st[i-1])/(0.5*(TILMedia.GasObjectFunctions.density_pTxi(p_out_nom_st[i-1],T_out_nom_st[i-1],xi_nom,GasPointerOutletNom[i-1]) + TILMedia.GasObjectFunctions.density_pTxi(p_out_st[i],T_out_st[i],inStream(inlet.
        xi_outflow),                                                                                                    GasPointerOutlet[i])));

 end for;

    Pi = product(Pi_st);
    tau = product(tau_st);

    kappa_aux = (kappa_in + kappa_out_st[N_stages])/2.0;
    der(kappa) = 1/Tau_aux*(kappa_aux-kappa);
    //kappa=1.4;
    eta_isen = (Pi^((kappa-1)/kappa)-1) / (tau - 1);

    Pi = outlet.p/inlet.p;
    tau = flueGas_outlet.T/flueGas_inlet.T;

 end if;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////END////END////END////END////END////END////END////END////END////END////END////END////END////END////END////END////END////END////END////END////END////END////END/////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Original functions for VIGV correction from [Goettlich1984]:
////eta_rel = (1 - 0.000111 * Delta_alpha_int[i]^2);
////psi_rel = (1 + 0.014 * Delta_alpha_int[i]);
////phi_rel = psi_rel = (1 + 0.014 * Delta_alpha_int[i]);
////These functions are not in use but newer spline fittings (same source new splines)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////SURGE AND CHOKE LINE////SURGE AND CHOKE LINE////SURGE AND CHOKE LINE////SURGE AND CHOKE LINE////SURGE AND CHOKE LINE////SURGE AND CHOKE LINE//////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for i in 1:N_stages loop
    if noEvent(phi_rel_st[i] <= phi_surge_rel_st[i]) then
      count_surge[i]=1;
    else
      count_surge[i]=0;
    end if;
   if noEvent(m_flow_corr_rel_st[i] >= 1.5 and time >0) then
      count_choke[i]=1;
    else
      count_choke[i]=0;
    end if;
  end for;

  if (useBoundaryAssert) then
     assert(sum(count_surge)<=4,
            "************StackedCompressorStages: Surge line exceeded in four or more stages************");
     assert(phi_rel_st[N_stages] > phi_surge_rel_st[N_stages],
            "************StackedCompressorStages: Surge line exceeded in last stage************");
     assert(sum(count_choke)<1,
             "************StackedCompressorStages: Choke line exceeded************");
  end if;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //___________________________________________________
  P_hyd = (flueGas_outlet.h - flueGas_inlet.h) * inlet.m_flow;
  P_shaft = P_hyd*1/eta_mech;
  Y=sum(Y_st[i]);
  V_flow = inlet.m_flow / flueGas_inlet.d;

  der(inlet.m_flow) = 1/0.01*(m_flow_aux-inlet.m_flow);
  //inlet.m_flow = m_flow_aux;

  inlet.m_flow + outlet.m_flow = 0.0;
  outlet.xi_outflow = inStream(inlet.xi_outflow);
  inlet.xi_outflow = inStream(outlet.xi_outflow);
  outlet.T_outflow = flueGas_outlet.T;
  inlet.T_outflow = inStream(outlet.T_outflow);

  //______________Eye port variable definition________________________
  eye_int.m_flow = -outlet.m_flow;
  eye_int.T = flueGas_outlet.T-273.15;
  eye_int.s = flueGas_outlet.s/1e3;
  eye_int.p = flueGas_outlet.p/1e5;
  eye_int.h = flueGas_outlet.h/1e3;

  connect(eye_int,eyeOut)  annotation (Line(
      points={{40,-60},{92,-60}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(shaft, getInputsRotary.rotatoryFlange)
    annotation (Line(points={{0,100},{0,100},{0,30}}, color={0,0,0}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})),  Icon(graphics),
    Documentation(info="<html>
<p><b>Model description: </b>A multi stage compressor model based on the stage stacking method of Gasparovic able to calculate off-design behaviour according to nominal values</p>
<p><b>Contact:</b> Lasse Nielsen, TLK-Thermo GmbH</p>
<p><b>FEATURES</b> </p>
<p><ul>
<li>This model uses TILMedia</li>
<li>Stationary mass and energy balances are used</li>
<li>Modeled according to the stage stacking method of Gasparovic</li>
<li>Supports multi stage VIGV&apos;s</li>
<li>Calculation of surge and choke margins (assert triggers when crossed)</li>
</ul></p>
<p><br/><b>NOTE: </b> If the design point is not known please determine all needed values with simplier compressor models to ensure the model to run.</p>
</html>"));
end CompressorGas_L1_stageStacked;
