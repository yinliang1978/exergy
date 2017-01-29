within Exergy.XClaRa.Components.TurboMachines.Compressors;
model CompressorGas_L1_affinity
  "A gas compressor or fan based on affinity laws"
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
import SI = ClaRa.Basics.Units;
import Modelica.Constants.pi;

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
    input SI.Power P_shaft "Hydraulic power";
    input Real Pi "Pressure ratio";
    input SI.PressureDifference Delta_p "Pressure difference";
    input SI.RPM rpm "Rotational speed";
    input Real eta "Hydraulic efficiency";
    input Real eta_mech "Mechanic efficiency";
  end Outline;

  model Summary
    extends ClaRa.Basics.Icons.RecordIcon;
    Outline outline;
    ClaRa.Basics.Records.FlangeGas  inlet;
    ClaRa.Basics.Records.FlangeGas  outlet;
  end Summary;

  inner parameter TILMedia.GasTypes.BaseGas    medium = simCenter.flueGasModel;

  final parameter Boolean allow_reverseFlow = false;

  ClaRa.Basics.Interfaces.GasPortIn inlet(Medium=medium, m_flow(min=if
          allow_reverseFlow then -Modelica.Constants.inf else 1e-5))
    "inlet flow"
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

  TILMedia.Gas_ph flueGas_outlet( gasType = medium, h = h_out, p = outlet.p,  xi = flueGas_inlet.xi)
    annotation (Placement(transformation(extent={{70,-12},{90,8}})));

  //__________________________/ Parameters \_____________________________
  parameter SI.RPM rpm_nom "|Characteristic field|Nomial rotational speed";
  parameter SI.VolumeFlowRate V_flow_max
    "|Characteristic field|Maximum volume flow rate at nominal speed";
  parameter SI.VolumeFlowRate V_flow_min= 0
    "|Characteristic field|V_flow(Delta_p_max, rpm_nom)";
  parameter SI.Pressure Delta_p_max
    "|Characteristic field|Maximum pressure difference at nominal speed";
  parameter Real exp_hyd= 0.5 "|Characteristic field|Exponent for affinity law";

  parameter Real eta = 0.85 "isentropic efficiency";
  parameter Real eta_mech = 0.99 "mechanical efficiency";
  parameter SI.Pressure Delta_p_eps= 100
    "|Expert Settings|Numerical Robustness|Small pressure difference for linearisation around zero mass flow";
  parameter Modelica.SIunits.Inertia J "Moment of Inertia" annotation(Dialog(group="Time Response Definitions", enable= not steadyStateTorque));
    parameter Boolean useMechanicalPort=false
    "|Fundamental Definitions|True, if a mechenical flange should be used";
  parameter Boolean steadyStateTorque=false
    "|Fundamental Definitions|True, if steady state mechanical momentum shall be used";
  parameter SI.RPM rpm_fixed = 60 "Constant rotational speed of pump" annotation (Dialog( group = "Fundamental Definitions", enable = not useMechanicalPort));
  parameter SI.Time Tau_aux=0.1 "Time constant of auxilliary kappa states"  annotation(Dialog(tab = "Advanced"));
  parameter Real kappa_initial = 1.3 "Initial value for kappas" annotation(Dialog(tab = "Advanced"));

  //________________________/ Variables \___________________________________
  SI.Pressure Delta_p(final start=100) "pressure increase";
  SI.Power P_hyd "Hydraulic power";
  SI.VolumeFlowRate V_flow;
  SI.VolumeFlowRate V_flow_max_aff;
  SI.Pressure Delta_p_max_aff;
  Modelica.SIunits.AngularAcceleration a "Angular acceleration of the shaft";
  SI.Power P_shaft "Mechanical power at shaft";
  SI.RPM rpm "Rotational speed";
  Modelica.SIunits.Torque tau_fluid "Fluid torque";
  SI.EnthalpyMassSpecific Delta_h;
  Real kappa;
  Real kappaA;
  Real kappaB;

protected
  SI.EnthalpyMassSpecific h_out;
  Real kappaB_aux;
  Real kappaA_aux;

public
  ClaRa.Basics.Interfaces.EyeOut eyeOut annotation (Placement(
        transformation(extent={{72,-78},{112,-42}}),
        iconTransformation(extent={{92,-70},{112,-50}})));
protected
  ClaRa.Basics.Interfaces.EyeIn eye_int annotation (Placement(
        transformation(extent={{48,-68},{32,-52}}),
        iconTransformation(extent={{90,-84},{84,-78}})));

public
  inner Summary summary(outline(
   V_flow = V_flow,
   P_hyd = P_hyd,
   P_shaft = P_shaft,
   Pi = outlet.p/inlet.p,
   Delta_p = outlet.p - inlet.p,
   rpm = rpm,
   eta = eta,
   eta_mech = eta_mech),
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
          xi= outlet.xi_outflow,
          H_flow = -outlet.m_flow* flueGas_outlet.h)) annotation (Placement(transformation(extent={{-100,
            -114},{-80,-94}})));
initial equation
  kappaB = kappa_initial;
  kappaA = kappa_initial;

equation

//____________________ Mechanics ___________________________
  if useMechanicalPort then
    der(getInputsRotary.rotatoryFlange.phi) = (2*pi*rpm/60);
    J*a*rpm = - tau_fluid*2*pi*rpm/60 + getInputsRotary.rotatoryFlange.tau*2*pi*rpm/60
      "Mechanical momentum balance";
  else
    rpm = rpm_fixed;
    getInputsRotary.rotatoryFlange.phi = 0.0;
  end if;

  if (steadyStateTorque) then
    a = 0;
  else
    a = (2*pi/60)^2*der(rpm);
  end if;
  tau_fluid = if noEvent(2*pi*rpm/60<1e-8) then 0 else P_shaft/(2*pi*rpm/60);

  //______________Affinity laws_______________________
  V_flow_max_aff = V_flow_max*rpm/rpm_nom;
  Delta_p_max_aff = Delta_p_max*(rpm/rpm_nom)^2;

  //______________Compressor_characteristic___________
  V_flow =  FluidDissipation.Utilities.Functions.General.SmoothPower((Delta_p_max_aff - Delta_p)/Delta_p_max, Delta_p_eps/Delta_p_max, exp_hyd)*(V_flow_max-V_flow_min) + V_flow_min;

  P_hyd = Delta_h*inlet.m_flow;

  P_shaft = P_hyd*1/eta*1/eta_mech;

  Delta_p = outlet.p - inlet.p;
  h_out = flueGas_inlet.h + Delta_h;

  V_flow = inlet.m_flow / flueGas_inlet.d;
  flueGas_inlet.cv * kappaA_aux = flueGas_inlet.cp;
  flueGas_outlet.cv * kappaB_aux = flueGas_outlet.cp;

  der(kappaB) = 1/Tau_aux*(kappaB_aux-kappaB); // auxilluary state
  der(kappaA) = 1/Tau_aux*(kappaA_aux-kappaA);

  if ( kappaA + kappaB) > 0.3 then
    kappa = ( kappaA + kappaB)/2.0;
  else
    kappa = 0.2;
  end if;

  eta * Delta_h =  kappa/(kappa - 1.0) * flueGas_inlet.p/flueGas_inlet.d * ((flueGas_outlet.p/flueGas_inlet.p)^((kappa -1.0)/kappa) - 1.0);

  inlet.m_flow + outlet.m_flow = 0.0;
  outlet.xi_outflow = inStream(inlet.xi_outflow);
  inlet.xi_outflow = inStream(outlet.xi_outflow);
  outlet.T_outflow = flueGas_outlet.T;
  inlet.T_outflow = actualStream(outlet.T_outflow);

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
    annotation (Line(points={{0,100},{0,30}}, color={0,0,0}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})),  Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}),
                                      graphics));
end CompressorGas_L1_affinity;
