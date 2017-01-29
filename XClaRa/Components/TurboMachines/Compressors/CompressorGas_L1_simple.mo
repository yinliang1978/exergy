within Exergy.XClaRa.Components.TurboMachines.Compressors;
model CompressorGas_L1_simple "Simple compressor or fan for gas"
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
//parameter Boolean allow_reverseFlow = true annotation(Evaluate=true, Dialog(tab="Advanced"));
extends ClaRa.Basics.Icons.Compressor;
parameter Boolean contributeToCycleSummary = simCenter.contributeToCycleSummary
    "True if component shall contribute to automatic efficiency calculation"                annotation(Dialog(tab="Summary and Visualisation"));
  ClaRa.Basics.Interfaces.Connected2SimCenter connected2SimCenter(
    powerIn=0,
    powerOut=-P_hyd,
    powerAux=-P_hyd + P_shaft) if                                                                                                  contributeToCycleSummary;

  model Outline
    extends ClaRa.Basics.Icons.RecordIcon;
    input SI.VolumeFlowRate V_flow "Volume flow rate";
    input SI.Power P_hyd "Hydraulic power";
    input SI.Power P_shaft "Hydraulic power";
    input Real Pi "Pressure ratio";
    input SI.PressureDifference Delta_p "Pressure difference";
    input Real eta "Hydraulic efficiency";
  end Outline;

  model Summary
    extends ClaRa.Basics.Icons.RecordIcon;
    Outline outline;
    ClaRa.Basics.Records.FlangeGas  inlet;
    ClaRa.Basics.Records.FlangeGas  outlet;
  end Summary;

final parameter Boolean allow_reverseFlow = false;
  ClaRa.Basics.Interfaces.GasPortIn inlet(Medium=medium, m_flow(min=if
          allow_reverseFlow then -Modelica.Constants.inf else 1e-5))
    "inlet flow"
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
  ClaRa.Basics.Interfaces.GasPortOut outlet(Medium=medium, m_flow(max=if
          allow_reverseFlow then Modelica.Constants.inf else -1e-5))
    "outlet flow"
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));

inner parameter TILMedia.GasTypes.BaseGas    medium = simCenter.flueGasModel;
 //parameter ClaRa.Basics.Units.MassFraction[medium.nc] mixingRatioInitial = medium.defaultMixingRatio
 //   "Initial value for mixing ratio" annotation(Dialog(group="Initial Values"));
  TILMedia.Gas_pT flueGas_inlet( p = inlet.p, T = inStream(inlet.T_outflow), xi = inStream(inlet.xi_outflow), gasType = medium)
    annotation (Placement(transformation(extent={{-90,-12},{-70,8}})));

  TILMedia.Gas_ph flueGas_outlet( gasType = medium, h = hOut, p = outlet.p,  xi = flueGas_inlet.xi)
    annotation (Placement(transformation(extent={{70,-12},{90,8}})));

Real kappa;
Real deltah;
Real hOut;
Real kappaA;
Real kappaB;
Real kappaB_aux;
Real kappaA_aux;

public
 SI.Pressure Delta_p(final start=100) "pressure increase";
 SI.Power P_hyd "Hydraulic power";
 SI.Power P_shaft "Drive power";
 SI.VolumeFlowRate V_flow;

parameter Real eta = 0.85 "isentropic efficiency";

//__________________________/ Inputs \_____________________________

  parameter
    Exergy.XClaRa.Components.TurboMachines.Compressors.Fundamentals.PresetVariableType
    presetVariableType="V_flow" "Specifies which variable is preset"
    annotation (Dialog(group="General Settings"));

parameter Boolean use_P_shaftInput=false "= true, if P_shaft defined by input"
    annotation(Dialog(enable=(presetVariableType=="P_shaft"), group="Mechanical shaft power"));

  parameter SI.Power P_shaft_fixed=5e3 "Fixed value for mechanical shaft power"
                                             annotation (Dialog(
        enable=(not use_P_shaftInput and presetVariableType ==
          "P_shaft"), group="Mechanical shaft power"));

   parameter Boolean use_Delta_p_input=false
    "= true, if Delta_p defined by input"
    annotation(Dialog(enable=(presetVariableType=="dp"), group="Pressure Increase"));
  parameter SI.Pressure Delta_p_fixed=0.1e5 "Fixed value for pressure increase"
                                        annotation (Dialog(enable=(
          not use_Delta_p_input and presetVariableType == "dp"),
        group="Pressure Increase"));

  parameter Boolean m_flowInput=false "= true, if m_flow defined by input"
    annotation(Dialog(enable=(presetVariableType=="m_flow"), group="Mass Flow Rate"));
  parameter SI.MassFlowRate m_flow_fixed=0.5
    "Fixed value for gas mass flow rate" annotation (Dialog(enable=(
          not m_flowInput and presetVariableType == "m_flow"), group=
          "Mass Flow Rate"));

  parameter Boolean V_flowInput=false "= true, if V_flow defined by input"
    annotation(Dialog(enable=(presetVariableType=="V_flow"), group="Volume Flow Rate"));
  parameter SI.VolumeFlowRate V_flow_fixed=0.5e-3
    "Fixed value for gas volume flow rate" annotation (Dialog(enable=
          (not V_flowInput and presetVariableType == "V_flow"), group=
         "Volume Flow Rate"));

parameter ClaRa.Basics.Units.Time Tau_aux=0.1
    "Time constant of auxilliary kappa states"  annotation(Dialog(tab = "Advanced"));

 parameter Real kappa_initial = 1.3 "Initial value for kappas" annotation(Dialog(tab = "Advanced"));

  Exergy.XClaRa.Components.TurboMachines.Compressors.Fundamentals.GetInputsHydraulic
    getInputs annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={0,40})));
  Modelica.Blocks.Interfaces.RealInput dp_in if
    use_Delta_p_input "Prescribed pressure increase"
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},  rotation=270,
        origin={80,110}),
                iconTransformation(extent={{-10,-10},{10,10}},
        rotation=270,
        origin={80,110})));
  Modelica.Blocks.Interfaces.RealInput P_shaft_in if
     use_P_shaftInput "Prescribed pressure increase"
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},  rotation=270,
        origin={34,110}),
                iconTransformation(extent={{-10,-10},{10,10}},
        rotation=270,
        origin={30,110})));
  Modelica.Blocks.Interfaces.RealInput V_flow_in if V_flowInput
    "Prescribed volume flow rate"
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},   rotation=270,
        origin={-32,110}),
               iconTransformation(extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-30,110})));
  Modelica.Blocks.Interfaces.RealInput m_flow_in if
                                                  m_flowInput
    "Prescribed mass flow rate"
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=270,
        origin={-80,110}),
                       iconTransformation(extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-80,110})));

protected
  Modelica.Blocks.Sources.Constant Delta_p_in_(k=0) if not use_Delta_p_input;
  Modelica.Blocks.Sources.Constant m_flow_in_(k=0) if not m_flowInput;
  Modelica.Blocks.Sources.Constant V_flow_in_(k=0) if not V_flowInput;
  Modelica.Blocks.Sources.Constant P_shaft_in_(k=0) if not use_P_shaftInput;

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
   eta = eta),
    inlet(m_flow = inlet.m_flow,
          T = inStream(inlet.T_outflow),
          p = inlet.p,
          h = flueGas_inlet.h,
          xi=inStream(inlet.xi_outflow),
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

  //____________________ Boundary equations _________________

  if presetVariableType == "dp" then
    if use_Delta_p_input then
      Delta_p = getInputs.dp_in;
    else
      Delta_p = Delta_p_fixed;
    end if;
  elseif presetVariableType == "m_flow" then
    if m_flowInput then
      inlet.m_flow = getInputs.m_flow_in;
    else
      inlet.m_flow = m_flow_fixed;
    end if;
  elseif presetVariableType == "V_flow" then
    if V_flowInput then
      V_flow = getInputs.V_flow_in;
    else
      V_flow = V_flow_fixed;
    end if;
  else
    if use_P_shaftInput then
      P_shaft = getInputs.P_shaft_in;
    else
      P_shaft = P_shaft_fixed;
    end if;
  end if;

Delta_p = inlet.p - outlet.p;
hOut = flueGas_inlet.h + deltah;
P_hyd = deltah*inlet.m_flow;
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

eta * deltah =  kappa/(kappa - 1.0) * flueGas_inlet.p/flueGas_inlet.d * ((flueGas_outlet.p/flueGas_inlet.p)^((kappa -1.0)/kappa) - 1.0);

inlet.m_flow + outlet.m_flow = 0.0;
  outlet.xi_outflow
                  = inStream(inlet.xi_outflow);
  inlet.xi_outflow
                 = inStream(outlet.xi_outflow);
outlet.T_outflow = flueGas_outlet.T;
inlet.T_outflow = actualStream(outlet.T_outflow);
inlet.m_flow * (deltah+1e-40) = P_shaft;

  //______________Eye port variable definition________________________
  eye_int.m_flow = -outlet.m_flow;
  eye_int.T = flueGas_outlet.T-273.15;
  eye_int.s = flueGas_outlet.s/1e3;
  eye_int.p = flueGas_outlet.p/1e5;
  eye_int.h = flueGas_outlet.h/1e3;

connect(Delta_p_in_.y, getInputs.dp_in);
  connect(V_flow_in_.y, getInputs.V_flow_in);
  connect(m_flow_in_.y, getInputs.m_flow_in);
  connect(P_shaft_in_.y, getInputs.P_shaft_in);
  connect(m_flow_in, getInputs.m_flow_in) annotation (Line(
      points={{-80,110},{-80,76},{-8,76},{-8,52}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(V_flow_in, getInputs.V_flow_in) annotation (Line(
      points={{-32,110},{-32,88},{-2,88},{-2,52}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(P_shaft_in, getInputs.P_shaft_in) annotation (Line(
      points={{34,110},{34,88},{2,88},{2,52}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(dp_in, getInputs.dp_in) annotation (Line(
      points={{80,110},{80,76},{8,76},{8,52}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(eye_int,eyeOut)  annotation (Line(
      points={{40,-60},{92,-60}},
      color={190,190,190},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                      graphics), Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}),
                                      graphics));
end CompressorGas_L1_simple;
