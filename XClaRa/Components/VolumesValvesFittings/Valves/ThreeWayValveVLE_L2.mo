within Exergy.XClaRa.Components.VolumesValvesFittings.Valves;
model ThreeWayValveVLE_L2 "A voluminous three way valve for VLE media"
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
  extends
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.ThreeWayValve_base;
  extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L2");

  outer ClaRa.SimCenter simCenter;
model Outline
  extends ClaRa.Basics.Icons.RecordIcon;
  input ClaRa.Basics.Units.Volume volume_tot "Total volume";
end Outline;

model Summary
  extends ClaRa.Basics.Icons.RecordIcon;
  Outline outline;
  ClaRa.Basics.Records.FlangeVLE           inlet;
  ClaRa.Basics.Records.FlangeVLE           outlet1;
  ClaRa.Basics.Records.FlangeVLE           outlet2;
  ClaRa.Basics.Records.FluidVLE_L2           fluid;
end Summary;

  replaceable model PressureLoss =
      Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.IdealSymetric_TWV
    constrainedby
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.TWV_L2
    "Pressure loss model"                                                                   annotation(choicesAllMatching);

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid   medium=simCenter.fluid1
    "Medium in the component"  annotation(Dialog(group="Fundamental Definitions"));

  parameter Boolean splitRatio_input=false
    "= true, if split ratio is defined by input" annotation(Dialog(group="Fundamental Definitions"));
  parameter Real splitRatio_fixed = 0.5 annotation(Dialog(enable=not splitRatio_input, group= "Fundamental Definitions"));

  parameter Boolean useHomotopy=simCenter.useHomotopy
    "True, if homotopy method is used during initialisation"  annotation(Dialog(tab="Initialisation"));
   parameter SI.Volume volume(min=1e-6)=0.1 "System Volume"                               annotation(Dialog(tab="General", group="Geometry"));
  parameter SI.MassFlowRate m_flowOut_nom[2]= {10, 10}
    "Nominal mass flow rates at outlet"  annotation(Dialog(tab="General", group="Nominal Values"));
  parameter SI.Pressure p_nom=1e5 "Nominal pressure"                    annotation(Dialog(group="Nominal Values"));
  parameter SI.EnthalpyMassSpecific h_nom=1e5 "Nominal specific enthalpy"                      annotation(Dialog(group="Nominal Values"));

  parameter SI.EnthalpyMassSpecific h_start= 1e5
    "Start value of sytsem specific enthalpy"
                                             annotation(Dialog(tab="Initialisation"));
  parameter SI.Pressure p_start= 1e5 "Start value of sytsem pressure"               annotation(Dialog(tab="Initialisation"));
  parameter ClaRa.Basics.Choices.Init initType=ClaRa.Basics.Choices.Init.noInit
    "Type of initialisation" annotation(Dialog(tab="Initialisation"), choicesAllMatching);
  parameter Boolean showExpertSummary=simCenter.showExpertSummary
    "|Summary and Visualisation||True, if expert summary should be applied";
  parameter Boolean showData=true
    "|Summary and Visualisation||True, if a data port containing p,T,h,s,m_flow shall be shown, else false";
  parameter Boolean preciseTwoPhase = true
    "|Expert Settings||True, if two-phase transients should be capured precisely";

protected
    parameter SI.DensityMassSpecific rho_nom= TILMedia.VLEFluidFunctions.density_phxi(medium, p_nom, h_nom)
    "Nominal density";
    SI.Power Hdrhodt =  if preciseTwoPhase then h*volume*drhodt else 0
    "h*V*drhodt";
public
  Real splitRatio;
  SI.EnthalpyFlowRate H_flow_in;
  SI.EnthalpyFlowRate H_flow_out[2];
  SI.EnthalpyMassSpecific h(start=h_start);
  SI.Mass mass "Total system mass";
  Real drhodt;//(unit="kg/(m3s)");
  SI.Pressure p(start=p_start, stateSelect=StateSelect.prefer)
    "System pressure";

public
   Summary summary(outline(volume_tot = volume),
                   inlet(showExpertSummary = showExpertSummary,m_flow=inlet.m_flow,  T=fluidIn.T, p=inlet.p, h=fluidIn.h,s=fluidIn.s, steamQuality=fluidIn.q, H_flow=fluidIn.h*inlet.m_flow, rho=fluidIn.d),
                   fluid(showExpertSummary = showExpertSummary, mass=mass, p=p, h=h, T=bulk.T,s=bulk.s, steamQuality=bulk.q, H=h*mass, rho=bulk.d, T_sat=bulk.VLE.T_l, h_dew=bulk.VLE.h_v, h_bub=bulk.VLE.h_l),
                   outlet1(showExpertSummary = showExpertSummary,m_flow = -outlet1.m_flow, T=fluidOut1.T, p=outlet1.p, h=fluidOut1.h, s=fluidOut1.s, steamQuality=fluidOut1.q, H_flow=-fluidOut1.h*outlet1.m_flow, rho=fluidOut1.d),
                   outlet2(showExpertSummary = showExpertSummary,m_flow = -outlet2.m_flow, T=fluidOut2.T, p=outlet2.p, h=fluidOut2.h, s=fluidOut2.s, steamQuality=fluidOut2.q, H_flow=-fluidOut2.h*outlet2.m_flow, rho=fluidOut2.d))
    annotation (Placement(transformation(extent={{-60,-102},{-40,-82}})));

protected
  TILMedia.VLEFluid_ph bulk(each vleFluidType = medium, p = p,h=h) annotation (Placement(transformation(extent={{-10,-12},
            {10,8}},rotation=0)));

  inner
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.ICom_TWV
    iCom(
    p_in=fluidIn.p,
    p_out1=fluidOut1.p,
    p_out2=fluidOut2.p,
    rho_out1=fluidOut1.d,
    rho_out2=fluidOut2.d,
    opening_=splitRatio,
    m_flow_in=inlet.m_flow,
    rho_in=fluidIn.d) annotation (Placement(transformation(extent={{-80,
            -100},{-60,-80}})));
PressureLoss pressureLoss annotation (Placement(transformation(extent={{-40,-40},{-20,-20}})));
equation
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Asserts ~~~~~~~~~~~~~~~~~~~
  assert(volume>0, "The system volume must be greater than zero!");
//~~~~~~~~~~~~~~~~~~~~~~~~~~~

// System definition ~~~~~~~~
   mass= if useHomotopy then volume*homotopy(bulk.d,rho_nom) else volume*bulk.d;

   drhodt*volume=inlet.m_flow + outlet1.m_flow + outlet2.m_flow "Mass balance";
   drhodt=der(p)*bulk.drhodp_hxi
                             + der(h)*bulk.drhodh_pxi;
                                                   //calculating drhodt from state variables

   der(h) = 1/mass*(sum(H_flow_out) + H_flow_in  + volume*der(p) -Hdrhodt)
    "Energy balance, decoupled from the mass balance to avoid heavy mass fluctuations during phase change or flow reversal. The term '-h*V*drhodt' is ommited";
//~~~~~~~~~~~~~~~~~~~~~~~~~
// Boundary conditions ~~~~

    H_flow_out[1]=if useHomotopy then homotopy(actualStream(outlet1.h_outflow)*outlet1.m_flow, -h*m_flowOut_nom[1]) else actualStream(outlet1.h_outflow)*outlet1.m_flow;
    H_flow_out[2]=if useHomotopy then homotopy(actualStream(outlet2.h_outflow)*outlet2.m_flow, -h*m_flowOut_nom[2]) else actualStream(outlet2.h_outflow)*outlet2.m_flow;
    outlet1.m_flow=-pressureLoss.m_flow_1;
    outlet1.h_outflow=h;
    outlet2.m_flow=-pressureLoss.m_flow_2;
    outlet2.h_outflow=h;

    H_flow_in= if useHomotopy then homotopy(actualStream(inlet.h_outflow)*inlet.m_flow, inStream(inlet.h_outflow)*sum(m_flowOut_nom)) else actualStream(inlet.h_outflow)*inlet.m_flow;
    inlet.p=p;
    inlet.h_outflow=h;

initial equation
  if initType == ClaRa.Basics.Choices.Init.steadyState then
    der(h)=0;
    der(p)=0;
  elseif initType == ClaRa.Basics.Choices.Init.steadyPressure then
    der(p)=0;
  elseif initType == ClaRa.Basics.Choices.Init.steadyEnthalpy then
    der(h)=0;
  end if;

equation

  annotation (Diagram(graphics));
end ThreeWayValveVLE_L2;
