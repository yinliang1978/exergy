within Exergy.XClaRa.Components.VolumesValvesFittings.Fittings;
model SplitVLE_L2_flex
  "A voluminous split for an arbitrary number of inputs NOT CAPABLE FOR PHASE-CHANGE"
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

 extends ClaRa.Basics.Interfaces.DataInterfaceVector(N_sets=N_ports_out);
  extends ClaRa.Basics.Icons.Adapter5_fw;

  extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L2");

  outer ClaRa.SimCenter simCenter;

model Outline
  extends ClaRa.Basics.Icons.RecordIcon;
  input ClaRa.Basics.Units.Volume volume_tot "Total volume";
end Outline;

model Summary
  parameter Integer N_ports_out;
  extends ClaRa.Basics.Icons.RecordIcon;
  Outline outline;
  ClaRa.Basics.Records.FlangeVLE           inlet;
  ClaRa.Basics.Records.FlangeVLE           outlet[N_ports_out];
  ClaRa.Basics.Records.FluidVLE_L2           fluid;
end Summary;

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid   medium=simCenter.fluid1
    "Medium in the component"  annotation(Dialog(group="Fundamental Definitions"));
  parameter Integer N_ports_out(min=1)=1 "Number of outlet  ports"
    annotation(Evaluate=true, Dialog(tab="General",group="Fundamental Definitions"));//connectorSizing=true,
  parameter Boolean useHomotopy=simCenter.useHomotopy
    "True, if homotopy method is used during initialisation"  annotation(Dialog(tab="Initialisation"));
   parameter Modelica.SIunits.Volume V(min=1e-6)=0.1 "System Volume"                               annotation(Dialog(tab="General", group="Geometry"));
  parameter Modelica.SIunits.MassFlowRate m_flow_out_nom[N_ports_out]= {10}
    "Nominal mass flow rates at inlet"  annotation(Dialog(tab="General", group="Nominal Values"));
  parameter Modelica.SIunits.Pressure p_nom=1e5 "Nominal pressure"                    annotation(Dialog(group="Nominal Values"));
  parameter Modelica.SIunits.SpecificEnthalpy h_nom=1e5
    "Nominal specific enthalpy"                                                                annotation(Dialog(group="Nominal Values"));

  parameter Modelica.SIunits.SpecificEnthalpy h_start= 1e5
    "Start value of sytsem specific enthalpy"
                                             annotation(Dialog(tab="Initialisation"));
  parameter Modelica.SIunits.Pressure p_start= 1e5
    "Start value of sytsem pressure"                                                annotation(Dialog(tab="Initialisation"));
  parameter ClaRa.Basics.Choices.Init initType=ClaRa.Basics.Choices.Init.noInit
    "Type of initialisation" annotation(Dialog(tab="Initialisation"), choicesAllMatching);
  final parameter Boolean  showSummary=false
    "True, if a summary shall be shown, else false"                                   annotation(Dialog(tab="Summary and Visualisation"));
    parameter Boolean showExpertSummary=false
    "True, if a summary shall be shown, else false"                                   annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean showData=true
    "True, if a data port containing p,T,h,s,m_flow shall be shown, else false"
                                                                                annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean preciseTwoPhase = true
    "|Expert Stettings||True, if two-phase transients should be capured precisely";
protected
    parameter Modelica.SIunits.Density rho_nom= TILMedia.VLEFluidFunctions.density_phxi(medium, p_nom, h_nom)
    "Nominal density";
    Modelica.SIunits.Power Hdrhodt =  if preciseTwoPhase then h*V*drhodt else 0
    "h*V*drhodt";
public
  Modelica.SIunits.EnthalpyFlowRate H_flow_in;
  Modelica.SIunits.EnthalpyFlowRate H_flow_out[N_ports_out];
  Modelica.SIunits.SpecificEnthalpy h(start=h_start);
  Modelica.SIunits.Mass mass "Total system mass";
  Real drhodt;//(unit="kg/(m3s)");
  Modelica.SIunits.Pressure p(start=p_start, stateSelect=StateSelect.prefer)
    "System pressure";

   Summary summary(N_ports_out=N_ports_out,outline(volume_tot = V),
                   inlet(showExpertSummary = showExpertSummary,m_flow=inlet.m_flow,  T=fluidIn.T, p=inlet.p, h=fluidIn.h,s=fluidIn.s, steamQuality=fluidIn.q, H_flow=fluidIn.h .* inlet.m_flow, rho=fluidIn.d),
                   fluid(showExpertSummary = showExpertSummary, mass=mass, p=p, h=h, T=bulk.T,s=bulk.s, steamQuality=bulk.q, H=h*mass, rho=bulk.d, T_sat=bulk.VLE.T_l, h_dew=bulk.VLE.h_v, h_bub=bulk.VLE.h_l),
                   outlet(each showExpertSummary = showExpertSummary,m_flow = -outlet.m_flow, T=fluidOut.T, p=outlet.p, h=fluidOut.h, s=fluidOut.s, steamQuality=fluidOut.q, H_flow=-fluidOut.h .* outlet.m_flow, rho=fluidOut.d))
    annotation (Placement(transformation(extent={{-60,-102},{-40,-82}})));

  ClaRa.Basics.Interfaces.FluidPortIn inlet(Medium=medium) "Inlet port"
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
  ClaRa.Basics.Interfaces.FluidPortOut outlet[N_ports_out](each Medium=medium)
    "Outlet port"
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));
protected
TILMedia.VLEFluid_ph bulk(vleFluidType =    medium,    p = p, h=h) annotation (Placement(transformation(extent={{70,-10},{90,10}}, rotation=0)));
TILMedia.VLEFluid_ph fluidIn(vleFluidType =    medium,    p = inlet.p, h=actualStream(inlet.h_outflow)) annotation (Placement(transformation(extent={{70,-10},{90,10}}, rotation=0)));
TILMedia.VLEFluid_ph fluidOut[N_ports_out](each vleFluidType =    medium,    p = outlet.p, h=actualStream(outlet.h_outflow)) annotation (Placement(transformation(extent={{70,-10},{90,10}}, rotation=0)));
equation
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Asserts ~~~~~~~~~~~~~~~~~~~
  assert(V>0, "The system volume must be greater than zero!");
//~~~~~~~~~~~~~~~~~~~~~~~~~~~
// System definition ~~~~~~~~
   mass= if useHomotopy then V*homotopy(bulk.d,rho_nom) else V*bulk.d;

   drhodt*V=sum(inlet.m_flow) + sum(outlet.m_flow) "Mass balance";
   drhodt=der(p)*bulk.drhodp_hxi
                             + der(h)*bulk.drhodh_pxi;
                                                   //calculating drhodt from state variables

   der(h) = 1/mass*(H_flow_in + sum(H_flow_out)  + V*der(p) -Hdrhodt)
    "Energy balance, decoupled from the mass balance to avoid heavy mass fluctuations during phase change or flow reversal. The term '-h*V*drhodt' is ommited";
//~~~~~~~~~~~~~~~~~~~~~~~~~
// Boundary conditions ~~~~
  for i in 1:N_ports_out loop
    outlet[i].h_outflow=h;
    H_flow_out[i]=if useHomotopy then homotopy(actualStream(outlet[i].h_outflow)*outlet[i].m_flow, -h*m_flow_out_nom[i]) else actualStream(outlet[i].h_outflow)*outlet[i].m_flow;
    outlet[i].p=p;
  end for;

    H_flow_in= if useHomotopy then homotopy(actualStream(inlet.h_outflow)*inlet.m_flow, inStream(inlet.h_outflow)*sum(m_flow_out_nom)) else actualStream(inlet.h_outflow)*inlet.m_flow;
    inlet.p=p;
    inlet.h_outflow=h;
  for i in 1:N_ports_out loop
    eye[i].m_flow=-outlet[i].m_flow;
    eye[i].T= bulk.T-273.15;
    eye[i].s=bulk.s/1e3;
    eye[i].p=bulk.p/1e5;
    eye[i].h=h/1e3;
  end for;

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

  annotation (Icon(graphics), Diagram(graphics));
end SplitVLE_L2_flex;
