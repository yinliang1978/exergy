within Exergy.XClaRa.Components.VolumesValvesFittings.Fittings;
model SplitGas_L2_flex "Adiabatic junction volume"
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

extends ClaRa.Basics.Icons.Adapter5_fw;
extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L2");
  outer ClaRa.SimCenter simCenter;

 model Gas
  extends ClaRa.Basics.Icons.RecordIcon;
  input ClaRa.Basics.Units.Mass m "Mass flow rate"
                                               annotation(Dialog);
  input ClaRa.Basics.Units.Temperature T "Temperature"
                                                   annotation(Dialog);
  input ClaRa.Basics.Units.Pressure p "Pressure"
                                             annotation(Dialog);
  input ClaRa.Basics.Units.EnthalpyMassSpecific h "Specific enthalpy"
                                                                  annotation(Dialog);
  input ClaRa.Basics.Units.Enthalpy H "Specific enthalpy"
                                                      annotation(Dialog);
  input ClaRa.Basics.Units.DensityMassSpecific rho "Specific enthalpy"
                                                                   annotation(Dialog);
 end Gas;

 inner model Summary
   parameter Integer N_ports_out;
   extends ClaRa.Basics.Icons.RecordIcon;
   Gas gas;
   ClaRa.Basics.Records.FlangeGas inlet;
   ClaRa.Basics.Records.FlangeGas outlet[N_ports_out];
 end Summary;

inner parameter ClaRa.Basics.Choices.Init initType=ClaRa.Basics.Choices.Init.noInit
    "Type of initialisation" annotation(Dialog(tab="Initialisation"), choicesAllMatching);

// ***************************** defintion of medium used in cell *************************************************
inner parameter TILMedia.GasTypes.BaseGas medium = simCenter.flueGasModel
    "Medium to be used in tubes"  annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));

  ClaRa.Basics.Interfaces.GasPortIn inlet(Medium=medium, m_flow)
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));
  ClaRa.Basics.Interfaces.GasPortIn outlet[N_ports_out](each Medium=medium, m_flow)
    annotation (Placement(transformation(extent={{90,-10},{110,10}}), iconTransformation(extent={{90,-10},{110,10}})));
 parameter Integer  N_ports_out(min=1)=1 "Number of inlet  ports"
    annotation(Evaluate=true, Dialog(tab="General",group="Fundamental Definitions"));//connectorSizing=true,

parameter ClaRa.Basics.Units.Volume volume=1 annotation(Dialog(tab="General",group="Geometry"));

protected
  TILMedia.Gas_pT     gasInlet(gasType = medium, p=inlet.p, T=noEvent(actualStream(inlet.T_outflow)), xi=noEvent(actualStream(inlet.xi_outflow)))
    annotation (Placement(transformation(extent={{-80,-12},{-60,8}})));
protected
  TILMedia.Gas_pT gasOutlet[N_ports_out](each
    gasType=medium,
    p=outlet.p,
    T=noEvent(actualStream(outlet.T_outflow)),
    xi=noEvent(actualStream(outlet.xi_outflow)))
    annotation (Placement(transformation(extent={{60,-14},{80,6}})));
protected
  inner TILMedia.Gas_ph     bulk(
    computeTransportProperties=false,
    gasType = medium,p=p,h=h,xi=xi,
    stateSelectPreferForInputs=true)
    annotation (Placement(transformation(extent={{-10,-12},{10,8}})));
  /****************** Nominal values *******************/
public
  parameter Modelica.SIunits.MassFlowRate m_flow_out_nom[N_ports_out]= {10}
    "Nominal mass flow rates at inlet"  annotation(Dialog(tab="General", group="Nominal Values"));
  parameter Modelica.SIunits.Pressure p_nom=1e5 "Nominal pressure"                    annotation(Dialog(group="Nominal Values"));
  parameter Modelica.SIunits.Temperature T_nom=293.15
    "Nominal specific enthalpy"                                                                annotation(Dialog(group="Nominal Values"));
  parameter ClaRa.Basics.Units.MassFraction xi_nom[medium.nc - 1]={0,0,0,0,0.76,0.23,0,0,0}  annotation(Dialog(group="Nominal Values"));

  final parameter Modelica.SIunits.Density rho_nom= TILMedia.GasFunctions.density_pTxi(medium, p_nom, T_nom, xi_nom)
    "Nominal density";
  /****************** Initial values *******************/
public
    parameter Boolean useHomotopy=simCenter.useHomotopy
    "True, if homotopy method is used during initialisation"  annotation(Dialog(tab="Initialisation"));
  parameter ClaRa.Basics.Units.Pressure p_start=1.013e5
    "Initial value for air pressure"
    annotation(Dialog(tab="Initialisation"));

  parameter ClaRa.Basics.Units.Temperature T_start=298.15
    "Initial value for air temperature"
    annotation(Dialog(tab="Initialisation"));

  parameter ClaRa.Basics.Units.MassFraction[medium.nc - 1]
                                                         mixingRatio_initial=zeros(medium.nc-1)
    "Initial value for mixing ratio" annotation(Dialog(tab="Initialisation"));

  final parameter Modelica.SIunits.SpecificEnthalpy h_start = TILMedia.GasFunctions.specificEnthalpy_pTxi(medium, p_start, T_start, mixingRatio_initial)
    "Start value for specific Enthalpy inside volume";

  ClaRa.Basics.Units.MassFraction xi[medium.nc - 1](start=mixingRatio_initial);
  Modelica.SIunits.SpecificEnthalpy h(start=h_start) "Specific enthalpy";
  ClaRa.Basics.Units.Pressure p(start=p_start) "Pressure";

  ClaRa.Basics.Units.Mass mass "Gas mass in control volume";

  Real drhodt "Density derivative";

public
  inner Summary    summary(N_ports_out=N_ports_out,inlet(m_flow=inlet.m_flow,  T=gasInlet.T, p=inlet.p, h=gasInlet.h, xi=gasInlet.xi, H_flow=inlet.m_flow*gasInlet.h),
                           outlet(m_flow=outlet.m_flow,  T=gasOutlet.T, p=outlet.p, h=gasOutlet.h, xi=gasOutlet.xi, each H_flow=outlet.m_flow*gasOutlet.h),
                   gas(m=mass, T=bulk.T, p=p, h=h, H=h*mass, rho=bulk.d))
    annotation (Placement(transformation(extent={{-60,-102},{-40,-82}})));

initial equation

  if initType == ClaRa.Basics.Choices.Init.steadyState then
    der(h)=0;
    der(p)=0;
  elseif initType == ClaRa.Basics.Choices.Init.steadyPressure then
    der(p)=0;
  elseif initType == ClaRa.Basics.Choices.Init.steadyEnthalpy then
    der(h)=0;
  elseif initType == ClaRa.Basics.Choices.Init.steadyDensity then
    drhodt=0;
  end if;

equation

  inlet.xi_outflow = xi;
  inlet.T_outflow = bulk.T;

  for i in 1:N_ports_out loop
   outlet[i].p = p;
   outlet[i].T_outflow = bulk.T;
   outlet[i].xi_outflow = xi;
  end for;

   //der(h) = 1/mass*(inlet.m_flow *(gasInlet.h - h) + sum(outlet.m_flow.*(gasOutlet.h - h*ones(N_ports_out))) + volume*der(p)) "Energy balance";
   der(h) =if useHomotopy then homotopy(1/mass*(inlet.m_flow *(gasInlet.h - h) + sum(outlet.m_flow.*(gasOutlet.h - h*ones(N_ports_out))) + volume*der(p)),  1/mass*(m_flow_out_nom[1] *(gasInlet.h - h) + sum(m_flow_out_nom.*(gasOutlet.h - h*ones(N_ports_out))) + volume*der(p))) else 1/mass*(inlet.m_flow *(gasInlet.h - h) + sum(outlet.m_flow.*(gasOutlet.h - h*ones(N_ports_out))) + volume*der(p))
    "Energy balance";

  for i in 1:medium.nc - 1 loop
    //der(xi[i]) = 1/mass.*(inlet.m_flow.*(gasInlet.xi[i]-xi[i]) + sum(outlet.m_flow.*(gasOutlet.xi[i]-xi[i]*ones(N_ports_out)))) "Mass balance";
   der(xi[i]) = if useHomotopy then homotopy(1/mass.*(inlet.m_flow.*(gasInlet.xi[i]-xi[i]) + sum(outlet.m_flow.*(gasOutlet.xi[i]-xi[i]*ones(N_ports_out)))), 1/mass.*(m_flow_out_nom[1].*(xi_nom[i]-xi[i]) +  sum(m_flow_out_nom.*(xi_nom[i]*ones(N_ports_out)-xi[i]*ones(N_ports_out))))) else 1/mass.*(inlet.m_flow.*(gasInlet.xi[i]-xi[i]) + sum(outlet.m_flow.*(gasOutlet.xi[i]-xi[i]*ones(N_ports_out))))
      "Mass balance";
  end for;

      //______________ Balance euqations _______________________

    mass = if useHomotopy then volume*homotopy(bulk.d,rho_nom) else volume*bulk.d
    "Mass in control volume";

    drhodt = bulk.drhodh_pxi*der(h) + bulk.drhodp_hxi*der(p) + sum({bulk.drhodxi_ph[i] * der(bulk.xi[i]) for i in 1:medium.nc-1});

    drhodt*volume = inlet.m_flow + sum(outlet.m_flow) "Mass balance";

    inlet.p = p "Momentum balance";

  annotation (Diagram(coordinateSystem(extent={{-100,-100},{100,100}},
          preserveAspectRatio=true)),
                                 Icon(coordinateSystem(extent={{-100,-100},{100,
            100}}, preserveAspectRatio=true),
        graphics));
end SplitGas_L2_flex;
