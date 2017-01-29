within Exergy.XClaRa.Components.TurboMachines.Pumps;
model PumpVLE_L1_simple
  "A pump for VLE mixtures with a volume flow rate depending on drive power and pressure difference only"
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

extends Exergy.XClaRa.Components.TurboMachines.Pumps.Pump_Base;
  ClaRa.Basics.Interfaces.Connected2SimCenter connected2SimCenter(
    powerIn=0,
    powerOut=inlet.m_flow*(fluidIn.h - fluidOut.h),
    powerAux=(P_drive - inlet.m_flow*(fluidOut.h - fluidIn.h))) if                                                                                                     contributeToCycleSummary;
extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L1");

  parameter Real eta_mech = 0.98 "Mechanic efficiency of the drive"
   annotation(Dialog(group="Part Load and Efficiency"));
parameter ClaRa.Basics.Units.Pressure Delta_p_eps=100
    "|Expert Settings| Numerical Robustnes|Small pressure difference for linearisation around zero";

  Modelica.Blocks.Interfaces.RealInput P_drive
    "Power input of the pump's motor" annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={0,120})));

public
    outer Exergy.Utilities.RefEnv refEnv "reference enviroment properties";

  Exergy.Utilities.ViewObjectNE viewObject(nEnergy={2,0,1,0});

  Exergy.Utilities.ViewPort viewOutput
    annotation (Placement(transformation(extent={{86,88},{106,108}})));
equation

  //  viewObject.h[1].E_flow = inlet.m_flow*noEvent(actualStream(inlet.h_outflow));
  viewObject.h[1].E_flow = inlet.m_flow*fluidIn.h;
  viewObject.h[1].Ex_flow = inlet.m_flow*(fluidIn.h-refEnv.T*fluidIn.s);

  viewObject.h[2].E_flow = outlet.m_flow*fluidOut.h;
  viewObject.h[2].Ex_flow = outlet.m_flow*(fluidOut.h-refEnv.T*fluidOut.s);

    viewObject.w[1].E_flow = P_drive;
    viewObject.w[1].Ex_flow = P_drive;

  connect(viewObject.viewOutput,viewOutput);

public
  Fundamentals.Summary
          summary(outline(V_flow=V_flow, P_hyd=P_hyd,Delta_p=Delta_p, head= Delta_p/(fluidIn.d*Modelica.Constants.g_n), NPSHa = (inlet.p - fluidIn.VLE.p_l)/(fluidIn.d*Modelica.Constants.g_n), eta_hyd=0, eta_mech=eta_mech),
                  inlet(showExpertSummary=showExpertSummary,m_flow=inlet.m_flow, T=fluidIn.T, p=inlet.p, h=fluidIn.h, s=fluidIn.s, steamQuality = fluidIn.q, H_flow= fluidIn.h*inlet.m_flow,  rho=fluidIn.d),
    outlet(
      showExpertSummary=showExpertSummary,
      m_flow=-outlet.m_flow,
      T=fluidOut.T,
      p=outlet.p,
      h=fluidOut.h,
      s=fluidOut.s,
      steamQuality=fluidOut.q,
      H_flow=-fluidOut.h*outlet.m_flow,
      rho=fluidOut.d))                                                                                                     annotation(Placement(transformation(
        extent={{-10,-11},{10,11}},
        origin={-70,-91})));
equation
  P_hyd=P_drive*eta_mech;
  V_flow= P_hyd/(Delta_p+Delta_p_eps);
  inlet.m_flow=V_flow * fluidIn.d;
  inlet.h_outflow=inStream(outlet.h_outflow); // This is a dummy - flow reversal is not supported!
//____________________ Balance equations ___________________
  inlet.m_flow + outlet.m_flow = 0.0 "Mass balance";
  Delta_p=outlet.p-inlet.p "Momentum balance";
//   inStream(inlet.h_outflow) + outlet.h_outflow + P_hyd/inlet.m_flow/eta_hyd = 0.0
//     "Energy balance";
  outlet.h_outflow = inStream(inlet.h_outflow)  + P_drive*eta_mech/(inlet.m_flow+1e-6)
    "Energy balance";

  annotation (Icon(graphics), Diagram(graphics));
end PumpVLE_L1_simple;
