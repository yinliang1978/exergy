within Exergy.XClaRa.Components.VolumesValvesFittings.Valves;
model ThreeWayValveVLE_L1
  "Three way valve for vle media, not suitable for back flows!"
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
  extends
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.ThreeWayValve_base;
  extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L1");

// record Outline
//   extends ClaRa.Basics.Icons.RecordIcon;
//     ClaRa.Basics.Units.Volume volume_tot "Total volume";
// end Outline;

  replaceable model PressureLoss =
      Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.QuadraticFrictionFlowAreaSymetric_TWV
    constrainedby
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.TWV_L1
    "Pressure loss model"                                                                   annotation(choicesAllMatching);

public
  parameter SI.Area effectiveFlowArea=7.85e-3 "Effective flow area for outlets"
    annotation(Dialog(group="Valve Parameters"));
//   parameter SI.Area effectiveFlowArea2=effectiveFlowArea1 "Effective flow area for outlet 2"
//     annotation(Dialog(group="Valve Parameters"));

  parameter SI.PressureDifference Delta_p_smooth = 100
    "Below this value, root function is approximated linearly"                   annotation(Dialog(tab = "Expert Settings", group="Numerical Robustness"));

parameter Boolean useStabilisedMassFlow=false
    "|Expert Settings|Numerical Robustness|";
    parameter SI.Time Tau= 0.001 "Time Constant of Stabilisation" annotation(Dialog(tab="Expert Settings", group = "Numerical Robustness", enable=useStabilisedMassFlow));

public
  model Summary
    extends ClaRa.Basics.Icons.RecordIcon;
    //     Outline outline;
    ClaRa.Basics.Records.FlangeVLE           inlet;
    ClaRa.Basics.Records.FlangeVLE           outlet1;
    ClaRa.Basics.Records.FlangeVLE           outlet2;
  end Summary;

   Summary summary(inlet(showExpertSummary = showExpertSummary,m_flow=inlet.m_flow,  T=fluidIn.T, p=inlet.p, h=fluidIn.h,s=fluidIn.s, steamQuality=fluidIn.q, H_flow=fluidIn.h*inlet.m_flow, rho=fluidIn.d),
                   outlet1(showExpertSummary = showExpertSummary,m_flow = -outlet1.m_flow, T=fluidOut1.T, p=outlet1.p, h=fluidOut1.h, s=fluidOut1.s, steamQuality=fluidOut1.q, H_flow=-fluidOut1.h*outlet1.m_flow, rho=fluidOut1.d),
                   outlet2(showExpertSummary = showExpertSummary,m_flow = -outlet2.m_flow, T=fluidOut2.T, p=outlet2.p, h=fluidOut2.h, s=fluidOut2.s, steamQuality=fluidOut2.q, H_flow=-fluidOut2.h*outlet2.m_flow, rho=fluidOut2.d))
    annotation (Placement(transformation(extent={{-60,-102},{-40,-82}})));

    inner
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.ICom_TWV
    iCom(
    p_in=inlet.p,
    p_out1=outlet1.p,
    p_out2=outlet2.p,
    rho_out1=fluidOut1.d,
    rho_out2=fluidOut2.d,
    opening_=splitRatio,
    m_flow_in=inlet.m_flow,
    rho_in=fluidIn.d) annotation (Placement(transformation(extent={{-80,
            -100},{-60,-80}})));

    PressureLoss pressureLoss(useStabilisedMassFlow=useStabilisedMassFlow,
                              Tau=Tau,
                              Delta_p_smooth=Delta_p_smooth,
                              effectiveFlowArea1=effectiveFlowArea) annotation (Placement(transformation(extent={{-40,-40},{-20,-20}})));

equation
  outlet1.m_flow = -pressureLoss.m_flow_1;
  outlet2.m_flow = -pressureLoss.m_flow_2;
  inlet.m_flow + outlet1.m_flow + outlet2.m_flow =0;

  // Isenthalpic state transformation (no storage and no loss of energy)
  inlet.h_outflow = (inStream(outlet1.h_outflow)*outlet1.m_flow + inStream(outlet2.h_outflow)*outlet2.m_flow)/(-inlet.m_flow);
  outlet1.h_outflow = inStream(inlet.h_outflow);
  outlet2.h_outflow = inStream(inlet.h_outflow);

  // No chemical reaction taking place:
  inlet.xi_outflow = (inStream(outlet1.xi_outflow)*outlet1.m_flow + inStream(outlet2.xi_outflow)*outlet2.m_flow)/(-inlet.m_flow);
  outlet1.xi_outflow = inStream(inlet.xi_outflow);
  outlet2.xi_outflow = inStream(inlet.xi_outflow);

annotation (
  Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,80}},
        grid={2,2}), graphics));
end ThreeWayValveVLE_L1;
