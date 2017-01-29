within Exergy.XClaRa.Components.MechanicalSeparation;
model SteamSeparatorVLE_L1
  extends ClaRa.Basics.Icons.Cyclone;
  extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L1");
  import SZT = ClaRa.Basics.Functions.SmoothZeroTransition;

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=simCenter.fluid1
    "Medium in the component"
    annotation (Dialog(group="Fundamental Definitions"), choicesAllMatching);

  parameter Real eta_vap(min=0.5, max=0.999) = 0.96
    "Efficiency of separating liquid from vapour (at steam outlet)"                                                  annotation (Dialog(group="Fundamental Definitions"));
  parameter Real eta_liq(min=0.5, max=0.999) = 0.98
    "Efficiency of separating vapour from liquid (at liquid outlet)"                                                  annotation (Dialog(group="Fundamental Definitions"));

  parameter Boolean showExpertSummary=simCenter.showExpertSummary
    "True, if expert summary should be applied"                                                               annotation(Dialog(group="Summary and Visualisation"));
  parameter Boolean showData=true
    "True, if a data port containing p,T,h,s,m_flow shall be shown, else false"
                                                                                              annotation(Dialog(group="Summary and Visualisation"));
  constant Real eps=1e-3;
 model Summary
  extends ClaRa.Basics.Icons.RecordIcon;

   ClaRa.Basics.Records.FlangeVLE inlet;
   ClaRa.Basics.Records.FlangeVLE outlet1;
   ClaRa.Basics.Records.FlangeVLE outlet2;
 end Summary;
    outer ClaRa.SimCenter simCenter;

  ClaRa.Basics.Interfaces.FluidPortOut outlet1(Medium=medium) "Liquid outlet"
                                                                             annotation (Placement(transformation(extent={{-10,-110},{10,-90}})));
  ClaRa.Basics.Interfaces.FluidPortIn inlet(Medium=medium) "Inlet port" annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
  ClaRa.Basics.Interfaces.FluidPortOut outlet2(Medium=medium) "Steam outlet"  annotation (Placement(transformation(extent={{-10,90},{10,110}})));
  TILMedia.VLEFluid_ph fluidIn(
    vleFluidType=medium,
    h=inStream(inlet.h_outflow),
    p=inlet.p,
    xi=inStream(inlet.xi_outflow)) annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));

  TILMedia.VLEFluid_ph fluidOut1(
    vleFluidType=medium,
    h=actualStream(outlet1.h_outflow),
    p=outlet1.p,
    xi=actualStream(outlet1.xi_outflow)) annotation (Placement(transformation(extent={{-10,-80},{10,-60}})));

  TILMedia.VLEFluid_ph fluidOut2(
    vleFluidType=medium,
    h=actualStream(outlet2.h_outflow),
    p=outlet2.p,
    xi=actualStream(outlet2.xi_outflow)) annotation (Placement(transformation(extent={{-10,60},{10,80}})));

  ClaRa.Basics.Interfaces.EyeOut eye_out2 if
                                            showData annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={40,-100}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={40,-110})));

  ClaRa.Basics.Interfaces.EyeOut eye_out1 if
                                            showData annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=270,
        origin={40,100}), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=270,
        origin={40,110})));

  Summary summary(
    inlet(
      showExpertSummary=showExpertSummary,
      m_flow=inlet.m_flow,
      p=inlet.p,
      h=fluidIn.h,
      T=fluidIn.T,
      s=fluidIn.s,
      steamQuality=fluidIn.q,
      H_flow=inlet.m_flow*fluidIn.h,
      rho=fluidIn.d),
    outlet1(
      showExpertSummary=showExpertSummary,
      m_flow=-outlet1.m_flow,
      p=outlet1.p,
      h=fluidOut1.h,
      T=fluidOut1.T,
      s=fluidOut1.s,
      steamQuality=fluidOut1.q,
      H_flow=-outlet1.m_flow*fluidOut1.h,
      rho=fluidOut1.d),
    outlet2(
      showExpertSummary=showExpertSummary,
      m_flow=-outlet2.m_flow,
      p=outlet2.p,
      h=fluidOut2.h,
      T=fluidOut2.T,
      s=fluidOut2.s,
      steamQuality=fluidOut2.q,
      H_flow=-outlet2.m_flow*fluidOut2.h,
      rho=fluidOut2.d))
                       annotation (Placement(transformation(extent={{-60,-60},{-80,-40}})));

protected
  ClaRa.Basics.Interfaces.EyeIn eye_int
    annotation (Placement(transformation(extent={{39,-53},{41,-51}})));

  ClaRa.Basics.Interfaces.EyeIn eye_int1
    annotation (Placement(transformation(extent={{39,77},{41,75}})));

equation
  outlet2.h_outflow= SZT(fluidIn.h,
                         SZT(SZT((eta_liq*inlet.m_flow*fluidIn.q*fluidIn.VLE.h_v + (1-eta_vap)*inlet.m_flow*(1-fluidIn.q)*fluidIn.VLE.h_l)/noEvent(max(1e-3, inlet.m_flow*(eta_liq*fluidIn.q + eta_vap*fluidIn.q  + 1 - eta_vap - fluidIn.q))),
                                 fluidIn.h,
                                 fluidIn.VLE.h_v-fluidIn.h,
                                 eps),
                             fluidIn.h,
                             fluidIn.h-fluidIn.VLE.h_l,
                             eps),
                         inlet.p-fluidIn.crit.p,
                         eps);
  outlet1.h_outflow= SZT(fluidIn.h,
                         SZT(SZT((eta_vap*inlet.m_flow*(1-fluidIn.q)*fluidIn.VLE.h_l + (1-eta_liq)*inlet.m_flow*fluidIn.q*fluidIn.VLE.h_v)/noEvent(max(1e-3, inlet.m_flow*(eta_vap + fluidIn.q - eta_vap*fluidIn.q - eta_liq*fluidIn.q))),
                                 fluidIn.h,
                                 fluidIn.VLE.h_v-fluidIn.h,
                                 eps),
                             fluidIn.h,
                             fluidIn.h-fluidIn.VLE.h_l,
                             eps),
                         inlet.p-fluidIn.crit.p,
                         1e-6);

  inlet.h_outflow= fluidIn.VLE.h_l; //backflow not supported!

  inlet. p = outlet2.p;
  outlet2.m_flow = - inlet.m_flow*(eta_liq*fluidIn.q + eta_vap*fluidIn.q + 1 - eta_vap - fluidIn.q);
  // outlet1.m_flow = - inlet.m_flow*fluidIn.q;
  0 = inlet.m_flow + outlet1.m_flow + outlet2.m_flow;

  eye_int.m_flow=-outlet2.m_flow;
  eye_int.T=summary.outlet2.T-273.15;
  eye_int.s=fluidOut2.s/1000;
  eye_int.h=summary.outlet2.h/1000;
  eye_int.p=summary.outlet2.p/100000;

  eye_int1.m_flow=-outlet1.m_flow;
  eye_int1.T=summary.outlet1.T-273.15;
  eye_int1.s=fluidOut1.s/1000;
  eye_int1.h=summary.outlet1.h/1000;
  eye_int1.p=summary.outlet1.p/100000;

  connect(eye_int, eye_out2) annotation (Line(
      points={{40,-52},{40,-100}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(eye_int1, eye_out1) annotation (Line(
      points={{40,76},{40,100}},
      color={190,190,190},
      smooth=Smooth.None));
end SteamSeparatorVLE_L1;
