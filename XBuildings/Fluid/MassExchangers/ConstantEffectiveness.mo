within Exergy.XBuildings.Fluid.MassExchangers;
model ConstantEffectiveness

 extends Buildings.Fluid.MassExchangers.ConstantEffectiveness(
      replaceable package Medium1 =
      Exergy.XBuildings.Media.Air,
      replaceable package Medium2 =
      Exergy.XBuildings.Media.Air);

  outer Exergy.Utilities.RefEnv refEnv "reference enviroment properties";

  Utilities.ViewObject viewObject(nEnergy={4,0,0,0});
 Utilities.ViewPort viewOutput
   annotation (Placement(transformation(extent={{86,88},{106,108}})));

  //show phi
//protected
  Real  phi_a1=min(1,Buildings.Utilities.Psychrometrics.Functions.phi_pTX(sta_a1.p,sta_a1.T,sta_a1.X[1])) if show_T
    "Block to compute water vapor phi at ports";
    Real  phi_b1=min(1,Buildings.Utilities.Psychrometrics.Functions.phi_pTX(sta_b1.p,sta_b1.T,sta_b1.X[1])) if show_T
    "Block to compute water vapor phi at ports";
  Real  phi_a2=min(1,Buildings.Utilities.Psychrometrics.Functions.phi_pTX(sta_a2.p,sta_a2.T,sta_a2.X[1])) if show_T
    "Block to compute water vapor phi at ports";
    Real  phi_b2=min(1,Buildings.Utilities.Psychrometrics.Functions.phi_pTX(sta_b2.p,sta_b2.T,sta_b2.X[1])) if show_T
    "Block to compute water vapor phi at ports";
equation

  // phi_a1=Buildings.Utilities.Psychrometrics.Functions.phi_pTX(sta_a1.p,sta_a1.T,sta_a1.X[1]) if show_T;
  // phi_b1=Buildings.Utilities.Psychrometrics.Functions.phi_pTX(sta_b1.p,sta_b1.T,sta_b1.X[1]) if show_T;
   //phi_a2=Buildings.Utilities.Psychrometrics.Functions.phi_pTX(sta_a2.p,sta_a2.T,sta_a2.X[1]) if show_T;
   //phi_b2=Buildings.Utilities.Psychrometrics.Functions.phi_pTX(sta_b2.p,sta_b2.T,sta_b2.X[1]) if show_T;
  //sta_a1

  viewObject.h[1].E_flow = port_a1.m_flow*noEvent(actualStream(port_a1.h_outflow));
  viewObject.h[1].Ex_flow = port_a1.m_flow*Medium1.specificExergy(
                                   Medium1.setState_phX(port_a1.p,
                                        noEvent(actualStream(port_a1.h_outflow)),
                                        noEvent(actualStream(port_a1.Xi_outflow))),
                                   refEnv,
                                   Exergy.Utilities.Types.ExTypes.Flow);

  viewObject.h[2].E_flow = port_b1.m_flow*noEvent(actualStream(port_b1.h_outflow));
  viewObject.h[2].Ex_flow = port_b1.m_flow*Medium1.specificExergy(
                                   Medium1.setState_phX(port_b1.p,
                                        noEvent(actualStream(port_b1.h_outflow)),
                                        noEvent(actualStream(port_b1.Xi_outflow))),
                                   refEnv,
                                   Exergy.Utilities.Types.ExTypes.Flow);

  viewObject.h[3].E_flow = port_a2.m_flow*noEvent(actualStream(port_a2.h_outflow));
  viewObject.h[3].Ex_flow = port_a2.m_flow*Medium2.specificExergy(
                                   Medium2.setState_phX(port_a2.p,
                                        noEvent(actualStream(port_a2.h_outflow)),
                                        noEvent(actualStream(port_a2.Xi_outflow))),
                                   refEnv,
                                   Exergy.Utilities.Types.ExTypes.Flow);

  viewObject.h[4].E_flow = port_b2.m_flow*noEvent(actualStream(port_b2.h_outflow));
  viewObject.h[4].Ex_flow = port_b2.m_flow*Medium2.specificExergy(
                                   Medium2.setState_phX(port_b2.p,
                                        noEvent(actualStream(port_b2.h_outflow)),
                                        noEvent(actualStream(port_b2.Xi_outflow))),
                                   refEnv,
                                   Exergy.Utilities.Types.ExTypes.Flow);

    //view out connect
   connect(viewObject.viewOutput,viewOutput);
end ConstantEffectiveness;
