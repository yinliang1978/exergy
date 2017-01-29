within Exergy.XBuildings.Fluid.HeatExchangers;
model ConstantEffectiveness
  extends Buildings.Fluid.HeatExchangers.ConstantEffectiveness(
      replaceable package Medium1 =
     Exergy.XModelica.Media.PartialMedium,
      replaceable package Medium2 =
     Exergy.XModelica.Media.PartialMedium);

  outer Exergy.Utilities.RefEnv refEnv "reference enviroment properties";

  Utilities.ViewObject viewObject(nEnergy={4,0,0,0});
 Utilities.ViewPort viewOutput
   annotation (Placement(transformation(extent={{86,88},{106,108}})));

equation
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
