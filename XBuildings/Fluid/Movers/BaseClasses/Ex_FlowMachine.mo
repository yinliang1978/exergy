within Exergy.XBuildings.Fluid.Movers.BaseClasses;
partial model Ex_FlowMachine

  outer Exergy.Utilities.RefEnv refEnv "reference enviroment properties";

  Utilities.ViewObject viewObject(nEnergy={2,1,1,1});

  Utilities.ViewPort viewOutput
    annotation (Placement(transformation(extent={{86,88},{106,108}})));
equation
  //heat

  viewObject.h[1].E_flow = port_a.m_flow*noEvent(actualStream(port_a.h_outflow));
  viewObject.h[1].Ex_flow = port_a.m_flow*Medium.specificEnthalpyExergy(
                                   Medium.setState_phX(port_a.p,
                                        noEvent(actualStream(port_a.h_outflow)),
                                        noEvent(actualStream(port_a.Xi_outflow))),
                                   refEnv.state,
                                   Exergy.Utilities.Types.ExTypes.EnthalpyExergy);

  viewObject.h[2].E_flow = port_b.m_flow*noEvent(actualStream(port_b.h_outflow));
  viewObject.h[2].Ex_flow = port_b.m_flow*Medium.specificEnthalpyExergy(
                                   Medium.setState_phX(port_b.p,
                                        noEvent(actualStream(port_b.h_outflow)),
                                        noEvent(actualStream(port_b.Xi_outflow))),
                                   refEnv.state,
                                   Exergy.Utilities.Types.ExTypes.EnthalpyExergy);

  if useSteadyStateTwoPort then
  viewObject.e[1].E_flow = 0;
  viewObject.e[1].Ex_flow = 0;
  else

  viewObject.e[1].E_flow = vol.dynBal.U;
  viewObject.e[1].Ex_flow = vol.dynBal.m*Medium.specificEnthalpyExergy(
                                   vol.dynBal.medium.state,
                                   refEnv.state,
                                   Exergy.Utilities.Types.ExTypes.InternalEnergyExergy);
  end if;

  viewObject.q[1].E_flow = heatPort.Q_flow;
  viewObject.q[1].Ex_flow = heatPort.Q_flow*(1-refEnv.state.T/heatPort.T);

  viewObject.w[1].E_flow = P;
  viewObject.w[1].Ex_flow = P;

    //view out connect
  connect(viewObject.viewOutput,viewOutput);
end Ex_FlowMachine;
