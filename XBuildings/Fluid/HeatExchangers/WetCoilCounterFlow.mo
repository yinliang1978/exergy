within Exergy.XBuildings.Fluid.HeatExchangers;
model WetCoilCounterFlow
  //  extends Buildings.Fluid.HeatExchangers.WetCoilCounterFlow(
  extends Exergy.XBuildings.Fluid.HeatExchangers.WetCoilCounterFlow_X(
      replaceable package Medium1 =
     Exergy.XBuildings.Media.Water,
      replaceable package Medium2 =
      Exergy.XBuildings.Media.Air);

  outer Exergy.Utilities.RefEnv refEnv "reference enviroment properties";

//  Utilities.ViewObject viewObject(nEnergy={4+nEle,1,0,0});
   Utilities.ViewObject viewObject(nEnergy={4,1,0,0});
 Utilities.ViewPort viewOutput
   annotation (Placement(transformation(extent={{86,88},{106,108}})));

equation
  viewObject.h[1].E_flow = port_a1.m_flow*noEvent(actualStream(port_a1.h_outflow));
  viewObject.h[1].Ex_flow = port_a1.m_flow*Medium1.specificExergy(
                                   Medium1.setState_phX(port_a1.p,
                                        noEvent(actualStream(port_a1.h_outflow)),
                                        noEvent(actualStream(port_a1.Xi_outflow))),
                                   refEnv,
                                   Exergy.Utilities.Types.ExTypes.Flow,
                                   Exergy.Utilities.Types.EqmTypes.TM);

  viewObject.h[2].E_flow = port_b1.m_flow*noEvent(actualStream(port_b1.h_outflow));
  viewObject.h[2].Ex_flow = port_b1.m_flow*Medium1.specificExergy(
                                   Medium1.setState_phX(port_b1.p,
                                        noEvent(actualStream(port_b1.h_outflow)),
                                        noEvent(actualStream(port_b1.Xi_outflow))),
                                   refEnv,
                                   Exergy.Utilities.Types.ExTypes.Flow,
                                   Exergy.Utilities.Types.EqmTypes.TM);

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

//  viewObject.h[5].E_flow = mWat_flow*Medium2.specificEnthalpy(
//                                   Medium2.setState_pTX(port_b2.p,
//                                        ele[3].vol2.heatPort.T));
//  viewObject.h[5].Ex_flow = 0;
/*
       for i in 1:nEle loop
  // viewObject.h[4+i].E_flow = mWat_ele_flow[i]*Medium2.enthalpyOfCondensingGas(ele[i].vol2.heatPort.T);
 // viewObject.h[4+i].E_flow=mWat_ele_flow[i]*Medium1.enthalpyOfLiquid(ele[i].vol2.heatPort.T);
//    viewObject.h[4+i].E_flow= mWat_ele_flow[i]*(Buildings.Utilities.Psychrometrics.Constants.cpSte*(ele[i].vol2.heatPort.T
  //                          -Exergy.XBuildings.Media.Air.reference_T)+Medium1.enthalpyOfLiquid(ele[i].vol2.heatPort.T));
 viewObject.h[4+i].E_flow=mWat_ele_flow[i]*Medium1.enthalpyOfLiquid(Exergy.XBuildings.Media.Air.reference_T);//mWat_ele_flow[i]*(Buildings.Utilities.Psychrometrics.Constants.cpSte*(ele[i].vol2.T-Exergy.XBuildings.Media.Air.reference_T));

                          //  -Exergy.XBuildings.Media.Air.reference_T)+Medium1.enthalpyOfLiquid(ele[i].vol2.heatPort.T));
  //(Medium2.specificEnthalpy(
    //                                  Medium2.setState_pTX(port_b2.p,
      //                                 ele[i].vol2.heatPort.T))-Medium2.specificEnthalpy(
        //                              Medium2.setState_pTX(port_b2.p,
          //                             Exergy.XBuildings.Media.Air.reference_T)));

  viewObject.h[4+i].Ex_flow = port_b2.m_flow*Medium2.specificExergy(
                                   Medium2.setState_pTX(port_b2.p,
                                        ele[i].vol2.heatPort.T),
                                   refEnv,
                                   Exergy.Utilities.Types.ExTypes.Flow);

  // mWat_ele_flow[i]=ele[i].vol2.mWat_flow;
 //  mWat_ele_T[i]=ele[i].vol2.heatPort.T;
  //  mWat_P
  end for;
  */

    viewObject.E[1].E = viewRoute.viewTotal.E.E;
  viewObject.E[1].Ex = viewRoute.viewTotal.E.Ex;

    //view out connect
//  connect(viewObject.viewOutput,viewOutput);
 // connect(viewRoute.viewTotal, viewOutput) annotation (Line(points={{81.8,91.9},
   //       {96,91.9},{96,98}}, color={28,108,200}));
  connect(viewRoute.viewTotal, viewOutput) annotation (Line(points={{
          81.8,91.9},{96,91.9},{96,98}}, color={28,108,200}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})));
end WetCoilCounterFlow;
