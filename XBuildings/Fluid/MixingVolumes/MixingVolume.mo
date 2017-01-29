within Exergy.XBuildings.Fluid.MixingVolumes;
model MixingVolume
  extends Buildings.Fluid.MixingVolumes.MixingVolume(redeclare replaceable
      package Medium =
       Exergy.XModelica.Media.PartialMedium);
  outer Exergy.Utilities.RefEnv refEnv "reference enviroment properties";

  Utilities.ViewObject viewObject(nEnergy={nPorts,1,0,1});

  Utilities.ViewPort viewOutput
    annotation (Placement(transformation(extent={{86,88},{106,108}})));

   parameter Boolean ExuseSteadyStateTwoPort=(nPorts == 2) and
      (prescribedHeatFlowRate or (not allowFlowReversal)) and (
      energyDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) and (
      massDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) and (
      substanceDynamics == Modelica.Fluid.Types.Dynamics.SteadyState) and (
      traceDynamics == Modelica.Fluid.Types.Dynamics.SteadyState);

//protected
        Real ex;
     //   parameter Real deltaT=0.1;
       // Real mWater;
      //  Real e2[nPorts];
      //  Real ex2[nPorts];
      //  Real T2[nPorts];
     //   Medium.ThermodynamicState sta_a2[nPorts];
        final parameter Boolean computeCSen2 = dynBal.CSen > Modelica.Constants.eps;

      //    Real  phis=min(1,Buildings.Utilities.Psychrometrics.Functions.phi_pTX(p,T,Xi[Medium.nXi]));
equation
  //heat
  for i in 1:nPorts loop
  viewObject.h[i].E_flow = ports[i].m_flow*noEvent(actualStream(ports[i].h_outflow));
  viewObject.h[i].Ex_flow = ports[i].m_flow*Medium.specificExergy(
                                   Medium.setState_phX(ports[i].p,
                                        noEvent(actualStream(ports[i].h_outflow)),
                                        noEvent(actualStream(ports[i].Xi_outflow))),
                                   refEnv,
                                   Exergy.Utilities.Types.ExTypes.Flow);
  //e2[i]=-viewObject.h[i].E_flow;
  //ex2[i]=-viewObject.h[i].Ex_flow;
  //sta_a2[i]=Medium.setState_phX(ports[i].p,
  //                        noEvent(actualStream(ports[i].h_outflow)),
  //                        noEvent(actualStream(ports[i].Xi_outflow)));
  //       T2[i]=       sta_a2[i].T;
 // T2[i]=(Medium.setState_phX(ports[i].p,
 //                         noEvent(actualStream(ports[i].h_outflow)),
  //                        noEvent(actualStream(ports[i].Xi_outflow)))).T;
  end for;

  if ExuseSteadyStateTwoPort then
  viewObject.E[1].E = 0;
  viewObject.E[1].Ex = 0;
  //mWater=0;
  else

  viewObject.E[1].E = dynBal.U;
 // mWater=0;

// viewObject.E[1].Ex = dynBal.m*Medium.specificEnthalpyExergy(
  //                                  dynBal.medium.state,
    //                                refEnv.state,
      //                              Exergy.Utilities.Types.ExTypes.InternalEnergyExergy);

                      //  if dynBal.computeCSen then
                        if computeCSen2 then
                                 ex =dynBal.m*Medium.specificExergy(
                                    dynBal.medium.state,
                                    refEnv,
                                    Exergy.Utilities.Types.ExTypes.Vol)+
                                    dynBal.CSen*(dynBal.medium.T-refEnv.state.T)
                                    -refEnv.state.T* dynBal.CSen*Modelica.Math.log(dynBal.medium.T/refEnv.state.T);
                        else

                   ex =dynBal.m*Medium.specificExergy(
                                    dynBal.medium.state,
                                    refEnv,
                                    Exergy.Utilities.Types.ExTypes.Vol);
                        end if;
     viewObject.E[1].Ex = ex;

 // exRate=(viewObject.E[1].Ex-delay(viewObject.E[1].Ex,0.001))/0.001;

  //  viewObject.e[1].Ex_flow = (ex-delay(ex,1e-5))/1e-5;
  //  viewObject.e[1].Ex_flow =  (3*ex-4*delay(ex,deltaT)+delay(ex,2*deltaT))/(2*deltaT);
  //    viewObject.e[1].Ex_flow =  (25*ex-48*delay(ex,deltaT)+36*delay(ex,2*deltaT)-delay(ex,3*deltaT)+3*delay(ex,4*deltaT))/(12*deltaT);

  end if;

  viewObject.q[1].E_flow = heatPort.Q_flow;
  viewObject.q[1].Ex_flow = heatPort.Q_flow*(1-refEnv.state.T/heatPort.T);

    //view out connect
  connect(viewObject.viewOutput,viewOutput);

end MixingVolume;
