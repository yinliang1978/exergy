within Exergy.XModelica.Thermal.HeatTransfer.Interfaces;
partial model Element1Port
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b  port;

  outer Exergy.Utilities.RefEnv refEnv "reference enviroment properties";

  Utilities.ViewObject viewObject(nEnergy={0,1,0,1});
  Utilities.ViewPort viewOutput
    annotation (Placement(transformation(extent={{86,88},{106,108}})));
equation
  //heat

  viewObject.q[1].E_flow = port.Q_flow;
  viewObject.q[1].Ex_flow = port.Q_flow*(1-refEnv.state.T/port.T);

  viewObject.e[1].E_flow = port.Q_flow;
  viewObject.e[1].Ex_flow = port.Q_flow*(1-refEnv.state.T/port.T);

    //view out connect
  connect(viewObject.viewOutput,viewOutput);
end Element1Port;
