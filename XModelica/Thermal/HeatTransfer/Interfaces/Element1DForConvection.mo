within Exergy.XModelica.Thermal.HeatTransfer.Interfaces;
partial model Element1DForConvection
  "Partial heat transfer element with two HeatPort connectors that does not store energy"

 // Modelica.SIunits.HeatFlowRate Q_flow "Heat flow rate from port_a -> port_b";
 // Modelica.SIunits.TemperatureDifference dT "port_a.T - port_b.T";
public
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a solid annotation (Placement(transformation(extent={{
            -110,-10},{-90,10}}, rotation=0)));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b fluid annotation (Placement(transformation(extent={{
            90,-10},{110,10}}, rotation=0)));

  outer Exergy.Utilities.RefEnv refEnv "reference enviroment properties";

  Utilities.ViewObject viewObject(nEnergy={0,0,0,2});

  Utilities.ViewPort viewOutput
    annotation (Placement(transformation(extent={{86,88},{106,108}})));
equation
  //heat

    viewObject.q[1].E_flow = solid.Q_flow;
  viewObject.q[1].Ex_flow = solid.Q_flow*(1-refEnv.state.T/solid.T);

  viewObject.q[2].E_flow = fluid.Q_flow;
  viewObject.q[2].Ex_flow = fluid.Q_flow*(1-refEnv.state.T/fluid.T);
    //view out connect
  connect(viewObject.viewOutput,viewOutput);
end Element1DForConvection;
