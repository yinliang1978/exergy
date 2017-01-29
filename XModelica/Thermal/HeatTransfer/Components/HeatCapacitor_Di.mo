within Exergy.XModelica.Thermal.HeatTransfer.Components;
model HeatCapacitor_Di
  "test using the difference form to replace the differential"
  extends Modelica.Thermal.HeatTransfer.Components.HeatCapacitor;

  outer Exergy.Utilities.RefEnv refEnv "reference enviroment properties";

  Utilities.ViewObject viewObject(nEnergy={0,1,0,1});

  Utilities.ViewPort viewOutput
    annotation (Placement(transformation(extent={{86,88},{106,108}})));
 // Modelica.SIunits.InternalEnergy x_U;
  // Modelica.SIunits.Entropy x_S;
          Real ex;
          Real e;
        parameter Real deltaT=0.1;
initial equation
  //der(e)=0;
  //der(ex)=0;

equation
  //heat

  viewObject.q[1].E_flow = port.Q_flow;
  viewObject.q[1].Ex_flow = port.Q_flow*(1-refEnv.state.T/port.T);

 // x_U = C*(T-refEnv.state.T);
 // x_S = C*Modelica.Math.log(T/refEnv.state.T);
  //viewObject.internalEnergy[1].E = der(x_U);
  //viewObject.internalEnergy[1].Ex = der(x_U)-refEnv.state.T*der(x_S);
 //  viewObject.E[1].E =  C*(T-refEnv.state.T);
   //viewObject.E[1].Ex = C*(T-refEnv.state.T)-refEnv.state.T* C*Modelica.Math.log(T/refEnv.state.T);

       e=C*(T-refEnv.state.T);
       ex=C*(T-refEnv.state.T)-refEnv.state.T* C*Modelica.Math.log(T/refEnv.state.T);
               viewObject.e[1].E_flow =  (e-delay(e,deltaT))/(deltaT);
        viewObject.e[1].Ex_flow =  (ex-delay(ex,deltaT))/(deltaT);
      // viewObject.e[1].E_flow =  (3*e-4*delay(e,deltaT)+delay(e,2*deltaT))/(2*deltaT);
      // viewObject.e[1].Ex_flow =  (3*ex-4*delay(ex,deltaT)+delay(ex,2*deltaT))/(2*deltaT);

    //view out connect
  connect(viewObject.viewOutput,viewOutput);

    annotation (Placement(transformation(extent={{88,84},{108,104}})), Icon(
        coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
              Diagram(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}})));
end HeatCapacitor_Di;
