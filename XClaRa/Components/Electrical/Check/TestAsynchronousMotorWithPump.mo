within Exergy.XClaRa.Components.Electrical.Check;
model TestAsynchronousMotorWithPump
 extends
    Exergy.XClaRa.Components.TurboMachines.Pumps.Check.TestPump_L1_WithEMotor;
 extends ClaRa.Basics.Icons.PackageIcons.ExecutableRegressiong100;
  model Regression
      extends ClaRa.Basics.Icons.RegressionSummary;

    Modelica.Blocks.Interfaces.RealInput V_flow;
    Modelica.Blocks.Interfaces.RealInput m_flow;
    Modelica.Blocks.Interfaces.RealInput U_term;
    Modelica.Blocks.Interfaces.RealInput rpm;
    Modelica.Blocks.Interfaces.RealInput tau;
    Real y_V_flow_min = timeExtrema.y_min;
    Real y_V_flow_max = timeExtrema.y_max;
    Real y_rpm_min = timeExtrema1.y_min;
    Real y_rpm_max = timeExtrema1.y_max;
    Real y_m_flow_int = integrator.y;
    Real y_U_term_int = integrator1.y;
    Real y_tau_int = integrator2.y;

  protected
    Utilities.Blocks.TimeExtrema timeExtrema(u=V_flow) annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
    Utilities.Blocks.Integrator integrator(u = m_flow) annotation (Placement(transformation(extent={{-40,-20},{-20,0}})));
    Utilities.Blocks.TimeExtrema timeExtrema1(u = rpm)  annotation (Placement(transformation(extent={{-40,60},{-20,80}})));
    Utilities.Blocks.Integrator integrator1(u = U_term)  annotation (Placement(transformation(extent={{-40,-50},{-20,-30}})));
    Utilities.Blocks.Integrator integrator2(u = tau) annotation (Placement(transformation(extent={{-40,-84},{-20,-64}})));

  end Regression;

  Regression regression(
    V_flow = pump.V_flow,
    m_flow = pump.summary.inlet.m_flow,
    U_term = motor.U_term,
    rpm = motor.summary.rpm,
    tau = motor.summary.tau_mech) annotation (Placement(transformation(extent={{-100,-40},{-80,-20}})));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={115,150,0},
          lineThickness=0.5)}));
end TestAsynchronousMotorWithPump;
