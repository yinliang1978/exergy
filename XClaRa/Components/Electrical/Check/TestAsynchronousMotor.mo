within Exergy.XClaRa.Components.Electrical.Check;
model TestAsynchronousMotor "A simple test for the simple motor"

  extends ClaRa.Basics.Icons.PackageIcons.ExecutableExampleb80;

  inner ClaRa.SimCenter simCenter(redeclare
      TILMedia.VLEFluidTypes.TILMedia_SplineWater                                       fluid1, showExpertSummary=false) annotation (Placement(transformation(extent={{-100,-100},{-80,-80}})));
  Modelica.Mechanics.Rotational.Components.Inertia inertia1(         w(start=10), J=50)
    annotation (Placement(transformation(extent={{-48,-80},{-28,-60}})));
  Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={0,-16})));
  Modelica.Blocks.Math.Gain              realExpression(k=0.5*2*Modelica.Constants.pi
        /60)
    annotation (Placement(transformation(extent={{26,28},{18,36}})));
  Modelica.Blocks.Sources.TimeTable
                               ramp1(table=[0,5100; 100,5100; 101,4600; 200,4600; 201,4100; 300,4100; 301,3600; 400,3600; 401,3100; 500,3100; 501,2950*2; 600,2950*2])
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=180,
        origin={62,32})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder(T=1, initType=Modelica.Blocks.Types.Init.SteadyState)
    annotation (Placement(transformation(extent={{38,28},{30,36}})));
  AsynchronousMotor_L2 motor(
    rpm_nom=2950,
    I_rotor_nom=10,
    J=50,
    cosphi=0.9,
    tau_bd_nom=550,
    P_nom=154.457e3,
    eta_stator=0.9,
    activateHeatPort=true,
    N_pp=1,
    U_term_nom=3e3)
            annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-62,-68})));
  Exergy.XClaRa.Components.Utilities.Blocks.LimPID PID(
    u_ref=500,
    initType=Modelica.Blocks.Types.InitPID.InitialOutput,
    y_max=100e3,
    Tau_d=500,
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    y_min=0.001,
    y_ref=3.8e3,
    y_start=40,
    k=50,
    Tau_i=0.05)
    annotation (Placement(transformation(extent={{10,22},{-10,42}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T=293.15)
    annotation (Placement(transformation(extent={{-9,-9},{9,9}},
        rotation=270,
        origin={-62,-45})));
  Modelica.Mechanics.Rotational.Components.IdealGear idealGear(ratio=1)
    annotation (Placement(transformation(extent={{-24,-80},{-4,-60}})));
  Modelica.Mechanics.Rotational.Sources.Torque torque
    annotation (Placement(transformation(extent={{36,-80},{16,-60}})));
  Modelica.Blocks.Sources.Ramp           Flow2(
    height=250,
    offset=-500,
    duration=100,
    startTime=2500) "-(0.6*200e5/(2*3.14*5100/60))*motor.rpm^2/3000"
    annotation (Placement(transformation(extent={{86,-80},{66,-60}})));
  Modelica.Blocks.Sources.Ramp           Flow1(
    startTime=250,
    duration=100,
    height=0,
    offset=3e3) "-(0.6*200e5/(2*3.14*5100/60))*motor.rpm^2/3000"
    annotation (Placement(transformation(extent={{-52,2},{-72,22}})));
  Modelica.Blocks.Sources.Ramp           Flow3(
    startTime=250,
    duration=100,
    height=0,
    offset=50) "-(0.6*200e5/(2*3.14*5100/60))*motor.rpm^2/3000"
    annotation (Placement(transformation(extent={{-52,-30},{-72,-10}})));
equation
  connect(firstOrder.y, realExpression.u) annotation (Line(
      points={{29.6,32},{26.8,32}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ramp1.y, firstOrder.u) annotation (Line(
      points={{51,32},{38.8,32}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(realExpression.y, PID.u_s) annotation (Line(
      points={{17.6,32},{12,32}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(speedSensor.w, PID.u_m) annotation (Line(
      points={{0,-5},{0,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(speedSensor.flange, idealGear.flange_b) annotation (Line(
      points={{0,-26},{0,-70},{-4,-70}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(inertia1.flange_b, idealGear.flange_a) annotation (Line(
      points={{-28,-70},{-24,-70}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(inertia1.flange_a, motor.shaft) annotation (Line(
      points={{-48,-70},{-50,-70},{-50,-68},{-52,-68}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(Flow2.y,torque. tau) annotation (Line(
      points={{65,-70},{38,-70}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(torque.flange, idealGear.flange_b) annotation (Line(
      points={{16,-70},{-4,-70}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(motor.heat, fixedTemperature.port) annotation (Line(
      points={{-62,-58},{-62,-54}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(Flow3.y, motor.f_term) annotation (Line(
      points={{-73,-20},{-84,-20},{-84,-64},{-74,-64}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(Flow1.y, motor.U_term) annotation (Line(
      points={{-73,12},{-94,12},{-94,-68},{-74,-68}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                               graphics={
                                  Text(
          extent={{-96,100},{102,60}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=10,
          textString="______________________________________________________________________________________________
PURPOSE:

______________________________________________________________________________________________
"),                    Text(
          extent={{-136,104},{64,84}},
          lineColor={0,128,0},
          fontSize=31,
          textString="TESTED -- YYYY-MM-DD //XX"),Text(
          extent={{-96,60},{68,46}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          textString="______________________________________________________________________________________________________________
Remarks: 
______________________________________________________________________________________________________________
",        fontSize=8),Text(
          extent={{-96,74},{104,56}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=10,
          textString="______________________________________________________________________________________________
Scenario:  

______________________________________________________________________________________________
")}),
    experiment(StopTime=600, __Dymola_NumberOfIntervals=50000),
    __Dymola_experimentSetupOutput);
end TestAsynchronousMotor;
