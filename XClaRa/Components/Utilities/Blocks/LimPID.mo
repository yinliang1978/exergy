within Exergy.XClaRa.Components.Utilities.Blocks;
block LimPID
  "P, PI, PD, and PID controller with limited output, anti-windup compensation and delayed, smooth activation"
  import Modelica.Blocks.Types.InitPID;
  import Modelica.Blocks.Types.SimpleController;

  output Real controlError = u_s - u_m
    "Control error (set point - measurement)";

//---------------------------------------
//General Design of the Controller ------
  parameter Modelica.Blocks.Types.SimpleController controllerType=
         Modelica.Blocks.Types.SimpleController.PID "Type of controller" annotation(Dialog(group="General Design of Controller"));
  parameter Real sign= 1
    "set to 1 if a positive control error leads to a positive control output, else -1"
                                                                                       annotation(Dialog(group="General Design of Controller"));
  parameter Boolean perUnitConversion= true
    "True, if input and output values should be normalised with respect to reference values"
                                                                                            annotation(Dialog(group="Normalisation of I/O Signals"));
  parameter Real u_ref = 1 "Reference value for controlled variable"
                                                                    annotation(Dialog(enable=perUnitConversion, group="Normalisation of I/O Signals"));
  parameter Real y_ref = 1 "Reference value for actuated variable"
                                                                  annotation(Dialog(enable=perUnitConversion, group="Normalisation of I/O Signals"));
  parameter Real y_max=1 "Upper limit of output" annotation(Dialog(group="Limiter for Controller Output"));
  parameter Real y_min=-y_max "Lower limit of output" annotation(Dialog(group="Limiter for Controller Output"));

//----------------------------------------
//Time Resononse of the Controller -------
  parameter Real k = 1 "Gain of Proportional block"
                                                   annotation(Dialog(group="Time Response of the Controller"));
  parameter Modelica.SIunits.Time Tau_i(min=Modelica.Constants.small)=0.5
    "1/Ti is gain of integrator block"
                                      annotation(Dialog(enable=controllerType==Modelica.Blocks.Types.SimpleController.PI or
                                controllerType==Modelica.Blocks.Types.SimpleController.PID,group="Time Response of the Controller"));
 parameter Modelica.SIunits.Time Tau_d(min=0)=0.1 "Gain of derivative block"
                              annotation(Dialog(enable=controllerType==Modelica.Blocks.Types.SimpleController.PD or
                                controllerType==Modelica.Blocks.Types.SimpleController.PID,group="Time Response of the Controller"));

  parameter Modelica.SIunits.Time Ni(min=100*Modelica.Constants.eps) = 0.9
    "1/Ni is gain of anti-windup compensation"
                                              annotation (Dialog(enable=controllerType==Modelica.Blocks.Types.SimpleController.PI or controllerType==Modelica.Blocks.Types.SimpleController.PID, group="Anti-Windup Compensation"));
  parameter Real Nd = 1
    "The smaller Nd, the more ideal the derivative block, setting Nd=0 introduces ideal derivative"
       annotation(Dialog(enable=controllerType==Modelica.Blocks.Types.SimpleController.PD or
                                controllerType==Modelica.Blocks.Types.SimpleController.PID,group="Derivative Filtering"));

//------------------- Controller activation --------------------

parameter Boolean use_activateInput = false
    "Provide Boolean input to switch controller on/off."
                                                    annotation(Dialog(tab="Controller activation"));
parameter ClaRa.Basics.Units.Time t_activation=0.0
    "Time when controller is switched on. For use_activateInput==true the controller is switched on if (time>t_activation AND activateController=true)."
    annotation (Dialog(tab="Controller activation"));
parameter ClaRa.Basics.Units.Time Tau_lag_I=0.0
    "Time lag for activation of integral part AFTER controller is being switched on "
    annotation (Dialog(tab="Controller activation"));

parameter Real y_inactive = 1 "Controller output if controller is not active" annotation(Dialog(tab="Controller activation"));

//Signal Smoothening---------------------------

public
  parameter Real Tau_in(min=0)=0
    "Time constant for input smoothening, Tau_in=0 refers to signal no smoothening"
      annotation(Dialog(tab="I/O Filters"));
  parameter Real Tau_out(min=0)=0
    "time constant for output smoothening, Tau_out=0 refers to signal no smoothening"
           annotation(Dialog(tab="I/O Filters"));

//Initialisation--------------------------
public
  parameter InitPID initType=Modelica.Blocks.Types.InitPID.DoNotUse_InitialIntegratorState
    "Type of initialization"         annotation (
      Dialog(tab="Initialization"));
  parameter Boolean limitsAtInit = true
    "= false, if limits are ignored during initializiation"
    annotation(Dialog(tab="Initialization",
                       enable=controllerType==SimpleController.PI or
                              controllerType==SimpleController.PID));

  parameter Real xi_start=0
    "Initial or guess value value for integrator output (= integrator state)"
    annotation (Dialog(enable= initType == Modelica.Blocks.Types.InitPID.InitialState or initType == Modelica.Blocks.Types.InitPID.NoInit,  tab="Initialization"));

  parameter Real xd_start=0
    "Initial or guess value for state of derivative block"
    annotation (Dialog(tab="Initialization",
                         enable=((controllerType==Modelica.Blocks.Types.SimpleController.PD or
                                controllerType==Modelica.Blocks.Types.SimpleController.PID) and initType == Modelica.Blocks.Types.InitPID.InitialState or initType == Modelica.Blocks.Types.InitPID.NoInit)));
  parameter Real y_start=0 "Initial value of output"
    annotation(Dialog(enable=(initType == Modelica.Blocks.Types.InitPID.InitialOutput or initType == Modelica.Blocks.Types.InitPID.SteadyState), tab=
          "Initialization"));

//Expert Settings---------------------------------------------------------------
  parameter Real Tau_add(min=0)=0
    "Set to >0 for additional state after add block in controller, if DAE-index reduction fails."
    annotation(Dialog(tab="Expert Settings", group="DAE Index Reduction"));

protected
  parameter Boolean with_I = controllerType==Modelica.Blocks.Types.SimpleController.PI or
                             controllerType==Modelica.Blocks.Types.SimpleController.PID annotation (HideResult=true);
  parameter Boolean with_D = controllerType==Modelica.Blocks.Types.SimpleController.PD or
                             controllerType==Modelica.Blocks.Types.SimpleController.PID annotation (HideResult=true);
  Real resetValueP "Input to P part before controller activation";
  Real resetValueID "Output of controller before activation";
//  parameter Real y_in_start(fixed=false) "Start value of inlet pseudo state";
//   parameter Real y_out_start(fixed=false) "Start value of outlet pseudo state";
//   parameter Real y_aux_start(fixed=false) "Start value of auxilliary pseudo state";
public
  Modelica.Blocks.Interfaces.RealInput u_s "Connector of setpoint input signal"
    annotation (Placement(transformation(extent={{-240.5,-20},{-200.5,20}},
          rotation=0), iconTransformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealInput u_m
    "Connector of measurement input signal"
    annotation (Placement(transformation(
        origin={0,-216},
        extent={{20,-20},{-20,20}},
        rotation=270), iconTransformation(
        extent={{20,-20},{-20,20}},
        rotation=270,
        origin={1,-120})));

    Modelica.Blocks.Interfaces.BooleanInput activateInput if use_activateInput
    "true, if controller is on" annotation (Placement(transformation(extent={{-239.5,140},{-199.5,180}}),
                                  iconTransformation(extent={{-140,-100},{-100,-60}})));

  Modelica.Blocks.Interfaces.RealOutput y "Connector of actuator output signal"
    annotation (Placement(transformation(extent={{260,-10},{280,10}}, rotation=0),
        iconTransformation(extent={{100,-10},{120,10}})));
  Modelica.Blocks.Math.Gain P(k=k)
                     annotation (Placement(transformation(extent={{-20,120},{0,140}}, rotation=0)));
  Exergy.XClaRa.Components.Utilities.Blocks.Integrator I(
    y_startInputIsActive=true,
    Tau_i_const=Tau_i,
    initType=if initType == InitPID.SteadyState then Modelica.Blocks.Types.Init.SteadyState
         else if initType == InitPID.InitialOutput then Modelica.Blocks.Types.Init.InitialOutput
         else if initType == InitPID.InitialState then Modelica.Blocks.Types.Init.InitialState
         else Modelica.Blocks.Types.Init.NoInit) if with_I
    annotation (Placement(transformation(extent={{-30,-92},{-10,-72}},
          rotation=0)));
  Exergy.XClaRa.Components.Utilities.Blocks.DerivativeClaRa D_approx(
    k=Tau_d,
    x_start=xd_start,
    initType=if initType == InitPID.SteadyState or initType ==
        InitPID.InitialOutput then Modelica.Blocks.Types.Init.SteadyState
         else if initType == InitPID.InitialState then Modelica.Blocks.Types.Init.InitialState
         else Modelica.Blocks.Types.Init.NoInit,
    Tau=Nd) if with_D annotation (Placement(transformation(extent={{-28,
            50},{-8,70}}, rotation=0)));
  Modelica.Blocks.Math.Add3 addPID(
    k1=1,
    k2=1,
    k3=1)                 annotation (Placement(transformation(extent={{40,55},{50,65}}, rotation=0)));
  Modelica.Blocks.Math.Add addI if with_I annotation (Placement(transformation(
        extent={{-5,-5},{5,5}},
        rotation=0,
        origin={-60,-105})));
  Modelica.Blocks.Math.Gain gainTrack(k=1/Ni) if   with_I
    annotation (Placement(transformation(extent={{3,-119},{-10,-132}},
                                                                     rotation=0)));
  Modelica.Blocks.Nonlinear.Limiter limiter(
    limitsAtInit=limitsAtInit,
    uMax=if perUnitConversion then y_max/y_ref else y_max,
    uMin=if perUnitConversion then y_min/y_ref else y_min)
    annotation (Placement(transformation(extent={{152,58},{172,78}},rotation=0)));

public
  Modelica.Blocks.Sources.Constant Dzero(k=0) if not with_D
    annotation (Placement(transformation(extent={{-15,34.5},{-8,41.5}},
                                                                     rotation=0)));
  Modelica.Blocks.Sources.Constant Izero(k=0) if not with_I
    annotation (Placement(transformation(
        extent={{-6,-4},{4,6}},
        rotation=0,
        origin={-15.5,-104.75})));
  Modelica.Blocks.Math.Gain toPU(k=if perUnitConversion then sign/u_ref else
        sign) "convert input values to \"per unit\""
                                           annotation (Placement(transformation(
          extent={{-106,-10},{-86,10}},
                                     rotation=0)));
  Modelica.Blocks.Math.Feedback feedback
    annotation (Placement(transformation(extent={{-190.5,-10},{-170.5,10}},
                                                                      rotation=0)));
  Modelica.Blocks.Math.Gain fromPU(k=if perUnitConversion then y_ref else 1)
    "convert output values to \"Real\""    annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
                                   rotation=0,
        origin={210,0})));

  Modelica.Blocks.Logical.Switch switch_OnOff_I if  with_I
    annotation (Placement(transformation(extent={{-49,-76},{-38,-87}})));
public
  Modelica.Blocks.Sources.Constant I_off_zero(k=0) if
                                                  with_I
    annotation (Placement(transformation(extent={{4.25,-4.5},{-4.25,4.5}},
                                                                     rotation=180,
        origin={-64.75,-67})));

  Modelica.Blocks.Logical.Switch switch_OnOff
    annotation (Placement(transformation(extent={{92,78},{112,58}})));
  Modelica.Blocks.Sources.RealExpression y_unlocked(y=if perUnitConversion
         then y_inactive/y_ref else y_inactive)
    annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={110,90})));
  Modelica.Blocks.Sources.BooleanExpression I_activation(y=time_lag_I_activation.y
         > Tau_lag_I)                                                   annotation (
      Placement(transformation(
        extent={{-9,7},{9,-7}},
        rotation=0,
        origin={-67,-82})));
  Exergy.XClaRa.Components.Utilities.Blocks.FirstOrderClaRa smoothPIDInput(Tau=
        Tau_in, initOption=if Tau_in > 0 then 1 else 4) annotation (
      Placement(transformation(extent={{-132,-10},{-112,10}})));
  Exergy.XClaRa.Components.Utilities.Blocks.FirstOrderClaRa smoothPIDOutput(Tau=
        Tau_out, initOption=if Tau_out > 0 then 1 else 4) annotation (
     Placement(transformation(extent={{230,-10},{250,10}})));
  Modelica.Blocks.Math.Feedback addSat
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},  rotation=180,
        origin={145,-125.5})));
  Modelica.Blocks.Logical.Timer time_lag_I_activation
    annotation (Placement(transformation(extent={{-76,175},{-63,188}})));
  Modelica.Blocks.Routing.BooleanPassThrough activate
    annotation (Placement(transformation(extent={{-39,160},{-19,180}})));
  Modelica.Blocks.Sources.BooleanExpression activate_(y=time >= t_activation)
    annotation (Placement(transformation(extent={{-154,133},{-134,153}})));
  Modelica.Blocks.Logical.And controllerActive if
                                                 use_activateInput
    annotation (Placement(transformation(extent={{-120,160},{-100,180}})));
  Modelica.Blocks.Routing.BooleanPassThrough activateIfNoSwitch if not use_activateInput
    annotation (Placement(transformation(extent={{-100,139.5},{-93,146.5}})));
  Exergy.XClaRa.Components.Utilities.Blocks.FirstOrderClaRa smoothPIDOutput1(Tau=
        Tau_add, initOption=if Tau_add > 0 then 1 else 4) annotation (
     Placement(transformation(extent={{120,57.5},{140,78}})));

  Modelica.Blocks.Sources.RealExpression y_start_I(y = y_start/y_ref)
    annotation (Placement(transformation(extent={{10,-64},{-8,-52}})));
  // The following alternative will truely allow to initialise the PID block at y_start if InitialOutput is chosen:
  // Modelica.Blocks.Sources.RealExpression y_start_I(y=if initType == InitPID.InitialOutput then y_start/y_ref - addPID.u1 - addPID.u2 elseif initType == InitPID.SteadyState then y_start/y_ref else xi_start)
  //   annotation (Placement(transformation(extent={{10,-64},{-8,-52}})));
  Modelica.Blocks.Math.Feedback resetP annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=90,
        origin={-80,110.5})));
  Modelica.Blocks.Sources.RealExpression y_unlocked1(y=resetValueP)
    annotation (Placement(transformation(extent={{10,10},{-10,-10}},
                                                                   rotation=180,
        origin={-115,111})));
  Modelica.Blocks.Math.Add resetPD annotation (Placement(transformation(
        extent={{-5,-5},{5,5}},
        rotation=0,
        origin={67,57})));
  Modelica.Blocks.Sources.RealExpression y_unlocked2(y=resetValueID) annotation (Placement(transformation(extent={{10,10},{-10,-10}},
                                                                   rotation=270,
        origin={54,30})));

initial equation
//  y_in_start = feedback.y;
//  y_out_start = fromPU.y;
//  y_aux_start = switch_OnOff.y;

equation

  assert(y_max >= y_min, "LimPID: Limits must be consistent. However, y_max (=" + String(y_max) +
                       ") < y_min (=" + String(y_min) + ")");
  if initType == InitPID.InitialOutput and (y_start < y_min or y_start > y_max) then
      Modelica.Utilities.Streams.error("LimPID: Start value y_start (=" + String(y_start) +
         ") is outside of the limits of y_min (=" + String(y_min) +") and y_max (=" + String(y_max) + ")");
  end if;
  assert(limitsAtInit or not limitsAtInit and y >= y_min and y <= y_max,
         "LimPID: During initialization the limits have been switched off.\n" +
         "After initialization, the output y (=" + String(y) +
         ") is outside of the limits of y_min (=" + String(y_min) +") and y_max (=" + String(y_max) + ")");

  when change(activate.y) then
    reinit(resetValueP, pre(toPU.y));
    reinit(resetValueID, y_unlocked.y - addPID.u2 - addPID.u3);

  end when;
  der(resetValueP)=0;
  der(resetValueID)=0;

  connect(P.y, addPID.u1) annotation (Line(points={{1,130},{1,130},{34,130},{34,64},{39,64}},
               color={0,0,127}));
  connect(D_approx.y, addPID.u2)
    annotation (Line(points={{-7,60},{-7,60},{39,60}},
                                              color={0,0,127}));
  connect(toPU.y, D_approx.u)
                       annotation (Line(points={{-85,0},{-80,0},{-80,60},{-56,60},{-30,60}},
                                                                    color={0,0,127}));
  connect(gainTrack.y, addI.u2) annotation (Line(points={{-10.65,-125.5},{-10.65,-126},{-80,-126},{-80,-108},{-66,-108}},
                           color={0,0,127}));
  connect(u_s, feedback.u1) annotation (Line(points={{-220.5,0},{-188.5,0}},
                                                                        color={0,
          0,127}));
  connect(u_m, feedback.u2) annotation (Line(points={{0,-216},{0,-160},{-180.5,-160},{-180.5,-8}},
                     color={0,0,127}));

  connect(switch_OnOff_I.y, I.u) annotation (Line(
      points={{-37.45,-81.5},{-37.45,-82},{-32,-82}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(I_off_zero.y, switch_OnOff_I.u3) annotation (Line(
      points={{-60.075,-67},{-55,-67},{-55,-77.1},{-50.1,-77.1}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(toPU.y, addI.u1) annotation (Line(
      points={{-85,0},{-80,0},{-80,-102},{-66,-102}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(addI.y, switch_OnOff_I.u1) annotation (Line(
      points={{-54.5,-105},{-54.5,-85.9},{-50.1,-85.9}},
      color={0,0,127},
      smooth=Smooth.None));

  connect(limiter.y, fromPU.u) annotation (Line(
      points={{173,68},{178,68},{178,0},{198,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(I_activation.y, switch_OnOff_I.u2) annotation (Line(
      points={{-57.1,-82},{-54,-82},{-54,-81.5},{-50.1,-81.5}},
      color={255,0,255},
      smooth=Smooth.None));
  connect(switch_OnOff.u3, y_unlocked.y) annotation (Line(
      points={{90,76},{90,90},{99,90}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(smoothPIDInput.y, toPU.u) annotation (Line(
      points={{-111,0},{-108,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(fromPU.y, smoothPIDOutput.u) annotation (Line(
      points={{221,0},{228,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(smoothPIDOutput.y, y) annotation (Line(
      points={{251,0},{270,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(addSat.y, gainTrack.u) annotation (Line(
      points={{136,-125.5},{4.3,-125.5}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(controllerActive.y, activate.u) annotation (Line(
      points={{-99,170},{-41,170}},
      color={255,0,255},
      smooth=Smooth.None));
  connect(activateInput, controllerActive.u1) annotation (Line(
      points={{-219.5,160},{-145,160},{-145,170},{-122,170}},
      color={255,0,255},
      smooth=Smooth.None));
  connect(activate_.y, controllerActive.u2) annotation (Line(
      points={{-133,143},{-125,143},{-125,162},{-122,162}},
      color={255,0,255},
      smooth=Smooth.None));
  connect(controllerActive.y, time_lag_I_activation.u) annotation (Line(
      points={{-99,170},{-94,170},{-94,181.5},{-77.3,181.5}},
      color={255,0,255},
      smooth=Smooth.None));
  connect(activate_.y, activateIfNoSwitch.u) annotation (Line(
      points={{-133,143},{-100.7,143}},
      color={255,0,255},
      smooth=Smooth.None));
  connect(activateIfNoSwitch.y, time_lag_I_activation.u) annotation (Line(
      points={{-92.65,143},{-83,143},{-83,181.5},{-77.3,181.5}},
      color={255,0,255},
      smooth=Smooth.None));
  connect(activateIfNoSwitch.y, activate.u) annotation (Line(
      points={{-92.65,143},{-41,143},{-41,170}},
      color={255,0,255},
      smooth=Smooth.None));
  connect(addPID.u2, Dzero.y) annotation (Line(
      points={{39,60},{0,60},{0,38},{-7.65,38}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(Izero.y, addPID.u3) annotation (Line(
      points={{-11,-103.75},{0,-103.75},{0,-82},{34,-82},{34,56},{39,56}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(switch_OnOff.y, smoothPIDOutput1.u) annotation (Line(
      points={{113,68},{116.85,68},{116.85,67.75},{118,67.75}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(smoothPIDOutput1.y, limiter.u) annotation (Line(
      points={{141,67.75},{150,67.75},{150,68}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(y_unlocked1.y,resetP. u2) annotation (Line(points={{-104,111},{-88,111},{-88,110.5}},
                                                                                          color={0,0,127}));
  connect(resetP.y, P.u) annotation (Line(points={{-80,119.5},{-80,130},{-22,130}},    color={0,0,127}));
  connect(toPU.y, resetP.u1) annotation (Line(points={{-85,0},{-80,0},{-80,102.5}},          color={0,0,127}));
  connect(limiter.y, addSat.u1) annotation (Line(points={{173,68},{178,68},{178,-125.5},{153,-125.5}},
                                                                                          color={0,0,127}));
  connect(y_unlocked2.y, resetPD.u2) annotation (Line(points={{54,41},{54,54},{61,54}}, color={0,0,127}));
  connect(y_start_I.y, I.y_start) annotation (Line(points={{-8.9,-58},{-20,-58},{-20,-70}}, color={0,0,127}));
  connect(activate.y, switch_OnOff.u2) annotation (Line(points={{-18,170},{75,170},{75,68},{90,68}}, color={255,0,255}));
  connect(smoothPIDOutput1.y, addSat.u2) annotation (Line(points={{141,67.75},{145,67.75},{145,-117.5}}, color={0,0,127}));
  connect(feedback.y, smoothPIDInput.u) annotation (Line(points={{-171.5,0},{-172,0},{-134,0}}, color={0,0,127}));
  connect(I.y, addPID.u3) annotation (Line(points={{-9,-82},{-9,-82},{34,-82},{34,56},{39,56}}, color={0,0,127}));
  connect(addPID.y, resetPD.u1) annotation (Line(points={{50.5,60},{56,60},{61,60}}, color={0,0,127}));
  connect(resetPD.y, switch_OnOff.u1) annotation (Line(points={{72.5,57},{74.25,57},{74.25,60},{90,60}}, color={0,0,127}));
  annotation (defaultComponentName="PID",
    Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={1,1}), graphics={
        Text(
          extent={{-100,130},{100,100}},
          lineColor={27,36,42},
          textString="%name"),  Rectangle(
          extent={{-100,-100},{100,100}},
          lineColor={221,222,223},
          fillColor={118,124,127},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-80,-90},{-80,80}},
          color={221,222,223},
          smooth=Smooth.Bezier),
        Polygon(
          points={{-80,90},{-88,68},{-72,68},{-80,90}},
          lineColor={221,222,223},
          fillColor={221,222,223},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-80,-80},{-80,40},{-80,6},{-60,-20},{30,60},{30,60},{80,60}},
          color={27,36,42},
          smooth=Smooth.Bezier),
        Text(
          extent={{-20,-20},{80,-60}},
          lineColor={221,222,223},
          fillColor={221,222,223},
          fillPattern=FillPattern.Solid,
          textString="PID"),
        Line(points={{-90,-80},{82,-80}}, color={221,222,223}),
        Polygon(
          points={{90,-80},{68,-72},{68,-88},{90,-80}},
          lineColor={221,222,223},
          fillColor={221,222,223},
          fillPattern=FillPattern.Solid)}),
    Documentation(info="<HTML>

</HTML>
", revisions="<html>
<ul>
  <li> 15.04.09 First revision, Boris Michaelsen, XRG Simulation GmbH</li>
  <li> 25.11.09 Update to independently set the gain and time constants of the PID, added a new parameter \"sign\" for case dependent control error evaluation, Friedrich Gottelt, XRG Simulation GmbH</li>
</ul>
</html>"),
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-200,-200},{260,200}},
        initialScale=0.1)));
end LimPID;
