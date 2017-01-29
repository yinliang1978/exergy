within Exergy.XClaRa.SubSystems.Boiler.Check;
model Compare_SteamGenerator_1and2
//___________________________________________________________________________//
// Component of the ClaRa library, version: 1.1.2                        //
//                                                                           //
// Licensed by the DYNCAP/DYNSTART research team under Modelica License 2.   //
// Copyright © 2013-2016, DYNCAP/DYNSTART research team.                     //
//___________________________________________________________________________//
// DYNCAP and DYNSTART are research projects supported by the German Federal //
// Ministry of Economic Affairs and Energy (FKZ 03ET2009/FKZ 03ET7060).      //
// The research team consists of the following project partners:             //
// Institute of Energy Systems (Hamburg University of Technology),           //
// Institute of Thermo-Fluid Dynamics (Hamburg University of Technology),    //
// TLK-Thermo GmbH (Braunschweig, Germany),                                  //
// XRG Simulation GmbH (Hamburg, Germany).                                   //
//___________________________________________________________________________//
  extends ClaRa.Basics.Icons.PackageIcons.ExecutableExampleb60;
  ClaRa.Components.Control.PredictorModels_3508.CoalSupplyBoiler_01_XRG
    Model_boiler(p_LS_nom=24000000)
    annotation (Placement(transformation(extent={{6,116},{48,150}})));
  Modelica.Blocks.Sources.Step ramp1(
    offset=1,
    height=-0.5,
    startTime=40000)
    annotation (Placement(transformation(extent={{-100,160},{-80,180}})));
  ClaRa.Components.Utilities.Blocks.LimPID PID(
    y_max=1,
    Tau_d=1000,
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    y_min=-1,
    perUnitConversion=true,
    y_ref=1,
    u_ref=1,
    Tau_i=500,
    initType=Modelica.Blocks.Types.InitPID.NoInit,
    k=2.5) annotation (Placement(transformation(extent={{-98,80},{-78,
            100}})));
  Modelica.Blocks.Sources.RealExpression SV_Pressure_LS(y=
        turbinesAndReheat_01_XRG.P_gen_) "Set value of live steam pressure"
    annotation (Placement(transformation(extent={{-132,80},{-112,100}})));
  Modelica.Blocks.Sources.RealExpression MV_Pressure_LS(y=homotopy(
        realPlantPower_.y, turbinesAndReheat_01_XRG.P_gen_))
    "Measurement value of live steam pressure"
    annotation (Placement(transformation(extent={{-132,64},{-112,84}})));
  Modelica.Blocks.Sources.Ramp ramp2(
    offset=1,
    height=-0.5,
    duration=600,
    startTime=6000)
    annotation (Placement(transformation(extent={{-100,120},{-80,140}})));
  Modelica.Blocks.Math.Add add
    annotation (Placement(transformation(extent={{-42,86},{-22,106}})));
  ClaRa.Components.TurboMachines.Turbines.SteamTurbineVLE_L1 HPTurbine(
    rho_nom=74.2585,
    Pi=28e5/240e5,
    p_nom=24000000) annotation (Placement(transformation(extent={{77,-74},{89,-58}})));
  ClaRa.Components.BoundaryConditions.BoundaryVLE_phxi pressureSink1(
      variable_p=false, p_const(displayUnit="Pa") = 3000) annotation (
     Placement(transformation(extent={{186,-90},{168,-70}})));
  inner ClaRa.SimCenter simCenter(redeclare replaceable
      TILMedia.VLEFluidTypes.TILMedia_InterpolatedWater fluid1)
    annotation (Placement(transformation(extent={{180,180},{200,200}})));
  ClaRa.Components.Control.PredictorModels_3508.TurbinesAndReheat_01_XRG
    turbinesAndReheat_01_XRG(CL_Deltah_p=[0.5000e7,0.9*0.1889e7;
        0.6000e7,0.9*0.1889e7; 0.8000e7,0.9*0.1910e7; 1.0000e7,0.9*
        0.1923e7; 1.2000e7,0.9*0.1930e7; 1.4000e7,0.9*0.1933e7;
        1.6000e7,0.9*0.1934e7; 1.8000e7,0.9*0.1933e7; 2.0000e7,0.9*
        0.1930e7; 2.2000e7,0.9*0.1926e7; 2.4000e7,0.9*0.1921e7;
        2.5000e7,0.9*0.1921e7]) annotation (Placement(transformation(
          extent={{82,90},{112,124}})));
  ClaRa.Components.TurboMachines.Turbines.SteamTurbineVLE_L1 LPTurbine(
    Pi=3000/28e5,
    rho_nom=8,
    p_nom=2800000) annotation (Placement(transformation(extent={{130,-80},{148,-58}})));
  Modelica.Blocks.Sources.RealExpression realPlantPower_(y=-(HPTurbine.P_t +
        LPTurbine.P_t)/turbinesAndReheat_01_XRG.P_G_nom)
    annotation (Placement(transformation(extent={{104,-26},{124,-6}})));
  Exergy.XClaRa.SubSystems.Boiler.SteamGenerator_L1 SG1(p_LS_nom=
        24000000) annotation (Placement(transformation(extent={{-6,-64},
            {20,-28}})));
  ClaRa.Components.BoundaryConditions.BoundaryVLE_hxim_flow massFlowSource_XRG(
    m_flow_const=419,
    h_const=500e3,
    variable_m_flow=true) annotation (Placement(transformation(extent=
           {{-94,-68},{-74,-48}})));
protected
  ClaRa.Basics.Interfaces.SteamSignal mediumData_b annotation (
      Placement(transformation(extent={{-134,-61},{-128,-55}})));
public
  Modelica.Blocks.Math.Gain gain(k=Model_boiler.m_flow_LS_nom) annotation (Placement(transformation(extent={{-122,-62},{-110,-50}})));
  ClaRa.Visualisation.Scope scope(
    color={255,255,0},
    hideInterface=false,
    t_end=15000) annotation (Placement(transformation(extent={{144,-44},
            {188,-4}})));
  ClaRa.Components.TurboMachines.Turbines.SteamTurbineVLE_L1 HPTurbine1(
    rho_nom=74.2585,
    Pi=28e5/240e5,
    p_nom=24000000) annotation (Placement(transformation(extent={{78,-204},{90,-188}})));
  ClaRa.Components.BoundaryConditions.BoundaryVLE_phxi pressureSink2(
      variable_p=false, p_const(displayUnit="Pa") = 3000) annotation (
     Placement(transformation(extent={{188,-220},{170,-200}})));
  ClaRa.Components.TurboMachines.Turbines.SteamTurbineVLE_L1 LPTurbine1(
    Pi=3000/28e5,
    rho_nom=8,
    p_nom=2800000) annotation (Placement(transformation(extent={{130,-210},{148,-188}})));
  Modelica.Blocks.Sources.RealExpression realPlantPower_1(
                                                         y=-(HPTurbine.P_t +
        LPTurbine.P_t)/turbinesAndReheat_01_XRG.P_G_nom)
    annotation (Placement(transformation(extent={{104,-156},{124,-136}})));
  Exergy.XClaRa.SubSystems.Boiler.SteamGenerator_L3 SG2(p_LS_nom=
        24000000) annotation (Placement(transformation(extent={{-8,-198},
            {20,-158}})));
  ClaRa.Components.BoundaryConditions.BoundaryVLE_hxim_flow massFlowSource_XRG1(
    m_flow_const=419,
    h_const=500e3,
    variable_m_flow=true) annotation (Placement(transformation(extent=
           {{-94,-198},{-74,-178}})));
  Modelica.Blocks.Math.Gain gain1(k=Model_boiler.m_flow_LS_nom) annotation (Placement(transformation(extent={{-122,-192},{-110,-180}})));
  ClaRa.Visualisation.Scope scope1(
    color={255,255,0},
    hideInterface=false,
    t_end=15000) annotation (Placement(transformation(extent={{144,-174},
            {188,-134}})));
protected
  ClaRa.Basics.Interfaces.SteamSignal mediumData_b1 annotation (
      Placement(transformation(extent={{-134,-191},{-128,-185}})));
equation

  connect(ramp1.y, Model_boiler.yT_)         annotation (Line(
      points={{-79,170},{46,170},{46,148.3},{43.968,148.3}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(SV_Pressure_LS.y, PID.u_s) annotation (Line(
      points={{-111,90},{-100,90}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(MV_Pressure_LS.y, PID.u_m)  annotation (Line(
      points={{-111,74},{-88,74},{-88,78},{-88,78}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ramp2.y, add.u1) annotation (Line(
      points={{-79,130},{-66,130},{-66,102},{-44,102}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(PID.y, add.u2) annotation (Line(
      points={{-77.1,90},{-44,90}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ramp2.y, Model_boiler.QF_setl_) annotation (Line(
      points={{-79,130},{6,130},{6,130.62}},
      color={0,0,127},
      smooth=Smooth.None));

  connect(turbinesAndReheat_01_XRG.inlet, Model_boiler.steamSignal) annotation (
     Line(
      points={{82.3,120.6},{54,120.6},{54,120.76},{48,120.76}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  connect(LPTurbine.outlet, pressureSink1.steam_a) annotation (Line(
      points={{148,-80},{168,-80}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(HPTurbine.outlet, SG1.reheat_in)                  annotation (Line(
      points={{89,-74},{89,-90},{15,-90},{15,-63.55},{14.8,-63.55}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(SG1.reheat_out, LPTurbine.inlet)                  annotation (Line(
      points={{14.8,-28},{130,-28},{130,-62.4}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(SG1.feedwater, massFlowSource_XRG.steam_a)
    annotation (Line(
      points={{7,-63.55},{-40,-63.55},{-40,-58},{-74,-58}},
      color={0,70,135},
      smooth=Smooth.None));
  connect(HPTurbine.inlet, SG1.livesteam)                  annotation (Line(
      points={{77,-61.2},{77,-18},{7,-18},{7,-28}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(gain.y, massFlowSource_XRG.m_flow) annotation (Line(
      points={{-109.4,-56},{-98,-56},{-98,-52},{-96,-52}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain.u, mediumData_b.p_) annotation (Line(
      points={{-123.2,-56},{-126,-56},{-126,-57.985},{-130.985,-57.985}},
      color={0,0,127},
      smooth=Smooth.None), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}}));
  connect(add.y, SG1.QF_setl_)                  annotation (Line(
      points={{-21,96},{-8.6,96},{-8.6,-52.75}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(realPlantPower_.y, scope.u) annotation (Line(
      points={{125,-16},{142.114,-16},{142.114,-15.0769}},
      color={0,0,127},
      smooth=Smooth.Bezier));
  connect(Model_boiler.steamSignal, mediumData_b) annotation (Line(
      points={{48,120.76},{48,186},{-140,186},{-140,-58},{-131,-58}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  connect(LPTurbine1.outlet, pressureSink2.steam_a)
                                                   annotation (Line(
      points={{148,-210},{170,-210}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(HPTurbine1.outlet, SG2.reheat_in)                 annotation (Line(
      points={{90,-204},{86,-204},{86,-216},{14,-216},{14,-197.5},{14.4,-197.5}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(SG2.reheat_out, LPTurbine1.inlet)                 annotation (Line(
      points={{14.4,-158},{130,-158},{130,-192.4}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(SG2.feedwater, massFlowSource_XRG1.steam_a)
    annotation (Line(
      points={{6,-197.5},{-40,-197.5},{-40,-188},{-74,-188}},
      color={0,70,135},
      smooth=Smooth.None));
  connect(HPTurbine1.inlet, SG2.livesteam)                 annotation (Line(
      points={{78,-191.2},{78,-146},{6,-146},{6,-158}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(gain1.y, massFlowSource_XRG1.m_flow) annotation (Line(
      points={{-109.4,-186},{-98,-186},{-98,-182},{-96,-182}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain1.u, mediumData_b1.p_)
                                   annotation (Line(
      points={{-123.2,-186},{-126,-186},{-126,-187.985},{-130.985,-187.985}},
      color={0,0,127},
      smooth=Smooth.None), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}}));
  connect(add.y, SG2.QF_setl_)                  annotation (Line(
      points={{-21,96},{-10.8,96},{-10.8,-185.5}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(realPlantPower_1.y, scope1.u)
                                      annotation (Line(
      points={{125,-146},{142.114,-146},{142.114,-145.077}},
      color={0,0,127},
      smooth=Smooth.Bezier));
  connect(Model_boiler.steamSignal, mediumData_b1)
                                                  annotation (Line(
      points={{48,120.76},{48,186},{-140,186},{-140,-188},{-131,-188}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-150,
            -250},{200,200}}),
                      graphics={
        Rectangle(extent={{-144,32},{198,-98}}, lineColor={0,0,0}),
        Rectangle(
          extent={{-144,200},{122,38}},
          lineColor={0,0,0},
          fillColor={236,236,236},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-32,50},{136,38}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="Model-based process control"),
        Text(
          extent={{48,30},{198,18}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="Process model"),
        Rectangle(extent={{-144,-98},{198,-228}},
                                                lineColor={0,0,0}),
        Text(
          extent={{-152,-222},{200,-240}},
          lineColor={0,128,0},
          lineThickness=0.5,
          fillColor={102,198,0},
          fillPattern=FillPattern.Solid,
          horizontalAlignment=TextAlignment.Left,
          textString="IDEA:
compare different boiler models
NOTE:
boiler SG1 calculates the mass flow from a characteristic line while boiler SG2 gets the mass flow rate from turbine models and calculates the pressure from a simple, 0D balance equations
 ")}),                           Icon(coordinateSystem(preserveAspectRatio=true,
          extent={{-100,-100},{100,100}})),
    experiment(StopTime=50000, Tolerance=1e-005),
    __Dymola_experimentSetupOutput);
end Compare_SteamGenerator_1and2;
