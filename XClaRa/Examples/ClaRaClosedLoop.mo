within Exergy.XClaRa.Examples;
model ClaRaClosedLoop
  "A closed steam cycle including single reheat, feedwater tank, LP and HP preheaters"
  extends ClaRa.Basics.Icons.PackageIcons.ExecutableRegressiong100;

model Regression
  extends ClaRa.Basics.Icons.RegressionSummary;

  Modelica.Blocks.Interfaces.RealInput V_flow_FWP "Volume flow of feedwater";
  Modelica.Blocks.Interfaces.RealInput p_HPT_in "Inlet pressure HP turbine";
  Modelica.Blocks.Interfaces.RealInput T_HPT_in "Inlet temperature HP turbine";
  Modelica.Blocks.Interfaces.RealInput level_FWT "level of FW tank";
  Modelica.Blocks.Interfaces.RealInput p_FWT "Pressure of FW tank";
  Modelica.Blocks.Interfaces.RealInput nom_m_flow_tapFWT;

  Real y_V_flow_FWP_min = timeExtrema_V_flow_FW.y_min;
  Real y_V_flow_FWP_max = timeExtrema_V_flow_FW.y_max;
  Real y_p_HPT_in_min = timeExtrema_p_HPTin.y_min;
  Real y_p_HPT_in_max = timeExtrema_p_HPTin.y_max;
  Real y_T_HPT_in_int = integratorT_HPTin.y;
  Real y_level_FWT_max = timeExtrema_level_FWT.y_max;
  Real y_level_FWT_min = timeExtrema_level_FWT.y_min;
  Real y_p_FWT_int = integrator_p_FWT.y;
  Real y_nom_m_flow_tapFWT = nom_m_flow_tapFWT;

  protected
  ClaRa.Components.Utilities.Blocks.TimeExtrema timeExtrema_V_flow_FW(u=
          V_flow_FWP, startTime=5000)
      annotation (Placement(transformation(extent={{-20,40},{0,60}})));
  ClaRa.Components.Utilities.Blocks.TimeExtrema timeExtrema_p_HPTin(u=p_HPT_in,
        startTime=5000)
      annotation (Placement(transformation(extent={{14,18},{34,38}})));
  ClaRa.Components.Utilities.Blocks.Integrator integratorT_HPTin(u=T_HPT_in)
      annotation (Placement(transformation(extent={{-20,-10},{0,10}})));
  ClaRa.Components.Utilities.Blocks.TimeExtrema timeExtrema_level_FWT(u=
          level_FWT, startTime=5000)
      annotation (Placement(transformation(extent={{-20,-40},{0,-20}})));
  ClaRa.Components.Utilities.Blocks.Integrator integrator_p_FWT(u=p_FWT)
      annotation (Placement(transformation(extent={{-20,-76},{0,-56}})));

end Regression;
  ClaRa.Components.TurboMachines.Turbines.SteamTurbineVLE_L1 turbine_HP1(
    p_nom=NOM.Turbine_HP.p_in,
    m_flow_nom=NOM.Turbine_HP.m_flow,
    Pi=NOM.Turbine_HP.p_out/NOM.Turbine_HP.p_in,
    rho_nom=TILMedia.VLEFluidFunctions.density_phxi(
        simCenter.fluid1,
        NOM.Turbine_HP.p_in,
        NOM.Turbine_HP.h_in),
    p_in_0=INIT.Turbine_HP.p_in,
    p_out_0=INIT.Turbine_HP.p_out,
    eta_mech=1,
    CL_eta_mflow=[0.0,NOM.efficiency_Turb_HP; 1,NOM.efficiency_Turb_HP],
    allowFlowReversal=true) annotation (Placement(transformation(extent={{-58,38},{-48,58}})));

  ClaRa.SubSystems.Boiler.SteamGenerator_L3 steamGenerator_1_XRG(
    p_LS_start=INIT.boiler.p_LS_out,
    p_RH_start=INIT.boiler.p_RS_out,
    Q_flow_F_nom=NOM.boiler.Q_nom,
    p_LS_nom=NOM.boiler.p_LS_out,
    p_RH_nom=NOM.boiler.p_RS_out,
    h_LS_nom=NOM.boiler.h_LS_out,
    h_RH_nom=NOM.boiler.h_RS_out,
    h_LS_start=INIT.boiler.h_LS_out,
    h_RH_start=INIT.boiler.h_RS_out,
    initHP=ClaRa.Basics.Choices.Init.noInit,
    initIP=ClaRa.Basics.Choices.Init.steadyState,
    CL_etaF_QF_=[0,1; 1,1],
    CL_yF_QF_=[0.4207,0.8341; 0.6246,0.8195; 0.8171,0.8049; 1,NOM.boiler.m_flow_feed*(NOM.boiler.h_LS_out - NOM.boiler.h_LS_in)/NOM.boiler.Q_nom],
    m_flow_nomLS=NOM.boiler.m_flow_nom,
    Delta_p_nomHP=NOM.Delta_p_LS_nom,
    Delta_p_nomIP=NOM.Delta_p_RS_nom,
    CL_Delta_pHP_mLS_=INIT.CharLine_Delta_p_HP_mLS_,
    CL_Delta_pIP_mLS_=INIT.CharLine_Delta_p_IP_mRS_)
                            annotation (Placement(transformation(extent={{-152,46},{-124,84}})));
  ClaRa.Components.TurboMachines.Turbines.SteamTurbineVLE_L1 Turbine_IP(
    p_in_0(displayUnit="Pa") = INIT.Turbine_IP.p_in,
    p_out_0(displayUnit="Pa") = INIT.Turbine_IP.p_out,
    p_nom=NOM.Turbine_IP.p_in,
    m_flow_nom=NOM.Turbine_IP.m_flow,
    Pi=NOM.Turbine_IP.p_out/NOM.Turbine_IP.p_in,
    rho_nom=TILMedia.VLEFluidFunctions.density_phxi(
        simCenter.fluid1,
        NOM.Turbine_IP.p_in,
        NOM.Turbine_IP.h_in),
    eta_mech=1,
    allowFlowReversal=true,
    CL_eta_mflow=[0.0,NOM.Turbine_IP.efficiency; 1,NOM.Turbine_IP.efficiency])
                            annotation (Placement(transformation(extent={{-12,38},{-2,58}})));

  ClaRa.Components.TurboMachines.Turbines.SteamTurbineVLE_L1 Turbine_LP2(
    p_out_0(displayUnit="Pa") = INIT.Turbine_LP2.p_out,
    p_in_0(displayUnit="Pa") = INIT.Turbine_LP2.p_in,
    p_nom=NOM.Turbine_LP2.p_in,
    m_flow_nom=NOM.Turbine_LP2.m_flow,
    Pi=NOM.Turbine_LP2.p_out/NOM.Turbine_LP2.p_in,
    rho_nom=TILMedia.VLEFluidFunctions.density_phxi(
        simCenter.fluid1,
        NOM.Turbine_LP2.p_in,
        NOM.Turbine_LP2.h_in),
    eta_mech=1,
    allowFlowReversal=true,
    CL_eta_mflow=[0.0,NOM.Turbine_LP2.efficiency; 1,NOM.Turbine_LP2.efficiency])
                            annotation (Placement(transformation(extent={{232,6},{242,26}})));

  Modelica.Blocks.Sources.RealExpression realPlantPower_(y=-(turbine_HP1.P_t +
        Turbine_IP.P_t + Turbine_LP1.P_t + Turbine_LP2.P_t)/550e6)
    annotation (Placement(transformation(extent={{-218,-8},{-198,12}})));
  Modelica.Blocks.Sources.Ramp PTarget(
    offset=INIT.P_target_,
    duration=1,
    startTime=10000,
    height=-0.2)
    annotation (Placement(transformation(extent={{-278,48},{-258,68}})));
  parameter Real k_PID=0.5;//1.305 "Gain of controller";
  parameter Modelica.SIunits.Time Ti_PID=650
    "Time constant of Integrator block";
                                            //216.667
  parameter Modelica.SIunits.Time startTime=2000;
  //Real Target;
  parameter Modelica.SIunits.Time Tu=127.469
    "equivalent dead time of steam generation";
                                        //127.469
  parameter Modelica.SIunits.Time Tg=204.966
    "balancing time of steam generation";
  parameter Modelica.SIunits.Time Ts=60.2459
    "Integration time of steam storage";
  ClaRa.Components.TurboMachines.Pumps.PumpVLE_L1_simple Pump_FW(eta_mech=0.8) annotation (Placement(transformation(extent={{-56,-128},{-76,-148}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=180,
        origin={272,-34})));

  Modelica.Blocks.Sources.RealExpression Q_cond(y=-(Turbine_LP2.summary.outlet.h
         - TILMedia.VLEFluidFunctions.bubbleSpecificEnthalpy_pxi(
        simCenter.fluid1, NOM.p_condenser))*Turbine_LP2.summary.outlet.m_flow)
    annotation (Placement(transformation(extent={{332,-48},{296,-20}})));
  ClaRa.Visualisation.Quadruple quadruple(decimalSpaces(p=2))
                                          annotation (Placement(transformation(
        extent={{-30,-10},{30,10}},
        rotation=0,
        origin={280,10})));
  ClaRa.Visualisation.Quadruple quadruple1
    annotation (Placement(transformation(extent={{-118,98},{-58,118}})));
  ClaRa.Visualisation.Quadruple quadruple2
    annotation (Placement(transformation(extent={{-218,98},{-158,118}})));
  ClaRa.Visualisation.Quadruple quadruple3
    annotation (Placement(transformation(extent={{12,58},{72,78}})));
  ClaRa.Visualisation.Quadruple quadruple4
    annotation (Placement(transformation(extent={{-38,98},{22,118}})));
  ClaRa.Visualisation.DynDisplay dynDisplay(
    decimalSpaces=3,
    varname="eta_el",
    unit="",
    x1=(abs(turbine_HP1.P_t + Turbine_IP.P_t + Turbine_LP1.P_t + Turbine_LP2.P_t) - abs(Pump_preheater_LP1.P_drive + Pump_cond.P_drive + Pump_FW.P_drive))/abs(steamGenerator_1_XRG.Q_flow_HP + steamGenerator_1_XRG.Q_flow_IP)) annotation (Placement(transformation(extent={{-218,26},{-186,38}})));
  ClaRa.Components.HeatExchangers.HEXvle_L3_2ph_BU condenser(
    height=5,
    width=5,
    length=10,
    initTypeShell=ClaRa.Basics.Choices.Init.steadyDensity,
    diameter_o=0.025,
    redeclare model HeatTransfer_Shell =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Constant_L3
        (alpha_nom={10000,10000}),
    redeclare model PressureLossShell =
        ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearParallelZones_L3
        (Delta_p_nom={100,100,100}),
    Tau_cond=0.3,
    Tau_evap=0.03,
    width_hotwell=4,
    length_hotwell=6,
    z_in_shell=4.9,
    z_out_shell=0.1,
    level_rel_start=0.5/6,
    m_flow_nom_shell=NOM.condenser.m_flow_in,
    p_nom_shell=NOM.condenser.p_condenser,
    p_start_shell=INIT.condenser.p_condenser,
    z_in_aux1=1,
    z_in_aux2=1)
    annotation (Placement(transformation(extent={{232,-46},{252,-22}})));

  ClaRa.Visualisation.Quadruple quadruple5(decimalSpaces(p=2))
    annotation (Placement(transformation(extent={{252,-72},{312,-52}})));

  Modelica.Blocks.Sources.RealExpression realPlantPower_1(y=-(turbine_HP1.P_t +
        Turbine_IP.P_t + Turbine_LP1.P_t + Turbine_LP2.P_t))
    annotation (Placement(transformation(extent={{-218,-28},{-198,-8}})));
  ClaRa.Components.MechanicalSeparation.FeedWaterTank_L3_advanced feedWaterTank(
    level_rel_start=0.5,
    diameter=5,
    orientation=ClaRa.Basics.Choices.GeometryOrientation.horizontal,
    p_start(displayUnit="bar") = INIT.feedwatertank.p_FWT,
    length=10,
    z_tapping=4.5,
    z_aux=4.5,
    z_vent=4.5,
    z_condensate=4.5,
    redeclare model PressureLoss =
        ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearParallelZones_L3
        (Delta_p_nom={1000,1000,1000}),
    initType=ClaRa.Basics.Choices.Init.steadyDensity,
    m_flow_cond_nom=NOM.feedwatertank.m_flow_cond,
    p_nom=NOM.feedwatertank.p_FWT,
    h_nom=NOM.feedwatertank.h_cond_in,
    m_flow_heat_nom=NOM.feedwatertank.m_flow_tap1 + NOM.feedwatertank.m_flow_tap2)
    "INIT.feedwatertank.h_cond_out"
    annotation (Placement(transformation(extent={{-12,-138},{48,-118}})));
  ClaRa.Components.TurboMachines.Pumps.PumpVLE_L1_simple Pump_cond(eta_mech=1, showExpertSummary=true) annotation (Placement(transformation(extent={{212,-112},{192,-132}})));
  Modelica.Blocks.Sources.Constant const3(k=0.5/6)
    annotation (Placement(transformation(extent={{262,-162},{242,-142}})));
  Modelica.Blocks.Sources.RealExpression condenser_relLevel(y=condenser.shell.phaseBorder.level_rel) annotation (Placement(transformation(extent={{316,-196},{276,-168}})));
  ClaRa.Components.Utilities.Blocks.LimPID PI_Pump_cond(
    sign=-1,
    y_ref=1e6,
    k=50,
    Tau_d=30,
    initType=Modelica.Blocks.Types.InitPID.InitialOutput,
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    Tau_i=300,
    y_max=NOM.Pump_cond.P_pump*10,
    y_min=NOM.Pump_cond.P_pump/200,
    y_start=INIT.Pump_cond.P_pump)
    annotation (Placement(transformation(extent={{232,-162},{212,-142}})));
  ClaRa.Visualisation.Quadruple quadruple6
    annotation (Placement(transformation(extent={{-14,-172},{46,-152}})));
  ClaRa.Components.VolumesValvesFittings.Fittings.Split_L2_Y join_IP(
    h_start=INIT.Turbine_IP.h_out,
    p_start=INIT.Turbine_IP.p_out,
    volume=0.1,
    initType=ClaRa.Basics.Choices.Init.noInit,
    p_nom=NOM.Turbine_IP.p_out,
    h_nom=NOM.Turbine_IP.h_out,
    m_flow_out_nom={NOM.feedwatertank.m_flow_cond,NOM.feedwatertank.m_flow_tap2})                                   annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={34,38})));

  ClaRa.Components.VolumesValvesFittings.Valves.ValveVLE_L1 valve_IP(redeclare
      model PressureLoss =
       ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint
        (                                                                                                    Delta_p_nom=NOM.valve_IP.Delta_p, m_flow_nom=NOM.valve_IP.m_flow))
    annotation (Placement(transformation(
        extent={{-10,6},{10,-6}},
        rotation=270,
        origin={34,-42})));
  ClaRa.Components.TurboMachines.Turbines.SteamTurbineVLE_L1 Turbine_LP1(
    p_in_0(displayUnit="Pa") = INIT.Turbine_LP1.p_in,
    p_out_0(displayUnit="Pa") = INIT.Turbine_LP1.p_out,
    p_nom=NOM.Turbine_LP1.p_in,
    m_flow_nom=NOM.Turbine_LP1.m_flow,
    Pi=NOM.Turbine_LP1.p_out/NOM.Turbine_LP1.p_in,
    rho_nom=TILMedia.VLEFluidFunctions.density_phxi(
        simCenter.fluid1,
        NOM.Turbine_LP1.p_in,
        NOM.Turbine_LP1.h_in),
    eta_mech=1,
    allowFlowReversal=true,
    CL_eta_mflow=[0.0,NOM.Turbine_LP1.efficiency; 1,NOM.Turbine_LP1.efficiency])
                            annotation (Placement(transformation(extent={{88,22},{98,42}})));

  ClaRa.Visualisation.Quadruple quadruple7
    annotation (Placement(transformation(extent={{-30,-10},{30,10}},
        rotation=0,
        origin={142,58})));
  ClaRa.Components.VolumesValvesFittings.Fittings.Split_L2_Y join_LP1(
    p_nom=INIT.Turbine_LP1.p_out,
    h_nom=INIT.Turbine_LP1.h_out,
    h_start=INIT.Turbine_LP1.h_out,
    p_start=INIT.Turbine_LP1.p_out,
    volume=0.1,
    initType=ClaRa.Basics.Choices.Init.noInit,
    m_flow_out_nom={-NOM.preheater_LP1.m_flow_cond,-NOM.preheater_LP1.m_flow_tap}) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={162,22})));

  ClaRa.Components.TurboMachines.Pumps.PumpVLE_L1_simple Pump_preheater_LP1(eta_mech=1) annotation (Placement(transformation(extent={{102,-158},{82,-178}})));
  ClaRa.Components.VolumesValvesFittings.Valves.ValveVLE_L1 valve_LP1(redeclare
      model PressureLoss =
       ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.QuadraticNominalPoint
        (                                                                                                    Delta_p_nom=NOM.valve_LP1.Delta_p_nom, m_flow_nom=NOM.valve_LP1.m_flow))
    annotation (Placement(transformation(
        extent={{-10,6},{10,-6}},
        rotation=270,
        origin={162,-42})));

  Modelica.Blocks.Sources.Constant const_reheater_LP1_relLevel(k=0.5)
    annotation (Placement(transformation(extent={{162,-202},{142,-182}})));
  Modelica.Blocks.Sources.RealExpression preheater_LP1_relLevel(y=preheater_LP1.shell.phaseBorder.level_rel) annotation (Placement(transformation(extent={{162,-236},{122,-208}})));
  ClaRa.Components.VolumesValvesFittings.Fittings.Split_L2_Y join_HP(
    volume=0.1,
    p_start=INIT.Turbine_HP.p_out,
    p_nom=NOM.Turbine_HP.p_out,
    h_nom=NOM.Turbine_HP.h_out,
    h_start=INIT.Turbine_HP.h_out,
    initType=ClaRa.Basics.Choices.Init.noInit,
    showExpertSummary=true,
    m_flow_out_nom={NOM.join_HP.m_flow_2,NOM.join_HP.m_flow_3})
                            annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={-86,18})));

  ClaRa.Components.HeatExchangers.HEXvle2vle_L3_2ph_CH_simple preheater_HP(
    redeclare replaceable model WallMaterial =
        TILMedia.SolidTypes.TILMedia_Steel,
    initTypeTubes=ClaRa.Basics.Choices.Init.noInit,
    m_flow_nom_shell=NOM.preheater_HP.m_flow_tap,
    p_nom_shell=NOM.preheater_HP.p_tap,
    h_nom_shell=NOM.preheater_HP.h_tap_out,
    m_flow_nom_tubes=NOM.preheater_HP.m_flow_cond,
    h_nom_tubes=NOM.preheater_HP.h_cond_out,
    h_start_tubes=INIT.preheater_HP.h_cond_out,
    N_passes=1,
    diameter_i=0.0189,
    diameter_o=0.0269,
    z_in_tubes=0.1,
    z_out_tubes=0.1,
    Q_flow_nom=2e8,
    z_out_shell=0.1,
    length=15,
    z_in_shell=preheater_HP.length,
    p_start_shell=INIT.preheater_HP.p_tap,
    N_tubes=1081,
    diameter=2.6,
    showExpertSummary=true,
    initTypeWall=ClaRa.Basics.Choices.Init.steadyState,
    Tau_cond=0.3,
    Tau_evap=0.03,
    alpha_ph=50000,
    redeclare model HeatTransferTubes =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.CharLine_L2
        (                                                                                                    alpha_nom=3500),
    initTypeShell=ClaRa.Basics.Choices.Init.steadyDensity,
    p_nom_tubes=NOM.preheater_HP.p_cond,
    p_start_tubes(displayUnit="bar") = INIT.preheater_HP.p_cond,
    redeclare model HeatTransfer_Shell =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Constant_L3
        (                                                                                                    alpha_nom={1650,10000}),
    redeclare model PressureLossShell =
        ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearParallelZones_L3
        (                                                                                                    Delta_p_nom={1000,1000,1000}),
    redeclare model PressureLossTubes =
        ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L2
        (                                                                                                    Delta_p_nom=10))
                    annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=90,
        origin={-140,-40})));
  ClaRa.Components.VolumesValvesFittings.Valves.ValveVLE_L1 valve1_HP(
    openingInputIsActive=false,
    showExpertSummary=true,
    redeclare model PressureLoss =
       ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.QuadraticNominalPoint
        (
        rho_in_nom=23,
        Delta_p_nom=NOM.valve1_HP.Delta_p_nom,
        m_flow_nom=NOM.valve1_HP.m_flow))
    annotation (Placement(transformation(
        extent={{10,-6},{-10,6}},
        rotation=90,
        origin={-88,-12})));
  ClaRa.Components.VolumesValvesFittings.Valves.ValveVLE_L1 valve2_HP(
    openingInputIsActive=true, redeclare model PressureLoss =
       ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.QuadraticNominalPoint
        (
        Delta_p_nom=NOM.valve2_HP.Delta_p,
        m_flow_nom=NOM.valve2_HP.m_flow,
        rho_in_nom=500))
    annotation (Placement(transformation(
        extent={{-10,-6},{10,6}},
        rotation=0,
        origin={-68,-82})));
  ClaRa.Visualisation.StatePoint_phTs statePoint_XRG
    annotation (Placement(transformation(extent={{-158,18},{-142,34}})));
  ClaRa.Components.HeatExchangers.HEXvle2vle_L3_2ph_CH_simple preheater_LP1(
    redeclare replaceable model WallMaterial =
        TILMedia.SolidTypes.TILMedia_Steel,
    m_flow_nom_shell=NOM.preheater_LP1.m_flow_tap,
    p_nom_shell=NOM.preheater_LP1.p_tap,
    h_nom_shell=NOM.preheater_LP1.h_tap_out,
    m_flow_nom_tubes=NOM.preheater_LP1.m_flow_cond,
    h_nom_tubes=NOM.preheater_LP1.h_cond_out,
    h_start_tubes=INIT.preheater_LP1.h_cond_out,
    p_start_shell=INIT.preheater_LP1.p_tap,
    diameter_i=0.017,
    diameter_o=0.020,
    N_passes=1,
    N_tubes=500,
    Q_flow_nom=2e8,
    initTypeShell=ClaRa.Basics.Choices.Init.steadyDensity,
    diameter=1.5,
    z_in_shell=preheater_LP1.length,
    z_in_tubes=preheater_LP1.diameter/2,
    z_out_tubes=preheater_LP1.diameter/2,
    length=13,
    z_out_shell=0.1,
    redeclare model HeatTransferTubes =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.CharLine_L2
        (                                                                                                    PL_alpha=[0,0.55; 0.5,0.65; 0.7,0.72; 0.8,0.77; 1,1], alpha_nom=3000),
    redeclare model PressureLossTubes =
        ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L2
        (                                                                                                    Delta_p_nom=1000),
    T_w_start={300,320,340},
    initTypeWall=ClaRa.Basics.Choices.Init.steadyState,
    Tau_cond=0.3,
    Tau_evap=0.03,
    redeclare model HeatTransfer_Shell =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Constant_L3
        (                                                                                                    alpha_nom={1500,8000}),
    redeclare model PressureLossShell =
        ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearParallelZones_L3
        (                                                                                                    Delta_p_nom={100,100,100}),
    p_nom_tubes=NOM.preheater_LP1.p_cond,
    p_start_tubes(displayUnit="bar") = INIT.preheater_LP1.p_cond,
    initTypeTubes=ClaRa.Basics.Choices.Init.noInit)
                   annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=180,
        origin={162,-124})));

  Modelica.Blocks.Sources.Constant const2(k=0.5)
    annotation (Placement(transformation(extent={{-208,-102},{-188,-82}})));
  Modelica.Blocks.Sources.RealExpression preheater_HP_relLevel2(y=preheater_HP.shell.phaseBorder.level_rel) annotation (Placement(transformation(extent={{-238,-76},{-198,-48}})));
  Modelica.Blocks.Continuous.LimPID PI_Pump_preheater2(
    Td=1,
    yMax=1,
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=100,
    initType=Modelica.Blocks.Types.InitPID.InitialOutput,
    y_start=1,
    yMin=0.01,
    Ti=10000)
    annotation (Placement(transformation(extent={{-188,-72},{-168,-52}})));

  Modelica.Blocks.Continuous.FirstOrder measurement(
    initType=Modelica.Blocks.Types.Init.InitialOutput,
    y_start=0.1,
    T=10)
    annotation (Placement(transformation(extent={{262,-192},{242,-172}})));
  ClaRa.StaticCycles.Check.StaticCycleExamples.InitSteamCycle_T_4_Pr_F1_C1_version2
    INIT(
    P_target_=1,
    p_condenser=NOM.p_condenser,
    preheater_HP_p_tap=NOM.preheater_HP_p_tap,
    preheater_HP_m_flow_tap=NOM.preheater_HP_m_flow_tap,
    preheater_LP1_p_tap=NOM.preheater_LP1_p_tap,
    preheater_LP1_m_flow_tap=NOM.preheater_LP1_m_flow_tap,
    p_FWT=NOM.p_FWT,
    valve1_HP_Delta_p_nom=NOM.valve1_HP_Delta_p_nom,
    valve_LP1_Delta_p_nom=NOM.valve_LP1_Delta_p_nom,
    valve_LP2_Delta_p_nom=NOM.valve_LP2_Delta_p_nom,
    T_LS_nom=NOM.T_LS_nom,
    T_RS_nom=NOM.T_RS_nom,
    p_LS_out_nom=NOM.p_LS_out_nom,
    p_RS_out_nom=NOM.p_RS_out_nom,
    Delta_p_LS_nom=NOM.Delta_p_LS_nom,
    Delta_p_RS_nom=NOM.Delta_p_RS_nom,
    CharLine_Delta_p_HP_mLS_=NOM.CharLine_Delta_p_HP_mLS_,
    CharLine_Delta_p_IP_mRS_=NOM.CharLine_Delta_p_IP_mRS_,
    efficiency_Pump_cond=NOM.efficiency_Pump_cond,
    efficiency_Pump_preheater_LP1=NOM.efficiency_Pump_preheater_LP1,
    efficiency_Pump_FW=NOM.efficiency_Pump_FW,
    tapping_IP_pressure=NOM.tapping_IP_pressure,
    efficiency_Turb_HP=NOM.efficiency_Turb_HP,
    efficiency_Turb_IP=NOM.efficiency_Turb_IP,
    efficiency_Turb_LP1=NOM.efficiency_Turb_LP1,
    efficiency_Turb_LP2=NOM.efficiency_Turb_LP2,
    m_flow_nom=NOM.m_flow_nom)
    annotation (Placement(transformation(extent={{-312,-236},{-292,-216}})));
  ClaRa.Components.VolumesValvesFittings.Valves.ValveVLE_L1 valve_LP2(
                                          Tau=1e-3, redeclare model
      PressureLoss =
       ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint
        (                                                                                                    m_flow_nom=NOM.valve_LP2.m_flow, Delta_p_nom=NOM.valve_LP2.Delta_p_nom))
    annotation (Placement(transformation(
        extent={{-10,-6},{10,6}},
        rotation=180,
        origin={102,-122})));
  ClaRa.Components.VolumesValvesFittings.Fittings.Join_L2_Y join_LP_main(
    initType=ClaRa.Basics.Choices.Init.noInit,
    useHomotopy=false,
    volume=0.2,
    m_flow_in_nom={NOM.join_LP_main.m_flow_1,NOM.join_LP_main.m_flow_2},
    p_nom=NOM.join_LP_main.p,
    h_nom=NOM.join_LP_main.h3,
    h_start=INIT.join_LP_main.h3,
    p_start=INIT.join_LP_main.p)
                annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={62,-122})));

  ClaRa.Components.Utilities.Blocks.LimPID PI_preheater1(
    sign=-1,
    Tau_d=30,
    initType=Modelica.Blocks.Types.InitPID.InitialOutput,
    y_max=NOM.Pump_preheater_LP1.P_pump*1.5,
    y_min=NOM.Pump_preheater_LP1.P_pump/100,
    y_ref=1e5,
    y_start=INIT.Pump_preheater_LP1.P_pump,
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=200,
    Tau_i=200)
    annotation (Placement(transformation(extent={{122,-202},{102,-182}})));
  ClaRa.Visualisation.Quadruple quadruple8
    annotation (Placement(transformation(extent={{-27,-8},{27,8}},
        rotation=0,
        origin={71,-62})));
  ClaRa.Visualisation.Quadruple quadruple9
    annotation (Placement(transformation(extent={{-30,-10},{30,10}},
        rotation=0,
        origin={202,-62})));
  ClaRa.Visualisation.Quadruple quadruple10
    annotation (Placement(transformation(extent={{-30,-10},{30,10}},
        rotation=0,
        origin={102,-102})));
  ClaRa.Visualisation.Quadruple quadruple11
    annotation (Placement(transformation(extent={{-30,-10},{30,10}},
        rotation=0,
        origin={-50,-32})));
  ClaRa.Visualisation.Quadruple quadruple12
    annotation (Placement(transformation(extent={{-30,-10},{30,10}},
        rotation=0,
        origin={-68,-104})));
  ClaRa.Visualisation.Quadruple quadruple13
    annotation (Placement(transformation(extent={{-30,-10},{30,10}},
        rotation=0,
        origin={202,-92})));
  ClaRa.Visualisation.Quadruple quadruple14
    annotation (Placement(transformation(extent={{-26,-8},{26,8}},
        rotation=0,
        origin={102,-142})));
  ClaRa.Visualisation.StatePoint_phTs statePoint_XRG2
    annotation (Placement(transformation(extent={{142,-162},{158,-146}})));
  ClaRa.Visualisation.StatePoint_phTs statePoint_XRG3
    annotation (Placement(transformation(extent={{142,-118},{158,-102}})));
  ClaRa.Visualisation.StatePoint_phTs statePoint_XRG1
    annotation (Placement(transformation(extent={{-168,-32},{-152,-16}})));
  Modelica.Blocks.Math.Gain Nominal_PowerFeedwaterPump1(k=NOM.Pump_FW.P_pump)
    annotation (Placement(transformation(extent={{-118,-172},{-98,-152}})));
  ClaRa.Visualisation.DynDisplay fuel1(
    unit="p.u.",
    decimalSpaces=1,
    x1=Nominal_PowerFeedwaterPump1.u,
    varname="feedwater")
    annotation (Placement(transformation(extent={{-124,-198},{-92,-186}})));
  ClaRa.Visualisation.DynDisplay fuel2(
    x1=time/3600,
    unit="h",
    decimalSpaces=2,
    varname="Time")
    annotation (Placement(transformation(extent={{-218,38},{-186,50}})));
  ClaRa.Visualisation.DynDisplay fuel3(
    decimalSpaces=2,
    varname="condenser level",
    unit="m",
    x1=condenser.shell.summary.outline.level_abs)
    annotation (Placement(transformation(extent={{188,-38},{228,-26}})));
  ClaRa.StaticCycles.Check.StaticCycleExamples.InitSteamCycle_T_4_Pr_F1_C1_version2
    NOM(
    final P_target_=1,
    preheater_HP_p_tap=51.95e5,
    Delta_p_RS_nom=4.91e5,
    tapping_IP_pressure=13e5,
    efficiency_Pump_cond=0.9,
    efficiency_Pump_preheater_LP1=0.8,
    efficiency_Pump_FW=0.9,
    efficiency_Turb_HP=0.93,
    efficiency_Turb_IP=0.93,
    efficiency_Turb_LP1=0.94,
    efficiency_Turb_LP2=0.94)
    annotation (Placement(transformation(extent={{-306,-200},{-286,-180}})));
  inner ClaRa.SimCenter simCenter(contributeToCycleSummary=true, redeclare
      TILMedia.VLEFluidTypes.TILMedia_InterpolatedWater fluid1)
    annotation (Placement(transformation(extent={{-280,-220},{-240,-200}})));
  Regression regression(
    V_flow_FWP = Pump_FW.summary.outline.V_flow,
    p_HPT_in = turbine_HP1.summary.inlet.p,
    T_HPT_in = turbine_HP1.summary.inlet.T,
    level_FWT = feedWaterTank.volume.summary.outline.level_abs,
    p_FWT = feedWaterTank.volume.summary.fluid.p[2],
    nom_m_flow_tapFWT = NOM.feedwatertank.tap_in2.m_flow) annotation (Placement(transformation(extent={{-320,-140},{-300,-120}})));
equation
  connect(steamGenerator_1_XRG.reheat_out, Turbine_IP.inlet)
                                                            annotation (Line(
      points={{-129.6,84},{-12,84},{-12,54}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));

  connect(Pump_FW.inlet, feedWaterTank.outlet)
                                            annotation (Line(
      points={{-56,-138},{-8,-138}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(turbine_HP1.inlet, steamGenerator_1_XRG.livesteam)
                                                           annotation (Line(
      points={{-58,54},{-58,90},{-138,90},{-138,84}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));

  connect(Pump_FW.outlet, preheater_HP.In2) annotation (Line(
      points={{-76,-138},{-138,-138},{-138,-50}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(preheater_HP.Out2, steamGenerator_1_XRG.feedwater) annotation (Line(
      points={{-138,-30},{-138,46.475}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(valve1_HP.outlet, preheater_HP.In1)
                                             annotation (Line(
      points={{-88,-22},{-88,-40},{-130.2,-40}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(preheater_HP.Out2, statePoint_XRG.port) annotation (Line(
      points={{-138,-30},{-138,18},{-158,18}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(preheater_LP1.In1, valve_LP1.outlet) annotation (Line(
      points={{162,-114.2},{162,-52}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(preheater_LP1.In2, Pump_cond.outlet) annotation (Line(
      points={{172,-122},{192,-122}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(preheater_LP1.Out1, Pump_preheater_LP1.inlet) annotation (Line(
      points={{162,-134},{162,-168},{102,-168}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(PI_Pump_preheater2.y,valve2_HP. opening_in) annotation (Line(
      points={{-167,-62},{-68,-62},{-68,-73}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const2.y, PI_Pump_preheater2.u_m) annotation (Line(
      points={{-187,-92},{-178,-92},{-178,-74}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const3.y, PI_Pump_cond.u_s) annotation (Line(
      points={{241,-152},{234,-152}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(measurement.u, condenser_relLevel.y) annotation (Line(
      points={{264,-182},{274,-182}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(measurement.y, PI_Pump_cond.u_m) annotation (Line(
      points={{241,-182},{221.9,-182},{221.9,-164}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(join_LP1.inlet, Turbine_LP1.outlet) annotation (Line(
      points={{152,22},{98,22}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(join_LP1.outlet1, Turbine_LP2.inlet) annotation (Line(
      points={{172,22},{232,22}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(join_LP1.outlet2, valve_LP1.inlet) annotation (Line(
      points={{162,12},{162,-32}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(Turbine_IP.outlet, join_IP.inlet) annotation (Line(
      points={{-2,38},{24,38}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(join_IP.outlet1, Turbine_LP1.inlet) annotation (Line(
      points={{44,38},{88,38}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(join_IP.outlet2, valve_IP.inlet) annotation (Line(
      points={{34,28},{34,-32}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(join_HP.inlet, turbine_HP1.outlet) annotation (Line(
      points={{-76,18},{-38,18},{-38,38},{-48,38}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(join_HP.outlet1, steamGenerator_1_XRG.reheat_in) annotation (Line(
      points={{-96,18},{-129.6,18},{-129.6,46.475}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(join_HP.outlet2, valve1_HP.inlet) annotation (Line(
      points={{-86,8},{-86,4},{-88,4},{-88,-2}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(valve_LP2.inlet, preheater_LP1.Out2) annotation (Line(
      points={{112,-122},{152,-122}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(valve_LP2.outlet, join_LP_main.inlet1) annotation (Line(
      points={{92,-122},{72,-122}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(PI_Pump_cond.y, Pump_cond.P_drive) annotation (Line(
      points={{211,-152},{202,-152},{202,-134}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(PI_preheater1.u_s, const_reheater_LP1_relLevel.y) annotation (Line(
      points={{124,-192},{141,-192}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(preheater_LP1_relLevel.y, PI_preheater1.u_m) annotation (Line(
      points={{120,-222},{111.8,-222},{111.8,-204},{111.9,-204}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(PI_preheater1.y, Pump_preheater_LP1.P_drive) annotation (Line(
      points={{101,-192},{92,-192},{92,-180}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(steamGenerator_1_XRG.eye_LS, quadruple2.eye) annotation (Line(
      points={{-143.6,84.7125},{-143.6,108},{-218,108}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(steamGenerator_1_XRG.eye_RH, quadruple1.eye) annotation (Line(
      points={{-124,84.7125},{-124,108},{-118,108}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(turbine_HP1.eye, quadruple4.eye) annotation (Line(
      points={{-47,42},{-38,42},{-38,108}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(Turbine_IP.eye, quadruple3.eye) annotation (Line(
      points={{-1,42},{2,42},{2,68},{12,68}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(Turbine_LP1.eye, quadruple7.eye) annotation (Line(
      points={{99,26},{102,26},{102,58},{112,58}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(Turbine_LP2.eye, quadruple.eye) annotation (Line(
      points={{243,10},{250,10}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(quadruple6.eye, feedWaterTank.eye) annotation (Line(
      points={{-14,-162},{-14,-150},{-4,-150},{-4,-139}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(quadruple8.eye, valve_IP.eye) annotation (Line(
      points={{44,-62},{38,-62},{38,-52}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(quadruple9.eye, valve_LP1.eye) annotation (Line(
      points={{172,-62},{166,-62},{166,-52}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(quadruple10.eye, valve_LP2.eye) annotation (Line(
      points={{72,-102},{68,-102},{68,-114},{92,-114},{92,-118}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(valve1_HP.eye, quadruple11.eye) annotation (Line(
      points={{-84,-22},{-84,-21.7},{-80,-21.7},{-80,-32}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(quadruple12.eye, valve2_HP.eye) annotation (Line(
      points={{-98,-104},{-32,-104},{-32,-86},{-58,-86}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(quadruple13.eye, Pump_cond.eye) annotation (Line(
      points={{172,-92},{172,-116},{191,-116}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(quadruple14.eye, Pump_preheater_LP1.eye) annotation (Line(
      points={{76,-142},{70,-142},{70,-162},{81,-162}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(statePoint_XRG2.port, Pump_preheater_LP1.inlet) annotation (Line(
      points={{142,-162},{142,-168},{102,-168}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(statePoint_XRG3.port, preheater_LP1.Out2) annotation (Line(
      points={{142,-118},{142,-122},{152,-122}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(preheater_HP.Out1, valve2_HP.inlet) annotation (Line(
      points={{-150,-40},{-158,-40},{-158,-82},{-78,-82}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(statePoint_XRG1.port, preheater_HP.Out1) annotation (Line(
      points={{-168,-32},{-160,-32},{-160,-40},{-150,-40}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));

  connect(Nominal_PowerFeedwaterPump1.y, Pump_FW.P_drive) annotation (Line(
      points={{-97,-162},{-66,-162},{-66,-150}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(valve_IP.outlet, feedWaterTank.aux) annotation (Line(
      points={{34,-52},{34,-122}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));
  connect(feedWaterTank.heatingSteam, valve2_HP.outlet) annotation (Line(
      points={{-2,-120},{-2,-82},{-58,-82}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(Pump_preheater_LP1.outlet, join_LP_main.inlet2) annotation (Line(
      points={{82,-168},{62,-168},{62,-132}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));
  connect(join_LP_main.outlet, feedWaterTank.condensate) annotation (Line(
      points={{52,-122},{38,-122}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));
  connect(PTarget.y, steamGenerator_1_XRG.QF_setl_) annotation (Line(
      points={{-257,58},{-210,58},{-210,57.875},{-154.8,57.875}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(PTarget.y, Nominal_PowerFeedwaterPump1.u) annotation (Line(
      points={{-257,58},{-248,58},{-248,-162},{-120,-162}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(condenser.outlet, Pump_cond.inlet)
    annotation (Line(
      points={{242,-46},{242,-122},{212,-122}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));
  connect(Turbine_LP2.outlet, condenser.inlet)
    annotation (Line(
      points={{242,6},{242,-22}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));
  connect(condenser.eye, quadruple5.eye) annotation (Line(
      points={{246,-47.2},{246,-62},{252,-62}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(condenser.heat, prescribedHeatFlow.port) annotation (Line(
      points={{252,-34},{262,-34}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(Q_cond.y, prescribedHeatFlow.Q_flow) annotation (Line(
      points={{294.2,-34},{282,-34}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(preheater_HP_relLevel2.y, PI_Pump_preheater2.u_s) annotation (Line(
      points={{-196,-62},{-190,-62}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-320,-240},
            {340,240}}),
                      graphics={
        Rectangle(
          extent={{-282,-160},{-240,-240}},
          lineColor={175,175,175},
          fillColor={215,215,215},
          fillPattern=FillPattern.Sphere),
        Rectangle(
          extent={{-322,-160},{-282,-240}},
          lineColor={175,175,175},
          fillColor={215,215,215},
          fillPattern=FillPattern.Sphere),
        Text(
          extent={{-316,-158},{-288,-178}},
          lineColor={75,75,75},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid,
          textStyle={TextStyle.Bold},
          textString="CycleINIT"),
        Text(
          extent={{-322,-198},{-282,-218}},
          lineColor={75,75,75},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid,
          textStyle={TextStyle.Bold},
          textString="CycleSettings"),
        Text(
          extent={{-278,-178},{-244,-196}},
          lineColor={75,75,75},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid,
          textStyle={TextStyle.Bold},
          textString="Model
Properties"),                     Text(
          extent={{-300,210},{-102,170}},
          lineColor={115,150,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=10,
          textString="______________________________________________________________________________________________
PURPOSE:
Example of a steam cycle with single reheat, feedwatertank  and feed water preheaters.
Initial and nominal values are computed with help of the static cycle package.
Models provide information for automatic efficiency calculation within SimCenter-model.
______________________________________________________________________________________________
"),                    Text(
          extent={{-312,234},{2,214}},
          lineColor={115,150,0},
          fontSize=31,
          textString="TESTED -- 2015-01-26 //TT"),Text(
          extent={{-300,138},{-136,124}},
          lineColor={115,150,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=8,
          textString="______________________________________________________________________________________________________________
Remarks: 
______________________________________________________________________________________________________________
"),                   Text(
          extent={{-300,162},{-100,144}},
          lineColor={115,150,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=10,
          textString="______________________________________________________________________________________________
Scenario:  
Reduction of the target output power of the burner by 20 percent after 10000s.
______________________________________________________________________________________________
"),     Rectangle(
          extent={{-320,240},{340,-240}},
          lineColor={115,150,0},
          lineThickness=0.5)}),  Icon(coordinateSystem(preserveAspectRatio=true,
          extent={{-100,-100},{100,100}})),
    experiment(
      StopTime=20000,
      __Dymola_NumberOfIntervals=5000,
      Tolerance=1e-006,
      __Dymola_Algorithm="Dassl"),
    __Dymola_experimentSetupOutput);
end ClaRaClosedLoop;
