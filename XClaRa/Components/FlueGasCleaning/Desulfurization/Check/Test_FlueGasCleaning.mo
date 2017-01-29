within Exergy.XClaRa.Components.FlueGasCleaning.Desulfurization.Check;
model Test_FlueGasCleaning
  extends ClaRa.Basics.Icons.PackageIcons.ExecutableRegressiong100;
model Regression
 extends ClaRa.Basics.Icons.RegressionSummary;

 Modelica.Blocks.Interfaces.RealInput deSO_ideal_L1_1_m_flow_out
      "Outlet mass flow of deSO_ideal_L1_1";
 Modelica.Blocks.Interfaces.RealInput deSO_ideal_L1_1_T_out
      "Outlet temperature of deSO_ideal_L1_1";
 Modelica.Blocks.Interfaces.RealInput deNOx_m_flow_out
      "Outlet mass flow of deNOx";
 Modelica.Blocks.Interfaces.RealInput deNOx_T_out "Outlet temperature of deNOx";
 Modelica.Blocks.Interfaces.RealInput deNOx_NH3_port_m_flow_out
      "Outlet mass flow of deNOx NH3 port";
 Modelica.Blocks.Interfaces.RealInput deNOx_NH3_port_T_out
      "Outlet temperature of deNOx NH3 port";
 Modelica.Blocks.Interfaces.RealInput e_Filter_simple_m_flow_out
      "Outlet mass flow of simple e-filter";
 Modelica.Blocks.Interfaces.RealInput e_Filter_simple_T_out
      "Outlet temperature of  simple e-filter";
 Modelica.Blocks.Interfaces.RealInput e_Filter_empirical_m_flow_out
      "Outlet mass flow of  empirical e-filter";
 Modelica.Blocks.Interfaces.RealInput e_Filter_empirical_T_out
      "Outlet temperature of  empirical e-filter";
 Modelica.Blocks.Interfaces.RealInput e_Filter_detailed_m_flow_out
      "Outlet mass flow of  detailed e-filter";
 Modelica.Blocks.Interfaces.RealInput e_Filter_detailed_T_out
      "Outlet temperature of  detailed e-filter";
 Modelica.Blocks.Interfaces.RealInput pipeFlowGas_L4_Simple_T_out
      "Outlet mass flow of flue gas pipe";
 Modelica.Blocks.Interfaces.RealInput pipeFlowGas_L4_Simple_p_out
      "Outlet temperature of flue gas pipe";

 Real y_deSO_ideal_L1_1_m_flow_out = integrator_deSO_ideal_L1_1_m_flow_out.y;
 Real y_deSO_ideal_L1_1_T_out = integrator_deSO_ideal_L1_1_T_out.y;
 Real y_deNOx_m_flow_out = integrator_deNOx_m_flow_out.y;
 Real y_deNOx_T_out = integrator_deNOx_T_out.y;
 Real y_deNOx_NH3_port_m_flow_out = integrator_deNOx_NH3_port_m_flow_out.y;
 Real y_deNOx_NH3_port_T_out = integrator_deNOx_NH3_port_T_out.y;
 Real y_e_Filter_simple_m_flow_out = integrator_e_Filter_simple_m_flow_out.y;
 Real y_e_Filter_simple_T_out = integrator_e_Filter_simple_T_out.y;
 Real y_e_Filter_empirical_m_flow_out = integrator_e_Filter_empirical_m_flow_out.y;
 Real y_e_Filter_empirical_T_out = integrator_e_Filter_empirical_T_out.y;
 Real y_e_Filter_detailed_m_flow_out = integrator_e_Filter_detailed_m_flow_out.y;
 Real y_e_Filter_detailed_T_out = integrator_e_Filter_detailed_T_out.y;
 Real y_pipeFlowGas_L4_Simple_T_out = integrator_pipeFlowGas_L4_Simple_T_out.y;
 Real y_pipeFlowGas_L4_Simple_p_out = integrator_pipeFlowGas_L4_Simple_p_out.y;

  protected
 Utilities.Blocks.Integrator integrator_deSO_ideal_L1_1_m_flow_out(u=
          deSO_ideal_L1_1_m_flow_out) annotation (Placement(
          transformation(extent={{-80,62},{-60,82}})));
 Utilities.Blocks.Integrator integrator_deSO_ideal_L1_1_T_out(u=
          deSO_ideal_L1_1_T_out) annotation (Placement(
          transformation(extent={{-80,-34},{-60,-14}})));
 Utilities.Blocks.Integrator integrator_deNOx_m_flow_out(u=
          deNOx_m_flow_out) annotation (Placement(transformation(
            extent={{-80,-2},{-60,18}})));
 Utilities.Blocks.Integrator integrator_deNOx_T_out(u=deNOx_T_out)
      annotation (Placement(transformation(extent={{-80,30},{-60,50}})));
 Utilities.Blocks.Integrator integrator_deNOx_NH3_port_m_flow_out(u=
          deNOx_NH3_port_m_flow_out) annotation (Placement(
          transformation(extent={{-80,-66},{-60,-46}})));
 Utilities.Blocks.Integrator integrator_deNOx_NH3_port_T_out(u=
          deNOx_NH3_port_T_out) annotation (Placement(
          transformation(extent={{-80,-98},{-60,-78}})));
 Utilities.Blocks.Integrator integrator_e_Filter_simple_m_flow_out(u=
          e_Filter_simple_m_flow_out) annotation (Placement(
          transformation(extent={{-18,62},{2,82}})));
 Utilities.Blocks.Integrator integrator_e_Filter_simple_T_out(u=
          e_Filter_simple_T_out) annotation (Placement(
          transformation(extent={{-18,30},{2,50}})));
 Utilities.Blocks.Integrator integrator_e_Filter_empirical_m_flow_out(u=
          e_Filter_empirical_m_flow_out) annotation (Placement(
          transformation(extent={{-18,-66},{2,-46}})));
 Utilities.Blocks.Integrator integrator_e_Filter_empirical_T_out(u=
          e_Filter_empirical_T_out) annotation (Placement(
          transformation(extent={{-18,-98},{2,-78}})));
 Utilities.Blocks.Integrator integrator_e_Filter_detailed_m_flow_out(u=
          e_Filter_detailed_m_flow_out) annotation (Placement(
          transformation(extent={{-18,-2},{2,18}})));
 Utilities.Blocks.Integrator integrator_e_Filter_detailed_T_out(u=
          e_Filter_detailed_T_out) annotation (Placement(
          transformation(extent={{-18,-34},{2,-14}})));
 Utilities.Blocks.Integrator integrator_pipeFlowGas_L4_Simple_T_out(u=
          pipeFlowGas_L4_Simple_T_out) annotation (Placement(
          transformation(extent={{-18,-2},{2,18}})));
 Utilities.Blocks.Integrator integrator_pipeFlowGas_L4_Simple_p_out(u=
          pipeFlowGas_L4_Simple_p_out) annotation (Placement(
          transformation(extent={{-18,-34},{2,-14}})));

end Regression;

  Desulfurization_L1_ideal deSO_ideal_L1_1(
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.HollowCylinder,
    m_flow_nom=530,
    xi_start={0,0,0.21,0.00099,0.7,0.0393,0,0.0367,0},
    SOx_separationRate=0.95,
    initType=ClaRa.Basics.Choices.Init.noInit,
    T_start=395.843,
    p_start=101800,
    useStabilisedMassFlow=false,
    redeclare model PressureLoss =
        Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L2)
    annotation (Placement(transformation(extent={{-38,26},{-18,46}})));

  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T(
    m_flow_const=551.153,
    xi_const={0,0,0.21,0.00099,0.7,0.0393,0,0.0367,0},
    T_const=395.843) annotation (Placement(transformation(extent={{-78,26},{-58,46}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT(
    xi_const={0,0,0.21,0.00099,0.7,0.0393,0,0.0367,0},
    p_const=101800,
    T_const=293.15) annotation (Placement(transformation(extent={{18,26},{-2,46}})));
  inner ClaRa.SimCenter simCenter(
    contributeToCycleSummary=true,
    redeclare TILMedia.GasTypes.FlueGasTILMedia flueGasModel,
    showExpertSummary=true) annotation (Placement(transformation(
          extent={{68,68},{88,88}})));
  Denitrification.Denitrification_L1
                     deNOx(
    separationRate=0.9,
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.GenericGeometry (
          volume=5),
    useHomotopy=simCenter.useHomotopy,
    use_dynamicMassbalance=true,
    redeclare model HeatTransfer =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Adiabat_L2,
    redeclare model PressureLoss =
        Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L2
        (                                                                                                    Delta_p_nom=100))
    annotation (Placement(transformation(extent={{4,-58},{24,-38}})));
  BoundaryConditions.BoundaryGas_Txim_flow                  idealGasFlowSource_XRG2(
    m_flow_const=10,
    variable_m_flow=true,
    variable_T=true,
    xi_const={0.01,0.01,0.73,0.01,0.065,0.036,0.01,0.13,0.0})
                                                 annotation (Placement(transformation(extent={{-30,-58},{-10,-38}})));
  Modelica.Blocks.Sources.Ramp massFlowRate2(
    offset=1e-3,
    height=1e-3,
    duration=10,
    startTime=100)
                 annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-68,-30})));
  Modelica.Blocks.Sources.Ramp Temperature2(
    duration=1,
    height=25,
    offset=273.15 + 200,
    startTime=150)      annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-68,-62})));
  BoundaryConditions.BoundaryGas_pTxi                  idealGasPressureSink_XRG1(p_const=100000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={52,-48})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperatureTop1(T=293.15)
                annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-68,0})));
  BoundaryConditions.BoundaryGas_Txim_flow                  idealGasFlowSource_XRG1(
    m_flow_const=10,
    variable_m_flow=true,
    variable_T=true,
    medium=simCenter.flueGasModel,
    xi_const={0.01,0.01,0.73,0.01,0.065,0.036,0.01,0.13,0.0})
                                                    annotation (Placement(transformation(extent={{-32,-184},{-12,-164}})));
  BoundaryConditions.BoundaryGas_pTxi                  idealGasPressureSink(medium=simCenter.flueGasModel, p_const=100000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={52,-174})));
  Modelica.Blocks.Sources.Ramp massFlowRate1(
    offset=1e-3,
    height=9e-3,
    startTime=5,
    duration=5)  annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-70,-156})));
  Modelica.Blocks.Sources.Ramp Temperature1(
    duration=1,
    height=20,
    startTime=1,
    offset=273.15 + 250)
                        annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-70,-188})));
  Denitrification.Denitrification_L1_NH3port deNOx_NH3_port(redeclare model
      HeatTransfer =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Adiabat_L2,                                                                                    redeclare
      model PressureLoss =
        Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L2
        (                                                                                                    Delta_p_nom=100)) annotation (Placement(transformation(extent={{4,-184},{24,-164}})));
  BoundaryConditions.BoundaryGas_Txim_flow                  idealGasFlowSource_XRG3(
    variable_T=false,
    m_flow_const=0.0010,
    variable_m_flow=true,
    medium=simCenter.flueGasModel,
    T_const=523.15,
    xi_const={0/9,0/9,0/9,0/9,0/9,0/9,0/9,0/9,9/9}) annotation (Placement(transformation(extent={{-32,-152},{-12,-132}})));
  Modelica.Blocks.Sources.Ramp massFlowRate3(
    offset=0.5e-4,
    height=5e-3,
    startTime=5,
    duration=30) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-70,-124})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperatureTop2(T=293.15)
                annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-70,-94})));
  BoundaryConditions.BoundaryGas_Txim_flow                  idealGasFlowSource_XRG(
    m_flow_const=10,
    variable_m_flow=true,
    variable_T=true,
    xi_const={0.01,0,0.73,0,0.065,0.036,0,0.13,0.0}) annotation (Placement(transformation(extent={{180,12},{200,32}})));
  BoundaryConditions.BoundaryGas_pTxi                  idealGasPressureSink1(
                                                                            p_const=100000, xi_const={0.0,0,0.73,0,0.065,0.036,0,0.13,0.0}) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={290,22})));
  Modelica.Blocks.Sources.Ramp massFlowRate(
    startTime=5,
    duration=1,
    height=-2,
    offset=1)    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={150,42})));
  Modelica.Blocks.Sources.Ramp Temperature(
    duration=1,
    startTime=1,
    height=50,
    offset=273.15 + 150)
                        annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={150,10})));
  E_Filter.E_Filter_L2_simple                                  e_Filter_dynamic(separationRate=0.9995, xi_start={0.0,0,0.73,0,0.065,0.036,0,0.13,0.0},
    redeclare model PressureLoss =
        Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.QuadraticNominalPoint_L2
        (                                                                                                    Delta_p_nom=100))              annotation (Placement(transformation(extent={{216,12},{236,32}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperatureTop3(T=293.15)
                annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={150,70})));
  BoundaryConditions.BoundaryGas_Txim_flow                  idealGasFlowSource_XRG4(
    m_flow_const=10,
    variable_m_flow=true,
    variable_T=true,
    xi_const={1/9,1/9,1/9,1/9,1/9,1/9,1/9,1/9,1/9})
                                                 annotation (Placement(transformation(extent={{180,-84},{200,-64}})));
  BoundaryConditions.BoundaryGas_pTxi                  idealGasPressureSink2(                                            p_const=100000, xi_const={1/9,1/9,1/9,1/9,1/9,1/9,1/9,1/9,1/9})
                                                                                        annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={290,-74})));
  Modelica.Blocks.Sources.Ramp massFlowRate4(
    offset=1,
    startTime=5,
    duration=10,
    height=-2)   annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={150,-54})));
  Modelica.Blocks.Sources.Ramp Temperature3(
    duration=1,
    startTime=1,
    height=50,
    offset=273.15 + 150)
                        annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={150,-86})));
  E_Filter.E_Filter_L2_empirical
                        e_Filter_dynamic1(
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.GenericGeometry,
    use_dynamicMassbalance=true,
    A_el=200,
    redeclare model PressureLoss =
        Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L2
        (                                                                                                    Delta_p_nom=100))
              annotation (Placement(transformation(extent={{216,-84},{236,-64}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperatureTop4(T=293.15)
                annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={150,-26})));
  BoundaryConditions.BoundaryGas_Txim_flow                  idealGasFlowSource_XRG5(
    variable_m_flow=false,
    variable_T=false,
    T_const=473.15,
    xi_const={0.01,0,0.73,0,0.065,0.036,0,0.13,0.0},
    m_flow_const=10)
                    annotation (Placement(transformation(extent={{80,-192},{100,-172}})));
  BoundaryConditions.BoundaryGas_pTxi                  idealGasPressureSink3(                                                            xi_const={0.01,0,0.73,0,0.065,0.036,0,0.13,0.0}, p_const=101300)
                                                                                        annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={250,-182})));
  E_Filter.E_Filter_L2_detailed
                       e_Filter_dynamic2(
                                        redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.GenericGeometry,                                                            use_dynamicMassbalance=true,
    T_start=473.15,
    xi_start={0.01,0,0.73,0,0.065,0.036,0,0.13,0.0},
    p_start=1.023e5,
    redeclare model PressureLoss =
        Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L2
        (                                                                                                    Delta_p_nom=100),
    m_flow_nom=10)                                                                                                     annotation (Placement(transformation(extent={{176,-192},{196,-172}})));
  Modelica.Blocks.Sources.Ramp U_applied(
    duration=10,
    height=20e3,
    startTime=10,
    offset=1000) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={150,-152})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperatureTop5(T=293.15)
                annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={150,-120})));
  VolumesValvesFittings.Pipes.PipeFlowGas_L4_Simple pipeFlowGas_L4_Simple(
    redeclare model PressureLoss =
        Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    initType=ClaRa.Basics.Choices.Init.noInit,
    xi_start={0.01,0,0.73,0,0.065,0.036,0,0.13,0.0},
    xi_nom={0.01,0,0.73,0,0.065,0.036,0,0.13,0.0},
    p_nom=1e5*ones(pipeFlowGas_L4_Simple.N_cv),
    T_nom=473.15*ones(pipeFlowGas_L4_Simple.N_cv),
    T_start=473.15*ones(pipeFlowGas_L4_Simple.N_cv),
    Delta_p_nom=100,
    frictionAtOutlet=true,
    redeclare model HeatTransfer =
        Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Constant_L4
        (                                                                                                   alpha_nom=10),
    length=1,
    diameter_i=0.5,
    p_start=1.033e5*ones(pipeFlowGas_L4_Simple.N_cv),
    m_flow_nom=10) annotation (Placement(transformation(extent={{124,-186},{152,-176}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperatureTop6(T=293.15)
                annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={90,-138})));
  Adapters.Scalar2VectorHeatPort            scalar2VectorHeatPort2(
    N=3,
    equalityMode="Equal Temperatures",
    length=10,
    Delta_x=ClaRa.Basics.Functions.GenerateGrid(
        {0,0},
        10,
        3))         annotation (Placement(transformation(extent={{106,-148},{126,-128}})));

  Regression regression(
    deSO_ideal_L1_1_m_flow_out = deSO_ideal_L1_1.summary.outlet.m_flow,
    deSO_ideal_L1_1_T_out = deSO_ideal_L1_1.summary.outlet.T,
    deNOx_m_flow_out = deNOx.summary.outlet.m_flow,
    deNOx_T_out = deNOx.summary.outlet.T,
    deNOx_NH3_port_m_flow_out = deNOx_NH3_port.summary.outlet.m_flow,
    deNOx_NH3_port_T_out = deNOx_NH3_port.summary.outlet.T,
    e_Filter_simple_m_flow_out = e_Filter_dynamic.summary.outlet.m_flow,
    e_Filter_simple_T_out = e_Filter_dynamic.summary.outlet.T,
    e_Filter_empirical_m_flow_out = e_Filter_dynamic1.summary.outlet.m_flow,
    e_Filter_empirical_T_out = e_Filter_dynamic1.summary.outlet.T,
    e_Filter_detailed_m_flow_out = e_Filter_dynamic2.summary.outlet.m_flow,
    e_Filter_detailed_T_out = e_Filter_dynamic2.summary.outlet.T,
    pipeFlowGas_L4_Simple_T_out = pipeFlowGas_L4_Simple.summary.outlet.T,
    pipeFlowGas_L4_Simple_p_out = pipeFlowGas_L4_Simple.summary.outlet.p) annotation (Placement(transformation(extent={{162,212},{182,232}})));

equation
  connect(gasFlowSource_T.gas_a, deSO_ideal_L1_1.inlet) annotation (Line(
      points={{-58,36},{-38,36}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(deSO_ideal_L1_1.outlet, gasSink_pT.gas_a) annotation (Line(
      points={{-18,36},{-2,36}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(massFlowRate2.y,idealGasFlowSource_XRG2. m_flow) annotation (Line(
      points={{-57,-30},{-48,-30},{-48,-42},{-30,-42}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(Temperature2.y,idealGasFlowSource_XRG2. T)
                                             annotation (Line(
      points={{-57,-62},{-48,-62},{-48,-48},{-30,-48}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(idealGasFlowSource_XRG2.gas_a,deNOx. inlet) annotation (Line(
      points={{-10,-48},{4,-48}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(deNOx.outlet,idealGasPressureSink_XRG1. gas_a) annotation (Line(
      points={{24,-48},{42,-48}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(fixedTemperatureTop1.port,deNOx. heat) annotation (Line(
      points={{-58,0},{8.8,0},{8.8,-38.4}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(massFlowRate1.y,idealGasFlowSource_XRG1. m_flow) annotation (Line(
      points={{-59,-156},{-50,-156},{-50,-168},{-32,-168}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(Temperature1.y,idealGasFlowSource_XRG1. T)
                                             annotation (Line(
      points={{-59,-188},{-50,-188},{-50,-174},{-32,-174}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowRate3.y,idealGasFlowSource_XRG3. m_flow) annotation (Line(
      points={{-59,-124},{-46,-124},{-46,-136},{-32,-136}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(idealGasFlowSource_XRG3.gas_a, deNOx_NH3_port.NH3_inlet) annotation (Line(
      points={{-12,-142},{14,-142},{14,-164}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(idealGasFlowSource_XRG1.gas_a, deNOx_NH3_port.inlet) annotation (Line(
      points={{-12,-174},{4,-174}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(deNOx_NH3_port.outlet, idealGasPressureSink.gas_a) annotation (Line(
      points={{24,-174},{42,-174}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(fixedTemperatureTop2.port, deNOx_NH3_port.heat) annotation (Line(
      points={{-60,-94},{8.8,-94},{8.8,-164.4}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(massFlowRate.y,idealGasFlowSource_XRG. m_flow) annotation (Line(
      points={{161,42},{170,42},{170,28},{180,28}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(Temperature.y,idealGasFlowSource_XRG. T)
                                             annotation (Line(
      points={{161,10},{170,10},{170,22},{180,22}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(idealGasFlowSource_XRG.gas_a,e_Filter_dynamic. inlet) annotation (
      Line(
      points={{200,22},{216,22}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(e_Filter_dynamic.outlet, idealGasPressureSink1.gas_a) annotation (Line(
      points={{236,22},{280,22}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(fixedTemperatureTop3.port,e_Filter_dynamic. heat) annotation (Line(
      points={{160,70},{226,70},{226,32}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(massFlowRate4.y, idealGasFlowSource_XRG4.m_flow) annotation (Line(
      points={{161,-54},{170,-54},{170,-68},{180,-68}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(Temperature3.y, idealGasFlowSource_XRG4.T) annotation (Line(
      points={{161,-86},{170,-86},{170,-74},{180,-74}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(idealGasFlowSource_XRG4.gas_a, e_Filter_dynamic1.inlet) annotation (Line(
      points={{200,-74},{216,-74}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(e_Filter_dynamic1.outlet, idealGasPressureSink2.gas_a) annotation (Line(
      points={{236,-74},{280,-74}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(fixedTemperatureTop4.port, e_Filter_dynamic1.heat) annotation (Line(
      points={{160,-26},{226,-26},{226,-64}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(U_applied.y, e_Filter_dynamic2.U_applied) annotation (Line(
      points={{161,-152},{178.6,-152},{178.6,-170.8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(e_Filter_dynamic2.outlet, idealGasPressureSink3.gas_a) annotation (Line(
      points={{196,-182},{240,-182}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(fixedTemperatureTop5.port, e_Filter_dynamic2.heat) annotation (Line(
      points={{160,-120},{186,-120},{186,-172}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(idealGasFlowSource_XRG5.gas_a, pipeFlowGas_L4_Simple.inlet) annotation (Line(
      points={{100,-182},{124,-182},{124,-181}},
      color={118,106,98},
      thickness=0.5));
  connect(pipeFlowGas_L4_Simple.outlet, e_Filter_dynamic2.inlet) annotation (Line(
      points={{152,-181},{166,-181},{166,-182},{176,-182}},
      color={118,106,98},
      thickness=0.5));
  connect(fixedTemperatureTop6.port, scalar2VectorHeatPort2.heatScalar) annotation (Line(points={{100,-138},{103,-138},{106,-138}}, color={191,0,0}));
  connect(scalar2VectorHeatPort2.heatVector, pipeFlowGas_L4_Simple.heat) annotation (Line(
      points={{126,-138},{138,-138},{138,-177}},
      color={167,25,48},
      thickness=0.5));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-220},{320,100}},
        initialScale=0.1),     graphics={
                                Text(
          extent={{-100,86},{-26,76}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=10,
          textString="________________________________________________________________
PURPOSE:
>>Regression tester for the flue gas cleaning components"),
                                Text(
          extent={{-100,100},{30,90}},
          lineColor={0,128,0},
          fontSize=34,
          textString="TESTED -- 2014-10-08 //LN"),
        Rectangle(
          extent={{-100,100},{320,-220}},
          lineColor={115,150,0},
          lineThickness=0.5)}),           Commands(file="../../plot_DeSO.mos"
        "plot_DeSO"),
    Icon(coordinateSystem(
        preserveAspectRatio=false,
        initialScale=0.1)),
    experiment(StopTime=10000));
end Test_FlueGasCleaning;
