within Exergy.XClaRa.Components.FlueGasCleaning.Desulfurization.Check;
model Test_Desulfurization_ideal
  extends ClaRa.Basics.Icons.PackageIcons.ExecutableExampleb50;
  Desulfurization_L1_ideal deSO_ideal_L1_1(
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.HollowCylinder,
    m_flow_nom=530,
    xi_start={0,0,0.21,0.00099,0.7,0.0393,0,0.0367,0},
    SOx_separationRate=0.95,
    initType=ClaRa.Basics.Choices.Init.noInit,
    T_start=395.843,
    p_start=101800,
    useStabilisedMassFlow=true,
    redeclare model PressureLoss =
        ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L2)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

  BoundaryConditions.BoundaryGas_Txim_flow gasFlowSource_T(
    m_flow_const=551.153,
    xi_const={0,0,0.21,0.00099,0.7,0.0393,0,0.0367,0},
    T_const=395.843) annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
  BoundaryConditions.BoundaryGas_pTxi gasSink_pT(
    xi_const={0,0,0.21,0.00099,0.7,0.0393,0,0.0367,0},
    p_const=101800,
    T_const=293.15) annotation (Placement(transformation(extent={{46,-10},{26,10}})));
  inner ClaRa.SimCenter simCenter(
    contributeToCycleSummary=true,
    redeclare TILMedia.GasTypes.FlueGasTILMedia flueGasModel,
    showExpertSummary=true) annotation (Placement(transformation(
          extent={{68,68},{88,88}})));
equation

  connect(gasFlowSource_T.gas_a, deSO_ideal_L1_1.inlet) annotation (Line(
      points={{-30,0},{-10,0}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(deSO_ideal_L1_1.outlet, gasSink_pT.gas_a) annotation (Line(
      points={{10,0},{26,0}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics={
                                Text(
          extent={{-100,86},{-26,76}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=12,
          textString="________________________________________________________________
PURPOSE:
>>Tester for the Desulfurization component"),
                                Text(
          extent={{-100,100},{30,90}},
          lineColor={0,128,0},
          fontSize=34,
          textString="TESTED -- 2014-10-08 //LN")}),
                                          Commands(file="../../plot_DeSO.mos"
        "plot_DeSO"));
end Test_Desulfurization_ideal;
