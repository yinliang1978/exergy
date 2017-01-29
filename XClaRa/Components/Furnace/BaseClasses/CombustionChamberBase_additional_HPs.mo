within Exergy.XClaRa.Components.Furnace.BaseClasses;
partial model CombustionChamberBase_additional_HPs
  import ClaRa;

 //________________________/ Connectors \_______________________________________________________
  ClaRa.Basics.Interfaces.HeatPort_a
                                   heat_CarrierTubes
    annotation (Placement(transformation(extent={{190,90},{210,110}})));
  ClaRa.Basics.Interfaces.HeatPort_a
                                   heat_TubeBundle
    annotation (Placement(transformation(extent={{190,-90},{210,-110}})));

 //________________________/ replacable modells for heat transfer, pressure loss and geometry \________________________
  replaceable model HeatTransfer_CarrierTubes =
      ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Gas_HT.Convection.Convection_carrierTubes_L2
  constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.HeatTransferBaseGas
    "1st: choose geometry definition | 2nd: edit corresponding record"
    annotation (Dialog(group="Heat Transfer"), choicesAllMatching=
        true);

  replaceable model HeatTransfer_TubeBundle =
      ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Gas_HT.Convection.Convection_tubeBank_L2
  constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.HeatTransferBaseGas
    "1st: choose geometry definition | 2nd: edit corresponding record"
    annotation (Dialog(group="Heat Transfer"), choicesAllMatching=
        true);

    inner HeatTransfer_CarrierTubes  heattransfer_CarrierTubes(heatSurfaceAlloc=3) annotation(Placement(transformation(extent={{10,-10},
            {-10,10}},
        rotation=180,
        origin={162,60})));
    inner HeatTransfer_TubeBundle heattransfer_TubeBundle(heatSurfaceAlloc=2) annotation(Placement(transformation(extent={{10,10},
            {-10,-10}},
        rotation=180,
        origin={162,-60})));

  ClaRa.Basics.Units.HeatFlowRate Q_flow_CarrierTubes
    "Heat flow from carrier tubes";
  ClaRa.Basics.Units.HeatFlowRate Q_flow_TubeBundle
    "Heat flow from tube bundle";

equation
  //____________/ Heat port temperatures and Q_flows \____________________________
   Q_flow_CarrierTubes = heat_CarrierTubes.Q_flow;
   Q_flow_TubeBundle =  heat_TubeBundle.Q_flow;

  //_____________/ Connections \______________________________________________
  connect(heattransfer_CarrierTubes.heat, heat_CarrierTubes) annotation (Line(
      points={{172,60},{200,60},{200,100}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(heattransfer_TubeBundle.heat, heat_TubeBundle) annotation (Line(
      points={{172,-60},{200,-60},{200,-100}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-300,-100},
            {300,100}}), graphics), Documentation(info="<html>
<p><b>Model description: </b>Base class for furnace sections with additional heat ports</p>
<p><b>Contact:</b> Lasse Nielsen, TLK-Thermo GmbH</p>
</html>"));
end CombustionChamberBase_additional_HPs;
