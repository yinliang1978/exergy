within Exergy.XClaRa.Components.Furnace.BaseClasses;
partial model CombustionChamberBase
  import ClaRa;

    //## S U M M A R Y   D E F I N I T I O N ###################################################################
  model Outline
    //  parameter Boolean showExpertSummary annotation(Dialog(hide));
    extends ClaRa.Basics.Icons.RecordIcon;
    input ClaRa.Basics.Units.Volume volume "Volume";
    input ClaRa.Basics.Units.Area A_cross "Free cross sectional area";
    input ClaRa.Basics.Units.Area A_wall "Wall area";
    input ClaRa.Basics.Units.Length height "Height of volume";
    input ClaRa.Basics.Units.Mass m "Mass inside volume";
    input ClaRa.Basics.Units.MassFlowRate m_flow_fuel_burned
      "Burned fuel mass flow rate";
    input ClaRa.Basics.Units.MassFlowRate m_flow_oxygen_burned
      "Burned oxygen mass flow rate";
    input ClaRa.Basics.Units.MassFlowRate m_flow_oxygen_req
      "Required O2 flow rate for stochiometric combustion";
    input ClaRa.Basics.Units.MassFlowRate m_flow_air_req
      "Required air flow rate for stochiometric combustion";
    input Real lambdaComb "Excess air";
    input Real NOx_fraction "NOx fraction at outlet";
    input Real CO_fraction "CO fraction at outlet";
    input ClaRa.Basics.Units.EnthalpyMassSpecific LHV "Lower heating value";
    input ClaRa.Basics.Units.HeatFlowRate Q_combustion "Combustion Heat";
    input ClaRa.Basics.Units.Velocity w_migration "Particle migration speed";
    input ClaRa.Basics.Units.Time t_dwell_flueGas "Flue gas dwelltime";
    input ClaRa.Basics.Units.Time burning_time "Burning time";
    input Real unburntFraction "Fuel diffusity";
    input ClaRa.Basics.Units.Temperature T_out "Outlet temperature";
    input ClaRa.Basics.Units.EnthalpyMassSpecific h_out
      "Flue gas enthalpy at outlet";
  end Outline;

  model Fuel
    extends ClaRa.Basics.Icons.RecordIcon;
    input ClaRa.Basics.Units.MassFlowRate m_flow "Mass flow rate"
      annotation (Dialog);
    input ClaRa.Basics.Units.Temperature T "Temperature" annotation (Dialog);
    input ClaRa.Basics.Units.Pressure p "Pressure" annotation (Dialog);
    input ClaRa.Basics.Units.HeatCapacityMassSpecific cp
      "Specific heat capacity" annotation (Dialog);
  end Fuel;

  model Slag
    extends ClaRa.Basics.Icons.RecordIcon;
    input ClaRa.Basics.Units.MassFlowRate m_flow "Mass flow rate"
      annotation (Dialog);
    input ClaRa.Basics.Units.Temperature T "Temperature" annotation (Dialog);
    input ClaRa.Basics.Units.Pressure p "Pressure" annotation (Dialog);
  end Slag;

  model Flow
    extends ClaRa.Basics.Icons.RecordIcon;
    ClaRa.Basics.Records.FlangeGas flueGas;
    Fuel fuel;
    Slag slag;
  end Flow;

  model Summary
    extends ClaRa.Basics.Icons.RecordIcon;
    Outline outline;
    Flow inlet;
    Flow outlet;
  end Summary;

//## P A R A M E T E R S #######################################################################################
   //__________________________/ Media definintions \______________________________________________
  outer ClaRa.SimCenter simCenter;
  inner parameter ClaRa.Basics.Media.Fuel.PartialFuel fuelType=simCenter.fuelModel1
    "Fuel elemental composition used for combustion"                                                                                 annotation(choices(choice=simCenter.fuelModel
        "Fuel model 1 as defined in simCenter"),                                            Dialog(group="Media Definitions"));
  parameter ClaRa.Basics.Media.Fuel.PartialSlag slagType=simCenter.slagModel
    "Slag properties"                                                                          annotation(choices(choice=simCenter.slagModel
        "Slag model 1 as defined in simCenter"),                                            Dialog(group="Media Definitions"));
  inner parameter TILMedia.GasTypes.BaseGas flueGas = simCenter.flueGasModel
    "Flue gas model used in component"                                                                          annotation(choicesAllMatching, Dialog(group="Media Definitions"));

  parameter Integer slagTemperature_calculationType=1
    "Calculation type of outflowing slag temperature"                                                   annotation (Dialog(group="Slag temperature definitions"), choices(
      choice=1 "Fixed slag temperature",
      choice=2 "Outlet flue gas temperature",
      choice=3 "Mean flue gas temperature",
      choice=4 "Inlet flue gas temperature"));

  inner parameter ClaRa.Basics.Units.Temperature T_slag=900
    "Constant slag outlet temperature"                                                         annotation (Dialog(enable=(slagTemperature_calculationType ==
          1),group="Slag temperature definitions"));

  //__________________________/ HeatTransfer \______________________________________________
    replaceable model HeatTransfer_Wall =
      ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Gas_HT.Radiation.Radiation_gas2Wall_L2
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.HeatTransferBaseGas
    "1st: choose geometry definition | 2nd: edit corresponding record" annotation (Dialog(group="Heat Transfer"), choicesAllMatching=true);

  replaceable model HeatTransfer_Top =
      ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Gas_HT.Radiation.Radiation_gas2Wall_L2
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.HeatTransferBaseGas
    "1st: choose geometry definition | 2nd: edit corresponding record" annotation (Dialog(group="Heat Transfer"), choicesAllMatching=true);
  inner parameter Modelica.SIunits.Time Tau_rad= 0.1 "Radiation time constant" annotation(Dialog(group="Heat Transfer"));

//________________________/Chemistry \________________________________________________________
  replaceable model ReactionZone =
      ClaRa.Components.Furnace.ChemicalReactions.CoalReactionZone
    constrainedby
    ClaRa.Components.Furnace.ChemicalReactions.PartialReactionZone
    "Model to regard chemical reactions" annotation (Dialog(group=
          "Combustion"), choicesAllMatching=true);

  replaceable model Burning_time =
      ClaRa.Components.Furnace.GeneralTransportPhenomena.BurningTime.ConstantBurningTime
    constrainedby
    ClaRa.Components.Furnace.GeneralTransportPhenomena.BurningTime.PartialBurningTime
    "Model for the buring time" annotation (Dialog(group="Combustion"),
      choicesAllMatching=true);

  replaceable model ParticleMigration =
      ClaRa.Components.Furnace.GeneralTransportPhenomena.ParticleMigration.MeanMigrationSpeed
    constrainedby
    ClaRa.Components.Furnace.GeneralTransportPhenomena.ParticleMigration.PartialMigrationSpeed
    "Model for the particle migration speed" annotation (Dialog(group=
          "Combustion"), choicesAllMatching=true);

  //__________________________/ Geometry \______________________________________________
  replaceable model Geometry =
  ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.HollowBlock constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.GenericGeometry(flowOrientation = ClaRa.Basics.Choices.GeometryOrientation.vertical, height_fill=-1)
    "1st: choose geometry definition | 2nd: edit corresponding record"
    annotation (Dialog(group="Geometry"), choicesAllMatching=true);

  //__________________/ Parameter \_______________________________________________

  inner parameter Modelica.SIunits.MassFlowRate m_flow_nom= 10
    "Nominal mass flow rates at inlet"                                                            annotation(Dialog(group="Nominal Values"));

  //_______________________/ Start values \_____________________________________________________________
   parameter ClaRa.Basics.Units.Pressure p_start_flueGas_out=1e5
    "Start pressure at outlet"                                                              annotation(Dialog(tab="Initialisation"));
  parameter ClaRa.Basics.Units.Temperature T_start_flueGas_out=700
    "Start temperature at outlet"                                                                annotation(Dialog(tab="Initialisation"));
  inner parameter Modelica.SIunits.Temperature T_top_initial= T_start_flueGas_out
    "Initial temperature of top volume"                                                                               annotation(Dialog(tab="Initialisation"));
  parameter ClaRa.Basics.Units.MassFraction xi_start_flueGas_out[flueGas.nc - 1]={0.01,0,0.1,0,0.74,0.13,0,0.02,0}
    "Start composition of flue gas"                                                         annotation(Dialog(tab="Initialisation"));
  //   parameter ClaRa.Basics.Units.VolumeFlowRate V_flow_flueGas_in_start=1 annotation(Dialog(tab="Initialisation"));
  parameter ClaRa.Basics.Units.VolumeFlowRate V_flow_flueGas_out_start=-15
    "Start volume flow at outlet"                                                                        annotation(Dialog(tab="Initialisation"));
  final parameter Modelica.SIunits.SpecificEnthalpy h_start = TILMedia.GasFunctions.specificEnthalpy_pTxi(flueGas, p_start_flueGas_out, T_start_flueGas_out, xi_start_flueGas_out)
    "Start flue gas enthalpy"                                                               annotation(Dialog(tab="Initialisation"));

  constant Real T_0=298.15 "Reference temperature";

//## V A R I A B L E   P A R T##################################################################################

protected
  inner ClaRa.Basics.Units.MassFraction xi_fuel
    "amount of fuel per flue gas mass";
  inner constant ClaRa.Basics.Units.Time Tau = 0.001
    "time constant for heat transfer temperature delay";

//________________/ FlueGas Composition \_____________________
public
 inner ClaRa.Basics.Units.MassFraction xi_flueGas[flueGas.nc - 1]
    "Flue gas composition ";
//________________/ Fuel Composition \_____________________
  ClaRa.Basics.Units.MassFraction xi_fuel_out[fuelType.nc - 1]
    "Fuel outlet composition";

//_____________/ Calculated quantities \_________________________________
  inner ClaRa.Basics.Units.Area A_cross=geo.A_front "Cross section";
  inner ClaRa.Basics.Units.VolumeFlowRate V_flow_flueGas_in "Inlet volume flow";
  inner ClaRa.Basics.Units.VolumeFlowRate V_flow_flueGas_out
    "Outlet volume flow";

  Modelica.SIunits.Mass mass "Gas mass";
  ClaRa.Basics.Units.HeatFlowRate Q_flow_top "Heat flow from top section";
  ClaRa.Basics.Units.HeatFlowRate Q_flow_bottom "Heat flow from bottom section";
  ClaRa.Basics.Units.HeatFlowRate Q_flow_wall "Heat flow from walls";

  inner ClaRa.Basics.Units.MassFraction xi_fuel_in[fuelType.nc - 1]
    "Fuel inlet composition";
  ClaRa.Basics.Units.EnthalpyMassSpecific h_flueGas_out
    "Gas outlet specific enthalpy";
  ClaRa.Basics.Units.EnthalpyMassSpecific h_flueGas_out_del
    "Gas outlet specific enthalpy - delayed";

  ClaRa.Basics.Units.MassFraction xi_flueGas_del[flueGas.nc - 1]
    "Flue gas outlet composition - dalayed";

//___________________/ Molar flow rates of educts and products \_____________
//_________/Educts\__________________
protected
  Modelica.SIunits.MolarFlowRate n_flow_C "Carbon molar flow";
  Modelica.SIunits.MolarFlowRate n_flow_H "Hydrogen molar flow";
  Modelica.SIunits.MolarFlowRate n_flow_O "Oxygen molar flow";
  Modelica.SIunits.MolarFlowRate n_flow_N "Nitrogen molar flow";
  Modelica.SIunits.MolarFlowRate n_flow_S "Sulphur molar flow";
  Modelica.SIunits.MolarFlowRate n_flow_Ash "Ash molar flow";
  Modelica.SIunits.MolarFlowRate n_flow_H2O "Water molar flow";

  ClaRa.Basics.Units.MassFraction xi_flueGas_id[flueGas.nc - 1]
    "Flue gas composition at ideal combustion";
  Real sum_xi "Sum of flue gas outlet's species mass fractions";
  Real sum_comp "Sum of flue gas species mass fractions after ideal combustion";
  Real m_flow_fuel_id "Fuel consumed at ideal combustion";
  Real m_flow_flueGas_id "Flue gas produced at ideal combustion";
  Real Delta_h_f "Formation enthalpy of used fuel";

public
  ClaRa.Basics.Units.EnthalpyMassSpecific LHV(start=(33907*fuelType.defaultComposition[1] + 142324*(fuelType.defaultComposition[2] - fuelType.defaultComposition[3]/8.) + 10465*fuelType.defaultComposition[5] - 2512*((1 - sum(fuelType.defaultComposition)) + 9*fuelType.defaultComposition[2]))*1000);
  Modelica.SIunits.SpecificHeatCapacity cp;
  ClaRa.Basics.Units.Time t_dwell_flueGas=geo.height/particleMigration.w
    "Flow time in z-direction";
  Real unburntFraction
    "Quantatity describes how much unburned fuel leaves control volume together with the fluegas to neighbor cell";
  ClaRa.Basics.Units.MassFlowRate m_flow_oxygen_req(min=1e-15)
    "Required O2 flow rate for stochiometric combustion";
  ClaRa.Basics.Units.MassFlowRate m_flow_air_req
    "Required combustion air flow rate for stochiometric combustion determined in dependence of m_flow_oxygen_req and actual Xi of gasInlet";

  ClaRa.Basics.Units.MassFlowRate m_flow_fuel_burned "Burned flue mass flow";
  ClaRa.Basics.Units.MassFlowRate m_flow_oxygen_burned
    "Oxygen used for combustion";

  Real lambdaComb "Excess air";

    //________________________/ Connectors \_______________________________________________________
  ClaRa.Basics.Interfaces.FuelSlagFlueGas_inlet inlet(
    flueGas(final Medium=flueGas),
    final fuelType=fuelType,
    final slagType=slagType)
    annotation (Placement(transformation(extent={{-130,-110},{-110,-90}}),
        iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-160,-100})));
  ClaRa.Basics.Interfaces.FuelSlagFlueGas_outlet outlet(
    flueGas(final Medium=flueGas, m_flow(start=-1)),
    final fuelType=fuelType,
    final slagType=slagType)                                                                        annotation (Placement(transformation(extent={{-130,90},
            {-110,110}}),
        iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-160,100})));
  ClaRa.Basics.Interfaces.HeatPort_a
                                   heat_wall
    annotation (Placement(transformation(extent={{290,-10},{310,10}})));
  ClaRa.Basics.Interfaces.HeatPort_a
                                   heat_top
    annotation (Placement(transformation(extent={{10,90},{30,110}})));
  ClaRa.Basics.Interfaces.HeatPort_a
                                   heat_bottom
     annotation (Placement(transformation(extent={{10,-90},{30,-110}})));

   //_____________________/ Media Objects \_________________________________
protected
    TILMedia.Gas_pT     flueGasInlet(p=inlet.flueGas.p, T= actualStream(inlet.flueGas.T_outflow), xi=actualStream(inlet.flueGas.xi_outflow),
      gasType=flueGas)
      annotation (Placement(transformation(extent={{-130,-88},{-110,-68}})));

public
  TILMedia.Gas_pT ideal_combustion(
    p=outlet.flueGas.p,
    gasType=flueGas,
    T=298.15)
    annotation (Placement(transformation(extent={{-100,-12},{-80,8}})));

public
  inner Geometry geo annotation(Placement(transformation(extent={{-250,50},{-230,
            70}})));

  inner HeatTransfer_Wall heattransfer_wall(heatSurfaceAlloc=1)  annotation(Placement(transformation(extent={{232,-10},
            {252,10}})));

  inner HeatTransfer_Top heattransfer_top annotation(Placement(transformation(extent={{94,50},
            {74,70}})));
  inner
    ClaRa.Components.Furnace.GeneralTransportPhenomena.ThermalCapacities.ThermalLowPass
                                                                                            radiationTimeConstant(T_out_initial=T_start_flueGas_out, Tau=Tau_rad)   annotation (Placement(transformation(extent={{32,50},
             {52,70}})));

  ReactionZone reactionZone(xi_fuel_in=xi_fuel_in, xi_flueGas=xi_flueGas)
    annotation (Placement(transformation(extent={{-130,-10},{-110,10}})));
  Burning_time burning_time
    annotation (Placement(transformation(extent={{-54,28},{-34,48}})));
  ParticleMigration particleMigration
    annotation (Placement(transformation(extent={{-54,-10},{-34,10}})));

  ClaRa.Basics.Interfaces.EyeOutGas
                           eyeOut annotation (Placement(transformation(extent={{-280,78},
            {-308,102}}),         iconTransformation(extent={{-290,70},{-310,90}})));
protected
           ClaRa.Basics.Interfaces.EyeInGas
                                   eye_int
                                annotation (Placement(transformation(extent={{-254,84},
            {-266,96}}),      iconTransformation(extent={{240,-64},{232,-56}})));

public
  parameter Boolean showData = false
    "True, if characteristic data shall be visualised in model icon"                                   annotation(Dialog(tab="Summary and Visualisation"));

initial equation

  h_flueGas_out_del = h_start;
  xi_flueGas_del = xi_start_flueGas_out;

equation

  der(h_flueGas_out_del) = 1/Tau*(h_flueGas_out-h_flueGas_out_del);
  der(xi_flueGas_del) = 1/Tau*(xi_flueGas - xi_flueGas_del);

  //____________/ Xi_outflow of Fuel and FlueGas \__________________
  //dummy values for inlets since flow reversal is not allowed
  outlet.fuel.xi_outflow = xi_fuel_out;
  inlet.fuel.xi_outflow = xi_fuel_out;

  //_____________/ Pressure \______________________________________________
  inlet.flueGas.p = outlet.flueGas.p;
  inlet.fuel.p = outlet.fuel.p;
  inlet.slag.p = outlet.slag.p;

  //____________/ Heat port temperatures and Q_flows \____________________________
   Q_flow_wall = heat_wall.Q_flow;
   Q_flow_top = heat_top.Q_flow;
   Q_flow_bottom = heat_bottom.Q_flow;

 //_____________/ Connections \______________________________________________
  connect(heat_top, radiationTimeConstant.heat_in) annotation (Line(
      points={{20,100},{20,60},{32,60}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(radiationTimeConstant.heat_out, heattransfer_top.heat) annotation (
      Line(
      points={{52,60},{74,60}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(heattransfer_wall.heat, heat_wall) annotation (Line(
      points={{252,0},{300,0}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(eye_int,eyeOut)
                         annotation (Line(
      points={{-260,90},{-294,90}},
      color={190,190,190},
      smooth=Smooth.None));
 annotation (Icon(coordinateSystem(preserveAspectRatio=false,extent={{-300,-100},
            {300,100}}),
                   graphics),               Diagram(coordinateSystem(
          preserveAspectRatio=false,extent={{-300,-100},{300,100}}),
                                                    graphics),
    Documentation(info="<html>
<p><b>Model description: </b>Base class for burner and furnace sections</p>
<p><b>Contact:</b> Lasse Nielsen, TLK-Thermo GmbH</p>
</html>"));
end CombustionChamberBase;
