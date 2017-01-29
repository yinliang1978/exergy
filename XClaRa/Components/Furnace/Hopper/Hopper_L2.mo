within Exergy.XClaRa.Components.Furnace.Hopper;
model Hopper_L2 "Model for a hopper section of a combustion chamber"
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

  import ClaRa;
  extends ClaRa.Basics.Icons.Hopper;
  //inner parameter ClaRa.Basics.Units.Temperature SlagTemperature=900;
extends ClaRa.Components.Furnace.BaseClasses.HopperBase(geo(
        flowOrientation=ClaRa.Basics.Choices.GeometryOrientation.vertical));

inner parameter Boolean useHomotopy=simCenter.useHomotopy
    "True, if homotopy method is used during initialisation"  annotation(Dialog(tab="Initialisation"));

  Real sum_xi "Sum of inlet components";
  Real drhodt "Density derivative";

  ClaRa.Basics.Units.MassFraction xi_flueGas_in_del[flueGas.nc - 1]
    "Flue gas mixture composition";

  ClaRa.Basics.Units.MassFlowRate m_flow_in_del
    "Pseudo state for inlet mass flow";
  ClaRa.Basics.Units.MassFlowRate m_flow_out_del
    "Pseudo state for outlet mass flow";
  ClaRa.Basics.Units.Temperature T_bulk_del "Pseudo state for bulk temperature";
  ClaRa.Basics.Units.DensityMassSpecific rho_bulk_del
    "Pseudo state for bulk density";

  model Outline
  //  parameter Boolean showExpertSummary annotation(Dialog(hide));
    extends ClaRa.Basics.Icons.RecordIcon;
    input ClaRa.Basics.Units.Volume
                                volume "Volume";
    input ClaRa.Basics.Units.Area
                              A_cross "Cross sectional area";
    input ClaRa.Basics.Units.Area
                              A_wall "Wall area";
    input ClaRa.Basics.Units.Length
                                height "Height of volume";
    input ClaRa.Basics.Units.Mass
                              mass "Mass inside volume";
    input ClaRa.Basics.Units.Temperature
                                     T_out "Outlet temperature";
    input ClaRa.Basics.Units.EnthalpyMassSpecific
                                              h_out
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

   Summary summary(
   outline(
    volume = geo.volume,
    A_cross = geo.width*geo.height,
    A_wall = geo.width*geo.length*2+geo.length*geo.height*2,
    height=geo.height,
    mass = mass,
    T_out = flueGasOutlet.T,
    h_out = flueGasOutlet.h),
    inlet(
      flueGas(
        m_flow=inlet.flueGas.m_flow,
        T=actualStream(inlet.flueGas.T_outflow),
        p=inlet.flueGas.p,
        h=flueGasInlet.h,
        xi=actualStream(inlet.flueGas.xi_outflow),
        H_flow=flueGasInlet.h*inlet.flueGas.m_flow),
      fuel(
        m_flow=inlet.fuel.m_flow,
        T=actualStream(inlet.fuel.T_outflow),
        p=inlet.fuel.p,
        cp=inlet.fuelType.cp),
      slag(
        m_flow=inlet.slag.m_flow,
        T=actualStream(inlet.slag.T_outflow),
        p=inlet.slag.p)),
    outlet(
      flueGas(
        m_flow=-outlet.flueGas.m_flow,
        T=actualStream(outlet.flueGas.T_outflow),
        p=outlet.flueGas.p,
        h=h_flueGas_out,
        xi=actualStream(outlet.flueGas.xi_outflow),
        H_flow=-h_flueGas_out*outlet.flueGas.m_flow),
      fuel(
        m_flow=-outlet.fuel.m_flow,
        T=actualStream(outlet.fuel.T_outflow),
        p=outlet.fuel.p,
        cp=outlet.fuelType.cp),
      slag(
        m_flow=outlet.slag.m_flow,
        T=actualStream(outlet.slag.T_outflow),
        p=outlet.slag.p))) annotation (Placement(transformation(extent={{274,-102},{300,-76}})));

//___________________/ Media Objects \_________
   TILMedia.Gas_ph bulk(
     p(start=p_start_flueGas_out) = outlet.flueGas.p,
     xi=xi_flueGas_del,
     gasType=flueGas,
     h=h_flueGas_out_del)
       annotation (Placement(transformation(extent={{-130,26},{-110,46}})));

//___________________/ iCom record \__________________
protected
  inner ClaRa.Basics.Records.IComGas_L2 iCom(
    m_flow_nom=m_flow_nom,
    T_bulk=flueGasOutlet.T,
    p_bulk=bulk.p,
    fluidPointer_in=flueGasInlet.gasPointer,
    fluidPointer_bulk=flueGasOutlet.gasPointer,
    fluidPointer_out=flueGasOutlet.gasPointer,
    mediumModel=flueGas,
    p_in=inlet.flueGas.p,
    T_in=T_bulk_del,
    m_flow_in=m_flow_in_del,
    V_flow_in=V_flow_flueGas_in,
    xi_in=xi_flueGas_del,
    p_out=outlet.flueGas.p,
    T_out=T_bulk_del,
    m_flow_out=m_flow_out_del,
    V_flow_out=V_flow_flueGas_out,
    xi_out=xi_flueGas_del) annotation (Placement(transformation(extent={{244,-102},{268,-76}})));

initial equation

h_flueGas_out = h_start;
xi_flueGas = xi_start_flueGas_out;

equation

   mass = geo.volume * bulk.d;

  //_______________/ Composition of fuel and gas \_____________________
  xi_fuel_in = inStream(inlet.fuel.xi_outflow);
  xi_fuel_out =  xi_fuel_in;

  //________________/ Mass balance - flue gas \______________________________________
  inlet.flueGas.m_flow + outlet.flueGas.m_flow  =  drhodt*geo.volume;

  der(xi_flueGas) =  1/mass * (inlet.flueGas.m_flow*(flueGasInlet.xi - xi_flueGas) + outlet.flueGas.m_flow*(flueGasOutlet.xi - xi_flueGas));
  drhodt = bulk.drhodh_pxi * der(bulk.h) + sum({bulk.drhodxi_ph[i] * der(bulk.xi[i]) for i in 1:flueGas.nc-1});

  //______________ / Mass balance - Slag \____________________________________________________________________________
  0 = inlet.slag.m_flow  + outlet.slag.m_flow;

  //______________/ Mass balance - Fuel \____________________________
  0 = outlet.fuel.m_flow + inlet.fuel.m_flow;

  //_______________/ Energy Balance for gas \__________________________
  der(h_flueGas_out) = (Q_flow_wall + Q_flow_top + Q_flow_bottom
                + inlet.flueGas.m_flow * (flueGasInlet.h - h_flueGas_out)
                + outlet.flueGas.m_flow * (flueGasOutlet.h - h_flueGas_out)
                + inlet.slag.m_flow * (inlet.slagType.cp * (actualStream(inlet.slag.T_outflow) - 298.15) - h_flueGas_out)
                + outlet.slag.m_flow * (inlet.slagType.cp * (actualStream(outlet.slag.T_outflow)  - 298.15) - h_flueGas_out)
                + inlet.fuel.m_flow *(inStream(inlet.fuel.cp_outflow) * (inStream(inlet.fuel.T_outflow)  - 298.15) - h_flueGas_out)
                + outlet.fuel.m_flow * (outlet.fuel.cp_outflow * (outlet.fuel.T_outflow - 298.15) - h_flueGas_out))/mass;

  V_flow_flueGas_in = 0;
  V_flow_flueGas_out = 0;

  m_flow_in_del = inlet.flueGas.m_flow;
  m_flow_out_del = outlet.flueGas.m_flow;
  xi_flueGas_in_del = flueGasInlet.xi;
  T_bulk_del = bulk.T;
  rho_bulk_del = bulk.d;

  sum_xi = sum(flueGasOutlet.xi);

  xi_fuel = 0; // amount of fuel per flue gas mass

 //___________/ T_outflows \______________
  outlet.fuel.T_outflow = bulk.T;
  outlet.flueGas.T_outflow = bulk.T;
  heat_bottom.T = bulk.T;

  if slagTemperature_calculationType==1 then
    inlet.slag.T_outflow = T_Slag;
    outlet.slag.T_outflow = T_Slag;
  elseif slagTemperature_calculationType==2 then
    inlet.slag.T_outflow = flueGasOutlet.T;
    outlet.slag.T_outflow = flueGasOutlet.T;
  elseif slagTemperature_calculationType==3 then
    inlet.slag.T_outflow = (flueGasOutlet.T + flueGasInlet.T)/2;
    outlet.slag.T_outflow = (flueGasOutlet.T + flueGasInlet.T)/2;
  elseif slagTemperature_calculationType==4 then
    inlet.slag.T_outflow = flueGasInlet.T;
    outlet.slag.T_outflow = flueGasInlet.T;
  else
    inlet.slag.T_outflow = T_Slag;
    outlet.slag.T_outflow = T_Slag;
    assert(slagTemperature_calculationType==1 or slagTemperature_calculationType==2 or slagTemperature_calculationType==3 or slagTemperature_calculationType==4, "Invalid slag temperature calculation type");
  end if;

  inlet.flueGas.xi_outflow  = xi_flueGas_del;
  outlet.flueGas.xi_outflow  = xi_flueGas_del;

    //___________/ LHV_outflows \__________________________________________
  outlet.fuel.LHV_outflow =inStream(inlet.fuel.LHV_outflow);
  inlet.fuel.LHV_outflow =inStream(outlet.fuel.LHV_outflow);
  outlet.fuel.LHV_calculationType = inlet.fuel.LHV_calculationType;

  outlet.fuel.cp_outflow =inStream(inlet.fuel.cp_outflow);
  inlet.fuel.cp_outflow =inStream(outlet.fuel.cp_outflow);

  //___________/ Dummy T_outflows \__________________________________________
  inlet.fuel.T_outflow = bulk.T;
  //outlet.slag.T_outflow = inStream(outlet.slag.T_outflow); //outlet.slag is inflowing slag
  inlet.flueGas.T_outflow  = bulk.T;

  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-300,
            -100},{300,100}}),
                      graphics), Icon(coordinateSystem(preserveAspectRatio=true,
          extent={{-300,-100},{300,100}}), graphics={
        Text(
          extent={{32,80},{240,46}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          visible=showData,
          textString=DynamicSelect("", "T_out="+String(bulk.T,format="1.0f") +" K")),
        Text(
          extent={{32,20},{242,-14}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          visible=showData,
          textString=DynamicSelect("", "Q_wall="+String(Q_flow_wall/1e6,format="1.0f")+" MW"))}),
    Documentation(info="<html>
<p><b>Model description: </b>A nonstationary model for the hopper furnace section</p>
<p><br/><b>Contact:</b> Lasse Nielsen, TLK-Thermo GmbH</p>
<p><b>FEATURES</b> </p>
<p><ul>
<li>This model uses TILMedia</li>
<li>Nonstationary mass and energy balances are used</li>
<li>Different heat transfer correlations can be chosen</li>
</ul></p>
</html>"));
end Hopper_L2;
