within Exergy.XClaRa.Components.Furnace.FlameRoom;
model FlameRoomWithTubeBundle_L2_Dynamic
  "Model for a combustion chamber section with inner tube bundle heating surfaces"
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

extends Exergy.XClaRa.Components.Furnace.BaseClasses.CombustionChamberBase(
      redeclare replaceable model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.HollowBlockWithTubesAndCarrierTubes,
      geo(flowOrientation=ClaRa.Basics.Choices.GeometryOrientation.vertical));
extends
    Exergy.XClaRa.Components.Furnace.BaseClasses.CombustionChamberBase_additional_HPs(
      redeclare replaceable model HeatTransfer_CarrierTubes =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Gas_HT.Convection.Convection_finnedTubes_L2);
extends ClaRa.Basics.Icons.FlameRoomTubeBundle;

//import ClaRa.Basics.Functions.Stepsmoother;
//## P A R A M E T E R S #######################################################################################
inner parameter Boolean useHomotopy=simCenter.useHomotopy
    "True, if homotopy method is used during initialisation"  annotation(Dialog(tab="Initialisation"));

//## V A R I A B L E   P A R T##################################################################################
protected
ClaRa.Basics.Units.MassFraction xi_flueGas_in_del[flueGas.nc - 1]
    "Pseudo state for flue gas mixture composition";

Real drhodt "Density derivative";

ClaRa.Basics.Units.MassFlowRate m_flow_in_del
    "Pseudo state for inlet mass flow";
ClaRa.Basics.Units.MassFlowRate m_flow_out_del
    "Pseudo state for outlet mass flow";

//_____________________/ Media Objects \_________________________________
protected
TILMedia.Gas_pT     flueGasOutlet(p=outlet.flueGas.p, T= actualStream(outlet.flueGas.T_outflow),xi=actualStream(outlet.flueGas.xi_outflow),
        gasType=flueGas)
        annotation (Placement(transformation(extent={{-130,74},{-110,94}})));
public
    TILMedia.Gas_ph bulk(
    p(start=p_start_flueGas_out) = outlet.flueGas.p,
    xi=xi_flueGas_del,
    gasType=flueGas,
    h=h_flueGas_out_del)
      annotation (Placement(transformation(extent={{-130,26},{-110,46}})));

//___________________/ iCom record \\__________________
protected
  inner ClaRa.Basics.Records.IComGas_L2 iCom(
    m_flow_nom=m_flow_nom,
    T_bulk=bulk.T,
    p_bulk=bulk.p,
    fluidPointer_in=flueGasInlet.gasPointer,
    fluidPointer_bulk=bulk.gasPointer,
    fluidPointer_out=flueGasOutlet.gasPointer,
    mediumModel=flueGas,
    p_in=inlet.flueGas.p,
    T_in=flueGasInlet.T,
    m_flow_in=m_flow_in_del,
    V_flow_in=V_flow_flueGas_in,
    xi_in=xi_flueGas_in_del,
    p_out=outlet.flueGas.p,
    T_out=bulk.T,
    m_flow_out=m_flow_out_del,
    V_flow_out=V_flow_flueGas_out,
    xi_out=xi_flueGas_del) annotation (Placement(transformation(extent={{244,-102},{268,-76}})));

//___________________/ Summary \\__________________
public
  Summary summary(
    outline(
      volume=geo.volume,
      A_cross=geo.A_front,
      A_wall=if geo.flowOrientation == ClaRa.Basics.Choices.GeometryOrientation.vertical then 2*(geo.width + geo.length)*geo.height else 2*(geo.width + geo.height)*geo.length,
      height=geo.height,
      m=mass,
      m_flow_fuel_burned=m_flow_fuel_burned,
      m_flow_oxygen_burned=m_flow_oxygen_burned,
      m_flow_oxygen_req=m_flow_oxygen_req,
      m_flow_air_req=m_flow_air_req,
      lambdaComb=lambdaComb,
      NOx_fraction=reactionZone.xi_NOx,
      CO_fraction=reactionZone.xi_CO,
      LHV=LHV,
      Q_combustion=m_flow_fuel_burned*LHV,
      w_migration=particleMigration.w,
      t_dwell_flueGas=t_dwell_flueGas,
      burning_time=burning_time.t,
      unburntFraction=unburntFraction,
      T_out=flueGasOutlet.T,
      h_out=flueGasOutlet.h),
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

initial equation

  h_flueGas_out = h_start;
  xi_flueGas = xi_start_flueGas_out;

equation

  if noEvent(t_dwell_flueGas < burning_time.t) then
    unburntFraction = (1.0 - t_dwell_flueGas/burning_time.t);
  else
    unburntFraction = 0;
  end if;

  mass = geo.volume * (bulk.d + flueGasInlet.d)/2;

  //____________/ Resulting Xi for entire fuel mass in the volume \______________
  xi_fuel_in = inStream(inlet.fuel.xi_outflow);

  //________________/ Mass balance - flue gas \______________________________________
  drhodt*geo.volume =m_flow_fuel_burned*(1 - xi_fuel_in[6]*reactionZone.xi_slag) + inlet.flueGas.m_flow + outlet.flueGas.m_flow;
  drhodt = bulk.drhodh_pxi * der(bulk.h) + sum({bulk.drhodxi_ph[i] * der(bulk.xi[i]) for i in 1:flueGas.nc-1});
  //______________ / Mass balance - Slag \____________________________________________________________________________
  0 =inlet.slag.m_flow + m_flow_fuel_burned*xi_fuel_in[6]*reactionZone.xi_slag + outlet.slag.m_flow;

  //______________/ Mass balance - Fuel \____________________________
  0 =outlet.fuel.m_flow + inlet.fuel.m_flow - m_flow_fuel_burned;

  //__________/ molar flow rates of combustable components (educts) into the whole burner system (maybe not all of it is burned) \________
  n_flow_C = xi_fuel_in[1]*inlet.fuel.m_flow /ClaRa.Basics.Constants.M_C;
  n_flow_H = xi_fuel_in[2]*inlet.fuel.m_flow /ClaRa.Basics.Constants.M_H;
  n_flow_O = xi_fuel_in[3]*inlet.fuel.m_flow /ClaRa.Basics.Constants.M_O;
  n_flow_N = xi_fuel_in[4]*inlet.fuel.m_flow /ClaRa.Basics.Constants.M_N;
                                                             // N not N2!!!!!!
  n_flow_S = xi_fuel_in[5]*inlet.fuel.m_flow /ClaRa.Basics.Constants.M_S;
  n_flow_Ash = xi_fuel_in[6]*inlet.fuel.m_flow /ClaRa.Basics.Constants.M_Ash;
  n_flow_H2O = (1-sum(xi_fuel_in))*inlet.fuel.m_flow /ClaRa.Basics.Constants.M_H2O;

  //_______________/ determination of lambda \_________________________
  // theoretically required oxygen mass flow rate to burn all the fuel
   m_flow_oxygen_req = (1-unburntFraction)*(n_flow_C + n_flow_H/4.0 + n_flow_S - n_flow_O/2)*ClaRa.Basics.Constants.M_O
                                                                                          *2.0;
   m_flow_air_req*max(1e-32,flueGasInlet.xi[6]) = m_flow_oxygen_req;

     if noEvent(m_flow_oxygen_req <= 0) then
    lambdaComb = 1.0e3;
     else
    lambdaComb = (inlet.flueGas.m_flow*flueGasInlet.xi[6])/max(1e-12, m_flow_oxygen_req);
     end if;

  //calculation of actual fuel and oxygen mass flow rates that are burned
  if noEvent(lambdaComb > 1) then
    m_flow_fuel_burned = (1 - unburntFraction)*inlet.fuel.m_flow;
    m_flow_oxygen_burned = m_flow_oxygen_req;
   else
    m_flow_fuel_burned = lambdaComb*(1 - unburntFraction)*inlet.fuel.m_flow;
    m_flow_oxygen_burned = lambdaComb*m_flow_oxygen_req;
   end if;

  //_____________/ Calculation of the LHV \______________________________________

   if inlet.fuel.LHV_calculationType == "predefined" then
     LHV = inStream(inlet.fuel.LHV_outflow);
   elseif inlet.fuel.LHV_calculationType == "Verbandsformel" then
     LHV =(33907*xi_fuel_in[1] + 142324*(xi_fuel_in[2] - xi_fuel_in[3]/8.) + 10465*xi_fuel_in[5] - 2512*((1 - sum(xi_fuel_in)) + 9*xi_fuel_in[2]))*1000;
   else
     LHV = inStream(inlet.fuel.LHV_outflow);
   end if;

   cp = inStream(inlet.fuel.cp_outflow);

//   if (inlet.fuel.m_flow < 0 or inlet.fuel.m_flow > 0) or (outlet.fuel.m_flow < 0 or outlet.fuel.m_flow > 0) then
//      LHV * (max(0,inlet.fuel.m_flow) + max(0,outlet.fuel.m_flow)) = (max(0,inlet.fuel.m_flow)*inStream(inlet.fuel.LHV_outflow) +  max(0,outlet.fuel.m_flow)*inStream(outlet.fuel.LHV_outflow));
//   else
//      LHV = inStream(inlet.fuel.LHV_outflow);
//   end if;

  //______________________________/ mass balance of flue gas components \__________________________
  der(xi_flueGas) =1/mass*(inlet.flueGas.m_flow*(flueGasInlet.xi - xi_flueGas) + outlet.flueGas.m_flow*(flueGasOutlet.xi - xi_flueGas) + m_flow_fuel_burned*(reactionZone.prod_comp - xi_flueGas));

  //_____________/ Calculation of fuel formation enthalpy with LHV for an ideal combustion\__________________
  m_flow_fuel_id = 1.0;
  m_flow_flueGas_id =(m_flow_fuel_id*(1 - xi_fuel_in[6]*reactionZone.xi_slag));           //ideal flue gas mass flow
   xi_flueGas_id =1/m_flow_flueGas_id*reactionZone.prod_comp;   //products of an ideal combustion

   sum_comp = sum(xi_flueGas_id);
  Delta_h_f - LHV =m_flow_flueGas_id*((ideal_combustion.h_i)*cat(
    1,
    xi_flueGas_id,
    {1 - sum(xi_flueGas_id)})) + xi_fuel_in[6]*reactionZone.xi_slag*outlet.slagType.cp*T_0;                  //formation enthalpy of used fuel

  //_______________/ Energy Balance flueGasCombustion \__________________________
  der(h_flueGas_out) =(Q_flow_wall + Q_flow_top + Q_flow_bottom + Q_flow_CarrierTubes + Q_flow_TubeBundle + inlet.flueGas.m_flow*(flueGasInlet.h - h_flueGas_out) + inlet.fuel.m_flow*((inStream(inlet.fuel.cp_outflow)*(inStream(inlet.fuel.T_outflow) - T_0) + Delta_h_f) - h_flueGas_out) + outlet.fuel.m_flow*((cp*(outlet.fuel.T_outflow - T_0) + Delta_h_f) - h_flueGas_out) + outlet.slag.m_flow*(outlet.slagType.cp*(actualStream(outlet.slag.T_outflow) - T_0) - h_flueGas_out) + inlet.slag.m_flow*(inlet.slagType.cp*(actualStream(inlet.slag.T_outflow) - T_0) - h_flueGas_out) + outlet.flueGas.m_flow*(flueGasOutlet.h - h_flueGas_out))/mass;

  //______________/Properties for heat transfer corellation\_________

  V_flow_flueGas_out = outlet.flueGas.m_flow/flueGasOutlet.d;
  V_flow_flueGas_in = inlet.flueGas.m_flow/flueGasInlet.d;

  m_flow_in_del = inlet.flueGas.m_flow;
  m_flow_out_del = outlet.flueGas.m_flow;

  xi_flueGas_in_del =flueGasInlet.xi;

  sum_xi = sum(flueGasOutlet.xi);

  xi_fuel_out = xi_fuel_in; //no change of fuel composition during combustion

  xi_fuel = (inlet.fuel.m_flow)/(inlet.flueGas.m_flow);// amount of fuel per flue gas mass

  //___________/ T_outflows \__________________________________________
  outlet.fuel.T_outflow = bulk.T;
  outlet.flueGas.T_outflow = bulk.T;
  //outlet.slag.T_outflow = inStream(outlet.slag.T_outflow); //outlet.slag is inflowing slag

  heat_bottom.T = bulk.T;

  inlet.flueGas.T_outflow  = bulk.T;
  inlet.fuel.T_outflow = bulk.T;

  if slagTemperature_calculationType==1 then
    inlet.slag.T_outflow = T_slag;
    outlet.slag.T_outflow = T_slag;
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
    inlet.slag.T_outflow = T_slag;
    outlet.slag.T_outflow = T_slag;
    assert(slagTemperature_calculationType==1 or slagTemperature_calculationType==2 or slagTemperature_calculationType==3 or slagTemperature_calculationType==4, "Invalid slag temperature calculation type");
  end if;

  //___________/ LHV_outflows \__________________________________________
  outlet.fuel.LHV_outflow =LHV;
  inlet.fuel.LHV_outflow =LHV;
  outlet.fuel.LHV_calculationType = inlet.fuel.LHV_calculationType;

  outlet.fuel.cp_outflow =cp;
  inlet.fuel.cp_outflow=inStream(outlet.fuel.cp_outflow);

  inlet.flueGas.xi_outflow  = xi_flueGas_del;
  outlet.flueGas.xi_outflow  = xi_flueGas_del;

  //______________Eye port variable definition________________________
  eye_int.m_flow = -outlet.flueGas.m_flow;
  eye_int.T = flueGasOutlet.T-273.15;
  eye_int.s = flueGasOutlet.s/1e3;
  eye_int.p = flueGasOutlet.p/1e5;
  eye_int.h = flueGasOutlet.h/1e3;
  eye_int.xi = flueGasOutlet.xi;

  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-300,-100},
            {300,100}}),
                      graphics), Icon(coordinateSystem(preserveAspectRatio=true,
          extent={{-300,-100},{300,100}}), graphics={
        Text(
          extent={{32,74},{240,44}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          visible=showData,
          textString=DynamicSelect("", "T_out="+String(bulk.T,format="1.0f") +" K")),
        Text(
          extent={{32,-6},{240,-36}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          visible=showData,
          textString=DynamicSelect("", "Q_wall="+String(Q_flow_wall/1e6,format="1.0f")+" MW")),
        Text(
          extent={{32,34},{240,4}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          visible=showData,
          textString=DynamicSelect("", "lambda="+String(min(99,lambdaComb),format="1.1f"))),
        Text(
          extent={{32,-46},{240,-76}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          visible=showData,
          textString=DynamicSelect("", "Q_bundle="+String(Q_flow_TubeBundle/1e6,format="1.0f")+" MW"))}),
    experiment(StopTime=10),
    __Dymola_experimentSetupOutput,
    Documentation(info="<html>
<p><b>Model description: </b>A stationary model for furnace sections with convective tube bundles inside (additional heat ports)</p>
<p><b>Contact:</b> Lasse Nielsen, TLK-Thermo GmbH</p>
<p><b>FEATURES</b> </p>
<p><ul>
<li>This model uses TILMedia</li>
<li>Stationary mass and energy balances are used</li>
<li>Tube bundles are regarded with replaceable geometry models</li>
<li>Additional heat transfer ports to regard heat transfer for tube banks and carrier tubes</li>
<li>The formation enthalpy of the used fuel is calculated with the given Lower heating value and and ideal combustion with the given elemental composition of the fuel</li>
<li>Lower heating can be regarded with a fixed value or calculated according to the &QUOT;Verbandsformel&QUOT;</li>
<li>Different heat transfer correlations can be chosen</li>
<li>Capable to burn unburnt fuel from lower burner sections</li>
<li>Amount of burnable fuel mass is calculated with particle diffusity which depends on models used to determine the mean migration speed according to the volume flow rates of the flue gas</li>
</ul></p>
</html>"));
end FlameRoomWithTubeBundle_L2_Dynamic;
