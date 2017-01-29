within Exergy.XClaRa.Components.Adapters;
model Scalar2VectorHeatPort
  "Connect a scalar heat port with a vectorised heat port"
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
  extends ClaRa.Basics.Icons.Adapter5_fw;
  import SI = ClaRa.Basics.Units;

  parameter String equalityMode = "Equal Temperatures"
    "Spacial equality of state or flow variable?"                                                       annotation(Dialog(group="Fundamental Definitions"),choices(choice="Equal Heat Flow Rates",  choice= "Equal Temperatures"));

  parameter Integer N = 3 "Number of axial elements" annotation(Dialog(group="Fundamental Definitions"));

  parameter SI.Length length=1 "Length of adapter" annotation (Dialog(
        group="Discretisation", enable=(equalityMode ==
          "Equal Heat Flow Rates")));
 parameter SI.Length Delta_x[N]=ClaRa.Basics.Functions.GenerateGrid(
              {0},
              length,
              N) "Discretisation scheme" annotation (Dialog(group=
          "Discretisation", enable=(equalityMode ==
          "Equal Heat Flow Rates")));

  parameter Boolean useStabiliserState= false
    "True, if a stabiliser state shall be used"                                           annotation(Dialog(tab="Expert Settings",enable = (equalityMode =="Equal Temperatures")));
  parameter ClaRa.Basics.Units.Time Tau=1 "Time Constant of Stabiliser State" annotation(Dialog(tab="Expert Settings",enable = (equalityMode =="Equal Temperatures") and (useStabiliserState)));

  ClaRa.Basics.Interfaces.HeatPort_a
                                   heatScalar
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
  ClaRa.Basics.Interfaces.HeatPort_b heatVector[N]
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));
equation

  if equalityMode ==  "Equal Heat Flow Rates" then
    heatScalar.T = sum(heatVector.T.*Delta_x)/sum(Delta_x);
    heatVector.Q_flow = -heatScalar.Q_flow.*Delta_x/sum(Delta_x);
  elseif equalityMode == "Equal Temperatures" then
    heatScalar.Q_flow = -sum(heatVector.Q_flow);
    if useStabiliserState then
       der(heatVector.T) = (ones(N).*heatScalar.T - heatVector.T)/Tau;
    else
      heatVector.T = ones(N).*heatScalar.T;
    end if;
  else
    assert(false, "Unknown equalityMode option in scalar2VectorHeatPort");
  end if;

initial equation
if equalityMode == "Equal Temperatures" and useStabiliserState ==true then
  heatVector.T= ones( N)*heatScalar.T;
end if;
  annotation (Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}}),
                   graphics), Diagram);
end Scalar2VectorHeatPort;
