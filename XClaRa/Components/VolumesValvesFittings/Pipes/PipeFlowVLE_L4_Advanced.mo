within Exergy.XClaRa.Components.VolumesValvesFittings.Pipes;
model PipeFlowVLE_L4_Advanced
  "A 1D tube-shaped control volume considering one-phase and two-phase heat transfer in a straight pipe with detailed dynamic momentum and energy balance."
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
  extends ClaRa.Basics.ControlVolumes.FluidVolumes.VolumeVLE_L4_Advanced(
      redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.PipeGeometry_N_cv (
        z_in=z_in,
        z_out=z_out,
        N_tubes=N_tubes,
        N_cv=N_cv,
        diameter=diameter_i,
        length=length,
        Delta_x=Delta_x,
        N_passes=N_passes));
  extends ClaRa.Basics.Icons.Pipe_L4;
  ClaRa.Basics.Interfaces.Connected2SimCenter connected2SimCenter(
    powerIn=noEvent(if sum(heat.Q_flow) > 0 then sum(heat.Q_flow) else 0),
    powerOut=if not heatFlowIsLoss then -sum(heat.Q_flow) else 0,
    powerAux=0) if  contributeToCycleSummary;
//## P A R A M E T E R S #######################################################################################

//____Geometric data_____________________________________________________________________________________
  parameter SI.Length length=1 "Length of the pipe (one pass)"
    annotation (Dialog(group="Geometry"));
  parameter SI.Length diameter_i=0.1 "Inner diameter of the pipe"
    annotation (Dialog(group="Geometry"));
  parameter SI.Length z_in=0.1 "Height of inlet above ground"
    annotation (Dialog(group="Geometry"));
  parameter SI.Length z_out=0.1 "Height of outlet above ground"
    annotation (Dialog(group="Geometry"));

  parameter Integer N_tubes= 1 "Number Of parallel pipes"
                                                         annotation(Dialog(group="Geometry"));
  parameter Integer N_passes=1 "Number of passes of the tubes" annotation(Dialog(group="Geometry"));

//____Discretisation_____________________________________________________________________________________
public
  parameter Integer N_cv(min=3)=3 "Number of finite volumes" annotation(Dialog(group="Discretisation"));

  parameter SI.Length Delta_x[N_cv]=
      ClaRa.Basics.Functions.GenerateGrid(
                {0},
                length*N_passes,
                N_cv) "Discretisation scheme"
    annotation (Dialog(group="Discretisation"));
//________Summary_________________
  parameter Boolean contributeToCycleSummary = simCenter.contributeToCycleSummary
    "True if component shall contribute to automatic efficiency calculation"                annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean heatFlowIsLoss = true
    "True if negative heat flow is a loss (not a process product)"                                       annotation(Dialog(tab="Summary and Visualisation"));
//### E Q U A T I O N P A R T #######################################################################################
//-------------------------------------------

public
  ClaRa.Basics.Interfaces.EyeOut eye if showData annotation (
      Placement(transformation(extent={{140,-50},{160,-30}}),
        iconTransformation(extent={{136,-44},{156,-24}})));
protected
  ClaRa.Basics.Interfaces.EyeIn eye_int annotation (Placement(
        transformation(extent={{95,-41},{97,-39}})));

equation
//-------------------------------------------
//Summary:
  assert(abs(z_out-z_in) <= length, "Length of pipe less than vertical height", AssertionLevel.error);
  eye_int.m_flow=-outlet.m_flow;
  eye_int.T= fluidOutlet.T-273.15;
  eye_int.s=fluidOutlet.s/1e3;
  eye_int.p=outlet.p/1e5;
  eye_int.h=actualStream(outlet.h_outflow)/1e3;

//-------------------------------------------

  connect(eye,eye_int)  annotation (Line(
      points={{150,-40},{96,-40}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false,extent={{-140,-50},
            {140,50}}),
                   graphics={
        Polygon(
          points={{-132,42},{-122,42},{-114,34},{-114,-36},{-122,-42},{-132,-42},
              {-132,42}},
          pattern=LinePattern.None,
          smooth=Smooth.None,
          fillColor= {0,131,169},
          fillPattern=FillPattern.Solid,
          visible=frictionAtInlet),
        Polygon(
          points={{132,42},{122,42},{114,34},{114,-36},{122,-42},{132,-42},
              {132,42}},
          pattern=LinePattern.None,
          smooth=Smooth.None,
          fillColor= {0,131,169},
          fillPattern=FillPattern.Solid,
          visible=frictionAtOutlet)}),
                              Diagram(coordinateSystem(preserveAspectRatio=false,
          extent={{-140,-50},{140,50}}),
                                      graphics));
end PipeFlowVLE_L4_Advanced;
