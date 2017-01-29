within Exergy.XClaRa.Components.VolumesValvesFittings.Pipes;
model PipeFlowGas_L4_Simple
  "A 1D tube-shaped control volume considering heat transfer in a straight pipe with static momentum balance and simple energy balance."
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

  extends ClaRa.Basics.ControlVolumes.GasVolumes.VolumeGas_L4(
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.PipeGeometry_N_cv (
        N_tubes=N_tubes,
        N_cv=N_cv,
        diameter=diameter_i,
        length=length,
        Delta_x=Delta_x));

  extends ClaRa.Basics.Icons.Pipe_L4;
  extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L4");
  ClaRa.Basics.Interfaces.Connected2SimCenter connected2SimCenter(
    powerIn=noEvent(if sum(heat.Q_flow) > 0 then sum(heat.Q_flow) else 0),
    powerOut=if not heatFlowIsLoss then -sum(heat.Q_flow) else 0,
    powerAux=0) if  contributeToCycleSummary;

//## P A R A M E T E R S #######################################################################################

//____Geometric data_____________________________________________________________________________________
  parameter SI.Length length=1 "|Geometry|Length of the pipe";
  parameter SI.Length diameter_i=0.1 "|Geometry|Inner diameter of the pipe";

  parameter Integer N_tubes= 1 "|Geometry|Number Of parallel pipes";

//____Discretisation_____________________________________________________________________________________
    parameter Integer N_cv(min=3)=3 "|Discretisation|Number of finite volumes";
public
  inner parameter SI.Length Delta_x[N_cv]=
      ClaRa.Basics.Functions.GenerateGrid(
                {0},
                length,
                N_cv) "|Discretisation|Discretisation scheme";

//________Summary_________________
  parameter Boolean contributeToCycleSummary = simCenter.contributeToCycleSummary
    "True if component shall contribute to automatic efficiency calculation"                annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean heatFlowIsLoss = true
    "True if negative heat flow is a loss (not a process product)"                                       annotation(Dialog(tab="Summary and Visualisation"));

protected
  ClaRa.Basics.Interfaces.EyeIn eye_int annotation (Placement(
        transformation(extent={{85,-41},{87,-39}})));
public
  ClaRa.Basics.Interfaces.EyeOut eye if showData annotation (
      Placement(transformation(extent={{130,-50},{150,-30}}),
        iconTransformation(extent={{136,-44},{156,-24}})));

//### E Q U A T I O N P A R T #######################################################################################
//-------------------------------------------
equation

  eye_int.m_flow=-outlet.m_flow;
  eye_int.T= fluidOutlet.T-273.15;
  eye_int.s=fluidOutlet.s/1e3;
  eye_int.p=outlet.p/1e5;
  eye_int.h=fluidOutlet.h/1e3;
         //fillColor={0,131,169};//DynamicSelect(if time > 0 then (if not FlowModel==FlowModelStructure.inlet_innerPipe_outlet and not FlowModel==FlowModelStructure.inlet_innerPipe_dp_outlet then {0,131,169} else {255,255,255}) else {255,255,255}),
  connect(eye_int,eye)  annotation (Line(
      points={{86,-40},{140,-40}},
      color={255,204,51},
      smooth=Smooth.None,
      thickness=0.5));
annotation (Icon(coordinateSystem(preserveAspectRatio=false,extent={{-140,-50},
            {140,50}}),        graphics={
        Polygon(
          points={{-132,42},{-122,42},{-114,34},{-114,-36},{-122,-42},{-132,-42},
              {-132,42}},
          pattern=LinePattern.None,
          smooth=Smooth.None,
          fillColor= {118,106,98},
          fillPattern=FillPattern.Solid,
          visible=frictionAtInlet),
        Polygon(
          points={{132,42},{122,42},{114,34},{114,-36},{122,-42},{132,-42},
              {132,42}},
          pattern=LinePattern.None,
          smooth=Smooth.None,
          fillColor= {118,106,98},
          fillPattern=FillPattern.Solid,
          visible=frictionAtOutlet)}), Diagram(coordinateSystem(preserveAspectRatio=false,
          extent={{-140,-50},{140,50}})),
    Documentation(info="<html>
<p><b>Model description: </b>A non-adiabatic 1D-tube model using a single pipe cell for the formulation</p>
<p><b>Contact:</b> Johannes Brunnemann, XRG Simulation GmbH</p>
<p>
<b>FEATURES</b>
<ul>
<li>This model uses TILMedia</li>
<li>Flow reversal is supported</li>

<li>distributed pressure loss, i.e. pressure loss occurs in first and second half cell, whereas the state is located in the cell center</li>
</ul></p>
<b>TODO</b>
<ul>

</ul>


<h4>Staggered Grid Approach</h4>


<p>

</p>

<h4>State Definitions</h4>
<p>

</p>


</html>"));
end PipeFlowGas_L4_Simple;
