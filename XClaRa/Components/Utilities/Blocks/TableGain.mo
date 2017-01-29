within Exergy.XClaRa.Components.Utilities.Blocks;
model TableGain "Table based gain"
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

  import Modelica.Blocks.Types;
  parameter Boolean tableOnFile=false
    "true, if table is defined on file or in function usertab"
    annotation(Dialog(group="table data definition"));
  parameter Real table[:, :]=fill(0.0,0,2)
    "table matrix (grid = first column; e.g., table=[0,2])"
       annotation(Dialog(group="table data definition", enable = not tableOnFile));
  parameter String tableName="NoName"
    "table name on file or in function usertab (see docu)"
       annotation(Dialog(group="table data definition", enable = tableOnFile));
  parameter String fileName="NoName" "file where matrix is stored"
       annotation(Dialog(group="table data definition", enable = tableOnFile,
                         __Dymola_loadSelector(filter="Text files (*.txt);;Matlab files (*.mat)",
                         caption="Open file in which table is present")));
  parameter Integer columns[:]=2:size(table, 2)
    "columns of table to be interpolated"
  annotation(Dialog(group="table data interpretation"));
  parameter Modelica.Blocks.Types.Smoothness smoothness=Types.Smoothness.LinearSegments
    "smoothness of table interpolation"
  annotation(Dialog(group="table data interpretation"));
  extends Modelica.Blocks.Interfaces.MIMOs(final n=size(columns, 1));

  parameter Boolean divide = false "True if y = u/table(u)";
protected
  Exergy.XClaRa.Components.Utilities.Blocks.ParameterizableTable1D table_block(
    tableOnFile=tableOnFile,
    table=table,
    tableName=tableName,
    fileName=fileName,
    columns=columns,
    smoothness=smoothness);
equation
  table_block.u = u;
  if divide then
    y = u ./ table_block.y;
  else
    y = u .* table_block.y;
  end if;
  annotation (
    Documentation(info="<html>

</HTML>
"), Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
                    graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={221,222,223},
          fillColor={118,124,127},
          fillPattern=FillPattern.Solid),
        Line(points={{-100,40},{-100,-40},{20,-40},{20,40},{-10,40},{-10,-40},{-70,-40},{-70,40},{-100,40},{-100,20},{20,20},{20,0},{-100,0},{-100,-20},{20,-20},{20,-40},{-100,-40},{-100,40},{20,40},{20,-40}},
                                                                      color={221,222,223}),
        Line(points={{-40,40},{-40,-40}},
                                      color={221,222,223}),
        Rectangle(
          extent={{-100,40},{-70,20}},
          lineColor={221,222,223},
          fillColor={235,183,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-100,20},{-70,0}},
          lineColor={221,222,223},
          fillColor={235,183,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-100,0},{-70,-20}},
          lineColor={221,222,223},
          fillColor={235,183,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-100,-20},{-70,-40}},
          lineColor={221,222,223},
          fillColor={235,183,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-100,100},{100,0},{-100,-100},{-100,100}},
          lineColor={221,222,223})}),
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}}), graphics={
        Rectangle(
          extent={{-60,60},{60,-60}},
          fillColor={235,235,235},
          fillPattern=FillPattern.Solid,
          lineColor={0,0,255}),
        Line(points={{-100,0},{-58,0}}, color={0,0,255}),
        Line(points={{60,0},{100,0}}, color={0,0,255}),
        Text(
          extent={{-100,100},{100,64}},
          textString="1 dimensional linear table interpolation",
          lineColor={0,0,255}),
        Line(points={{-54,40},{-54,-40},{54,-40},{54,40},{28,40},{28,-40},{-28,
              -40},{-28,40},{-54,40},{-54,20},{54,20},{54,0},{-54,0},{-54,-20},
              {54,-20},{54,-40},{-54,-40},{-54,40},{54,40},{54,-40}}, color={
              0,0,0}),
        Line(points={{0,40},{0,-40}}, color={0,0,0}),
        Rectangle(
          extent={{-54,40},{-28,20}},
          lineColor={0,0,0},
          fillColor={255,255,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-54,20},{-28,0}},
          lineColor={0,0,0},
          fillColor={255,255,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-54,0},{-28,-20}},
          lineColor={0,0,0},
          fillColor={255,255,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-54,-20},{-28,-40}},
          lineColor={0,0,0},
          fillColor={255,255,0},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-50,54},{-32,42}},
          textString="u[1]/[2]",
          lineColor={0,0,255}),
        Text(
          extent={{-24,54},{0,42}},
          textString="y[1]",
          lineColor={0,0,255}),
        Text(
          extent={{-2,-40},{30,-54}},
          textString="columns",
          lineColor={0,0,255}),
        Text(
          extent={{2,54},{26,42}},
          textString="y[2]",
          lineColor={0,0,255})}));
end TableGain;
