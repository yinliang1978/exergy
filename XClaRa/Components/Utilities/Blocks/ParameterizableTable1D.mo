within Exergy.XClaRa.Components.Utilities.Blocks;
model ParameterizableTable1D
  "Table look-up in one dimension (matrix/file) with n inputs and n outputs "
//___________________________________________________________________________//
// Component of the ClaRa library, version: 1.1.2                            //
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
//
// This is a copy of the MSL 3.2 CombiTable1D.
//
// If the value of the "table" parameter depends on values calculated in the initial equation system,
// OpenModelica and SimulationX do not suppport the usage of MSL 3.2.1 CombiTables.
//
// If this usecase is intended, ClaRa.Components.Utilities.Blocks.ParameterizableTable1D
// must be used instead of Modelica.Blocks.Tables.CombiTable1D,
// otherwise it will not be possible to simulate the model with SimulationX or OpenModelica.
//
// For a basic example see: ClaRa.Components.Utilities.Blocks.Check.TestParameterizableTable1D.
// For a realworld example see: ClaRa.Examples.ClaRaClosedLoop, which uses ClaRa.SubSystems.Boiler.SteamGenerator_L3
//
// This is a workaround. There is a request in the modelica ticket system, which asks to find a better solution:
// see https://trac.modelica.org/Modelica/ticket/2023
//
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

protected
  Integer tableID;

  function tableInit
    "Initialize 1-dim. table defined by matrix (for details see: Modelica/Resources/C-Sources/ModelicaTables.h)"
    input String tableName;
    input String fileName;
    input Real table[ :, :];
    input Modelica.Blocks.Types.Smoothness smoothness;
    output Integer tableID;
  external "C" tableID=  ModelicaTables_CombiTable1D_init(
                 tableName, fileName, table, size(table, 1), size(table, 2),
                 smoothness) annotation(Library="ModelicaExternalC");
  end tableInit;

  function tableIpo
    "Interpolate 1-dim. table defined by matrix (for details see: Modelica/Resources/C-Sources/ModelicaTables.h)"
    input Integer tableID;
    input Integer icol;
    input Real u;
    output Real value;
  external "C" value =
                     ModelicaTables_CombiTable1D_interpolate(tableID, icol, u) annotation(Library="ModelicaExternalC");
  end tableIpo;
equation
  if tableOnFile then
    assert(tableName<>"NoName", "tableOnFile = true and no table name given");
  end if;
  if not tableOnFile then
    assert(size(table,1) > 0 and size(table,2) > 0, "tableOnFile = false and parameter table is an empty matrix");
  end if;

  for i in 1:n loop
    y[i] = if not tableOnFile and size(table,1)==1 then
             table[1, columns[i]] else tableIpo(tableID, columns[i], u[i]);
  end for;
  when initial() then
    tableID=tableInit(if tableOnFile then tableName else "NoName",
                      if tableOnFile then fileName else "NoName", table, smoothness);
  end when;
  annotation (
    Documentation(info="<html>
This is a copy of the MSL 3.2 CombiTable1D.

If the value of the \"table\" parameter depends on values calculated in the initial equation system,
OpenModelica and SimulationX do not suppport the usage of MSL 3.2.1 CombiTables.

If this usecase is intended, \"ClaRa.Components.Utilities.Blocks.ParameterizableTable1D\"
should be used instead of \"Modelica.Blocks.Tables.CombiTable1D\",
as this enables the model to be simulated with SimulationX or OpenModelica.
</html>
"), Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
            100}}), graphics={     Rectangle(
          extent={{-100,-100},{100,100}},
          lineColor={0,0,127},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid), Text(
          extent={{-150,150},{150,110}},
          textString="%name",
          lineColor={0,0,255}),
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={221,222,223},
          fillColor={118,124,127},
          fillPattern=FillPattern.Solid),
        Line(points={{-60,40},{-60,-40},{60,-40},{60,40},{30,40},{30,-40},{-30,-40},{-30,40},{-60,40},{-60,20},{60,20},{60,0},{-60,0},{-60,-20},{60,-20},{60,-40},{-60,-40},{-60,40},{60,40},{60,-40}},
                                                                      color={221,222,223}),
        Line(points={{0,40},{0,-40}}, color={221,222,223}),
        Rectangle(
          extent={{-60,20},{-30,0}},
          lineColor={221,222,223},
          fillColor={235,183,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-60,0},{-30,-20}},
          lineColor={221,222,223},
          fillColor={235,183,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-60,-20},{-30,-40}},
          lineColor={221,222,223},
          fillColor={235,183,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-60,40},{-30,20}},
          lineColor={221,222,223},
          fillColor={235,183,0},
          fillPattern=FillPattern.Solid)}),
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
end ParameterizableTable1D;
