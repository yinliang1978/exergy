within Exergy.XClaRa.Components.Utilities.Blocks;
model ReadMatrixFromFile
  "Read a 2D matrix from file  || *.csv and *.mat(-v4) are supported"
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
  import ModelicaServices.ExternalReferences.loadResource;
  ////////// this is Dymola specific ///////////
  import DataFiles;
  ////////// MSL 3.2.2 conform Version / going to be used in the future ///////
  // import Modelica.Utilities.Streams.*;
  //////////////////////////////////////////////
  extends Modelica.Blocks.Interfaces.BlockIcon;

  parameter String fileName="NoName" "file where matrix is stored"
       annotation(Dialog(group="table data definition", enable = tableOnFile,
                         __Dymola_loadSelector(filter="Text files (*.txt);;Matlab files (*.mat)",
                         caption="Open file in which table is present")));
  parameter String matrixName="NoName"
    "table name on file or in function usertab (see docu)"
       annotation(Dialog(group="table data definition", enable = tableOnFile));

  final parameter Integer matrixSize[2]=readMatrixSize(loadResource(fileName), matrixName);
  ////////// this is Dymola specific ///////////
protected
  parameter Real M[:,:] = readMatrix(loadResource(fileName), matrixName, matrixSize[1], matrixSize[2]);
  ////////// MSL 3.2.2 conform Version / going to be used in the future ///////
  // parameter Real M[:,:] = readRealMatrix(loadResource(fileName), matrixName, matrixSize[1], matrixSize[2]);
  //////////////////////////////////////////////

public
  Modelica.Blocks.Interfaces.RealOutput y[matrixSize[1], matrixSize[2]]
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));
equation
 y=M;
  annotation (
    Documentation(info="<html>
<p>
<b>Linear interpolation</b> in <b>two</b> dimensions of a <b>table</b>.
The grid points and function values are stored in a matrix \"table[i,j]\",
where:
</p>
<ul>
<li> the first column \"table[2:,1]\" contains the u[1] grid points,</li>
<li> the first row \"table[1,2:]\" contains the u[2] grid points,</li>
<li> the other rows and columns contain the data to be interpolated.</li>
</ul>
<p>
Example:
</p>
<pre>
           |       |       |       |
           |  1.0  |  2.0  |  3.0  |  // u2
       ----*-------*-------*-------*
       1.0 |  1.0  |  3.0  |  5.0  |
       ----*-------*-------*-------*
       2.0 |  2.0  |  4.0  |  6.0  |
       ----*-------*-------*-------*
     // u1
   is defined as
      table = [0.0,   1.0,   2.0,   3.0;
               1.0,   1.0,   3.0,   5.0;
               2.0,   2.0,   4.0,   6.0]
   If, e.g. the input u is [1.0;1.0], the output y is 1.0,
       e.g. the input u is [2.0;1.5], the output y is 3.0.
</pre>
<ul>
<li> The interpolation is <b>efficient</b>, because a search for a new interpolation
     starts at the interval used in the last call.</li>
<li> If the table has only <b>one element</b>, the table value is returned,
     independent of the value of the input signal.</li>
<li> If the input signal <b>u1</b> or <b>u2</b> is <b>outside</b> of the defined <b>interval</b>,
     the corresponding value is also determined by linear
     interpolation through the last or first two points of the table.</li>
<li> The grid values (first column and first row) have to be <b>strict</b>
     monotonically increasing.</li>
</ul>
<p>
The table matrix can be defined in the following ways:
</p>
<ol>
<li> Explicitly supplied as <b>parameter matrix</b> \"table\",
     and the other parameters have the following values:
<pre>
   tableName is \"NoName\" or has only blanks,
   fileName  is \"NoName\" or has only blanks.
</pre></li>
<li> <b>Read</b> from a <b>file</b> \"fileName\" where the matrix is stored as
      \"tableName\". Both ASCII and binary file format is possible.
      (the ASCII format is described below).
      It is most convenient to generate the binary file from Matlab
      (Matlab 4 storage format), e.g., by command
<pre>
   save tables.mat tab1 tab2 tab3 -V4
</pre>
      when the three tables tab1, tab2, tab3 should be
      used from the model.</li>
<li>  Statically stored in function \"usertab\" in file \"usertab.c\".
      The matrix is identified by \"tableName\". Parameter
      fileName = \"NoName\" or has only blanks.</li>
</ol>
<p>
Table definition methods (1) and (3) do <b>not</b> allocate dynamic memory,
and do not access files, whereas method (2) does. Therefore (1) and (3)
are suited for hardware-in-the-loop simulation (e.g. with dSpace hardware).
When the constant \"NO_FILE\" is defined, all parts of the
source code of method (2) are removed by the C-preprocessor, such that
no dynamic memory allocation and no access to files takes place.
</p>
<p>
If tables are read from an ASCII-file, the file need to have the
following structure (\"-----\" is not part of the file content):
</p>
<pre>
-----------------------------------------------------
#1
double table2D_1(3,4)   # comment line
0.0  1.0  2.0  3.0  # u[2] grid points
1.0  1.0  3.0  5.0
2.0  2.0  4.0  6.0

double table2D_2(4,4)   # comment line
0.0  1.0  2.0  3.0  # u[2] grid points
1.0  1.0  3.0  5.0
2.0  2.0  4.0  6.0
3.0  3.0  5.0  7.0
-----------------------------------------------------
</pre>
<p>
Note, that the first two characters in the file need to be
\"#1\". Afterwards, the corresponding matrix has to be declared
with type, name and actual dimensions. Finally, in successive
rows of the file, the elements of the matrix have to be given.
Several matrices may be defined one after another.
The matrix elements are interpreted in exactly the same way
as if the matrix is given as a parameter. For example, the first
column \"table2D_1[2:,1]\" contains the u[1] grid points,
and the first row \"table2D_1[1,2:]\" contains the u[2] grid points.
</p>

</html>
"), Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
            100}}), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={221,222,223},
          fillColor={118,124,127},
          fillPattern=FillPattern.Solid),
        Line(points={{-60,40},{-60,-40},{60,-40},{60,40},{30,40},{30,-40},{-30,
              -40},{-30,40},{-60,40},{-60,20},{60,20},{60,0},{-60,0},{-60,-20},
              {60,-20},{60,-40},{-60,-40},{-60,40},{60,40},{60,-40}}, color={221,222,223}),
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
          extent={{-30,40},{0,20}},
          lineColor={221,222,223},
          fillColor={235,183,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{0,40},{30,20}},
          lineColor={221,222,223},
          fillColor={235,183,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{30,40},{60,20}},
          lineColor={221,222,223},
          fillColor={235,183,0},
          fillPattern=FillPattern.Solid),
        Line(points={{-60,40},{-30,20}}, color={221,222,223}),
        Line(points={{-30,40},{-60,20}}, color={221,222,223})}),
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
            100}}),     graphics={
        Rectangle(
          extent={{-60,60},{60,-60}},
          fillColor={235,235,235},
          fillPattern=FillPattern.Solid,
          lineColor={0,0,255}),
        Line(points={{60,0},{100,0}}, color={0,0,255}),
        Text(
          extent={{-100,100},{100,64}},
          textString="2 dimensional linear table interpolation",
          lineColor={0,0,255}),
        Line(points={{-54,40},{-54,-40},{54,-40},{54,40},{28,40},{28,-40},{-28,
              -40},{-28,40},{-54,40},{-54,20},{54,20},{54,0},{-54,0},{-54,-20},
              {54,-20},{54,-40},{-54,-40},{-54,40},{54,40},{54,-40}}, color={
              0,0,0}),
        Line(points={{0,40},{0,-40}}, color={0,0,0}),
        Rectangle(
          extent={{-54,20},{-28,0}},
          lineColor={0,0,0},
          fillColor={255,158,2},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-54,0},{-28,-20}},
          lineColor={0,0,0},
          fillColor={255,158,2},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-54,-20},{-28,-40}},
          lineColor={0,0,0},
          fillColor={255,158,2},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-28,40},{0,20}},
          lineColor={0,0,0},
          fillColor={255,158,2},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{0,40},{28,20}},
          lineColor={0,0,0},
          fillColor={255,158,2},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{28,40},{54,20}},
          lineColor={0,0,0},
          fillColor={255,158,2},
          fillPattern=FillPattern.Solid),
        Line(points={{-54,40},{-28,20}}, color={0,0,0}),
        Line(points={{-28,40},{-54,20}}, color={0,0,0}),
        Text(
          extent={{-54,-40},{-30,-56}},
          textString="u1",
          lineColor={0,0,255}),
        Text(
          extent={{28,58},{52,44}},
          textString="u2",
          lineColor={0,0,255}),
        Text(
          extent={{-2,12},{32,-22}},
          textString="y",
          lineColor={0,0,255})}));
end ReadMatrixFromFile;
