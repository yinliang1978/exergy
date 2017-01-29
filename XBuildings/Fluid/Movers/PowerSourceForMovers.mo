within Exergy.XBuildings.Fluid.Movers;
model PowerSourceForMovers
 // input Real Power;
   outer Exergy.Utilities.RefEnv refEnv "reference enviroment properties";

  Utilities.ViewObject viewObject(nEnergy={0,1,1,0})
    annotation (Placement(transformation(extent={{40,60},{60,80}})));
  Utilities.ViewPort viewPort
    annotation (Placement(transformation(extent={{86,84},{106,104}})));
  Modelica.Blocks.Interfaces.RealInput Power
    annotation (Placement(transformation(extent={{-120,56},{-80,96}})));
equation

     viewObject.e[1].E_flow =  -Power;
  viewObject.e[1].Ex_flow = -Power;
    viewObject.w[1].E_flow =  -Power;
  viewObject.w[1].Ex_flow = -Power;
  connect(viewObject.viewOutput, viewPort)
    annotation (Line(points={{60,70},{96,70},{96,94}}, color={28,108,200}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})), Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}), graphics={Line(
          points={{-76,-2},{-56,48},{-34,0},{-10,-40},{-4,-40},{12,0},{32,40},{42,
              26},{60,0}},
          color={0,0,0},
          smooth=Smooth.Bezier)}));
end PowerSourceForMovers;
