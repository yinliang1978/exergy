within Exergy.XClaRa.SubSystems.Boiler;
partial model CoalSupplyBoiler_base "The coal mills and the boiler"
  extends ClaRa.Basics.Icons.Boiler;
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid   medium=simCenter.fluid1
    "Medium in the component"            annotation(choicesAllMatching=true, Dialog(group="Fundamental Definitions"));

  outer ClaRa.SimCenter simCenter;
  ClaRa.Basics.Interfaces.FluidPortOut livesteam(Medium=medium)
    annotation (Placement(transformation(extent={{-10,184},{10,204}}),
        iconTransformation(extent={{-10,210},{10,230}})));
   ClaRa.Basics.Interfaces.FluidPortOut reheat_out(Medium=medium)
     annotation (Placement(transformation(extent={{48,184},{68,204}}),
         iconTransformation(extent={{50,210},{70,230}})));
  Modelica.Blocks.Interfaces.RealInput QF_setl_
    "Set value of thermal output in p.u." annotation (Placement(transformation(
          extent={{-120,-20},{-80,20}}),iconTransformation(extent={{-140,-20},{
            -100,20}})));
  ClaRa.Basics.Interfaces.FluidPortIn feedwater(Medium=medium)
    annotation (Placement(transformation(extent={{-10,-110},{10,-90}}),
        iconTransformation(extent={{-10,-106},{10,-86}})));
   ClaRa.Basics.Interfaces.FluidPortIn reheat_in(Medium=medium)
     annotation (Placement(transformation(extent={{50,-110},{70,-90}}),
         iconTransformation(extent={{50,-106},{70,-86}})));
  annotation (Icon(coordinateSystem(extent={{-100,-100},{100,200}},
          preserveAspectRatio=false),
                   graphics), Diagram(coordinateSystem(extent={{-100,-100},{100,
            200}}, preserveAspectRatio=true),
                                      graphics));
end CoalSupplyBoiler_base;
