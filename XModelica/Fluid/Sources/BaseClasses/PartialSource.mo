within Exergy.XModelica.Fluid.Sources.BaseClasses;
partial model PartialSource
  extends Modelica.Fluid.Sources.BaseClasses.PartialSource(redeclare
      replaceable package Medium =
       Exergy.XBuildings.Media.Interfaces.PartialMedium);

  outer Exergy.Utilities.RefEnv refEnv "reference enviroment properties";

  //Utilities.ViewObject viewObject(nEnergy={nPorts,0,0,0});
 // Utilities.ViewPort viewOutput
  //  annotation (Placement(transformation(extent={{86,88},{106,108}})));
equation
/*
 for i in 1:nPorts loop
  viewObject.h[i].E_flow = ports[i].m_flow*noEvent(actualStream(ports[i].h_outflow));
  viewObject.h[i].Ex_flow = ports[i].m_flow*Medium.specificEnthalpyExergy(
                                   Medium.setState_phX(ports[i].p,
                                        noEvent(actualStream(ports[i].h_outflow)),
                                        noEvent(actualStream(ports[i].Xi_outflow))),
                                   refEnv.state,
                                   Exergy.Utilities.Types.ExTypes.EnthalpyExergy);
  end for;
  //heat

algorithm 
  viewObject.e[1].E_flow:=0;
  viewObject.e[1].Ex_flow:=0;
  viewObject.e[2].E_flow:=0;
  viewObject.e[2].Ex_flow:=0;

 for i in 1:nPorts loop
   if noEvent(ports[i].m_flow > 0) then
  viewObject.e[1].E_flow := viewObject.e[1].E_flow + viewObject.h[i].E_flow;
  viewObject.e[1].Ex_flow :=viewObject.e[1].Ex_flow + viewObject.h[i].Ex_flow;
   else
  viewObject.e[2].E_flow := viewObject.e[2].E_flow + viewObject.h[i].E_flow;
  viewObject.e[2].Ex_flow :=viewObject.e[2].Ex_flow + viewObject.h[i].Ex_flow;
   end if;
 end for;
*/
    //view out connect
 // connect(viewObject.viewOutput,viewOutput);
end PartialSource;
