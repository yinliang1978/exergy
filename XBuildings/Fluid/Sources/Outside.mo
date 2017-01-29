within Exergy.XBuildings.Fluid.Sources;
model Outside
  extends Buildings.Fluid.Sources.Outside( replaceable package Medium =
     Exergy.XModelica.Media.PartialMedium);

       parameter Boolean exhaustLoss = false
    "Get the trace substances from the input connector"
    annotation(Evaluate=true);

  outer Exergy.Utilities.RefEnv refEnv "reference enviroment properties";

  Utilities.ViewObject viewObject(nEnergy={nPorts,2,0,0});
 Utilities.ViewPort viewOutput
   annotation (Placement(transformation(extent={{86,88},{106,108}})));

equation
 for i in 1:nPorts loop
  viewObject.h[i].E_flow = ports[i].m_flow*noEvent(actualStream(ports[i].h_outflow));
  viewObject.h[i].Ex_flow = ports[i].m_flow*Medium.specificExergy(
                                   Medium.setState_phX(ports[i].p,
                                        noEvent(actualStream(ports[i].h_outflow)),
                                        noEvent(actualStream(ports[i].Xi_outflow))),
                                   refEnv,
                                   Exergy.Utilities.Types.ExTypes.Flow);
  end for;
  //heat

algorithm

     viewObject.e[1].E_flow := 0;
  viewObject.e[1].Ex_flow :=0;
    viewObject.e[2].E_flow := 0;
  viewObject.e[2].Ex_flow :=0;
 for i in 1:nPorts loop
   if noEvent(ports[i].m_flow > 0) then
  viewObject.e[1].E_flow := viewObject.e[1].E_flow + viewObject.h[i].E_flow;
       if exhaustLoss == false then
       viewObject.e[1].Ex_flow :=viewObject.e[1].Ex_flow + viewObject.h[i].Ex_flow;
       else
       viewObject.e[1].Ex_flow :=viewObject.e[1].Ex_flow + 0;
       end if;
   else
  viewObject.e[2].E_flow := viewObject.e[2].E_flow + viewObject.h[i].E_flow;
  viewObject.e[2].Ex_flow :=viewObject.e[2].Ex_flow + viewObject.h[i].Ex_flow;
   end if;
 end for;

equation
    //view out connect
   connect(viewObject.viewOutput,viewOutput);
end Outside;
