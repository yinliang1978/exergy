within Exergy.XClaRa.Components.TurboMachines.Pumps.Fundamentals;
model Outline
 extends ClaRa.Basics.Icons.RecordIcon;
 input SI.VolumeFlowRate V_flow "Volume flow rate";
 input SI.Power P_hyd "Hydraulic power";
 input SI.PressureDifference Delta_p "Pressure difference";
 input SI.Length head "Head";
 input SI.Length NPSHa "Net positive suction head available";
 input Real eta_hyd "Hydraulic efficiency";
 input Real eta_mech "Mechanic efficiency";

end Outline;
