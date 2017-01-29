within Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals;
record ICom
  extends ClaRa.Basics.Icons.IComIcon;
  import SI = ClaRa.Basics.Units;

//  SI.Pressure Delta_p "Pressure difference p_in - p_out" annotation(Dialog);
  SI.Pressure p_in "Inlet pressure" annotation(Dialog);
  SI.Pressure p_out "Outlet pressure" annotation(Dialog);
  Real opening_ "Opening of the valve (p.u.)" annotation(Dialog);
  SI.DensityMassSpecific rho_in "Inlet density" annotation(Dialog);
  Real gamma_in "Heat capacity ratio at inlet"   annotation(Dialog);
  Real gamma_out "Heat capacity ratio at outlet"    annotation(Dialog);
  Real opening_leak_ annotation(Dialog);
end ICom;
