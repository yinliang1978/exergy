within Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals;
partial model GenericPressureLoss

  extends ClaRa.Basics.Icons.Delta_p;
  import SI = ClaRa.Basics.Units;
  parameter Real CL_valve[:, :]=[0, 0; 1, 1]
    "|Valve Characteristics|Effective apperture as function of valve position in p.u.";

  outer Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.ICom
    iCom;
  outer Boolean checkValve;
  outer Boolean useHomotopy;

  SI.MassFlowRate m_flow;
  Real aperture_ "Effective apperture";

  Modelica.Blocks.Tables.CombiTable1D ValveCharacteristics(table=CL_valve,
      columns={2});
  SI.PressureDifference Delta_p "Pressure difference p_in - p_out";
  Real gamma "Heat capacity ratio at actual inlet";
equation
  ValveCharacteristics.u[1] = noEvent(min(1, max(iCom.opening_, iCom.opening_leak_)));
  aperture_=noEvent(max(0,ValveCharacteristics.y[1]));
  Delta_p = iCom.p_in - iCom.p_out
  annotation (Icon(graphics));

end GenericPressureLoss;
