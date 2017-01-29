within Exergy.XModelica.Media.Test;
model PsychrometricData1 "Produces plot data for psychrometric charts"
  extends Modelica.Icons.Example;

  package Medium = Exergy.XBuildings.Media.Air "Used medium package";

  Exergy.Utilities.RefEnv refEnv;

  parameter Modelica.SIunits.Pressure p_const=1e5 "Pressure";
  parameter Integer n_T=11 "Number of isotherms";
  parameter Modelica.SIunits.Temperature T_min=253.15 "Lowest isotherm";
  parameter Modelica.SIunits.Temperature T_step=10
    "Temperature step between two isotherms";
  parameter Integer n_h=16 "Number of lines with constant specific enthalpy";
  parameter Modelica.SIunits.SpecificEnthalpy h_min=-20e3
    "Lowest line of constant enthalpy";
  parameter Modelica.SIunits.SpecificEnthalpy h_step=1e4
    "Enthalpy step between two lines of constant enthalpy";
  parameter Integer n_phi=10 "Number of lines with constant relative humidity";
  parameter Real phi_min=0.1 "Lowest line of constant humidity";
  parameter Real phi_step=0.1 "Step between two lines of constant humidity";

  parameter Integer n_ex=1 "Number of lines with constant exergy";
  parameter Real ex_min=100 "Lowest line of constant exergy";
  parameter Real ex_step=1000 "Step between two lines of constant exergy";

  parameter Modelica.SIunits.MassFraction x_min=0.00
    "Minimum diagram absolute humidity";
  parameter Modelica.SIunits.MassFraction x_max=0.03
    "Maximum diagram absolute humidity";
  parameter Modelica.SIunits.Time t=1 "Simulation time";

  final parameter Modelica.SIunits.Temperature[n_T] T_const={T_min -
      T_step + i*T_step for i in 1:n_T} "Constant temperatures";
  final parameter Modelica.SIunits.SpecificEnthalpy[n_h] h_const={(i - 1)*
      h_step + h_min for i in 1:n_h} "Constant enthalpies";
  final parameter Real[n_phi] phi_const={(i - 1)*phi_step + phi_min for i in
          1:n_phi} "Constant relative humidities";

  final parameter Real[n_ex] ex_const={(i - 1)*ex_step + ex_min for i in
          1:n_ex} "Constant relative humidities";

  final parameter Real diagSlope=Medium.enthalpyOfVaporization(273.15)
    "Rotation of diagram that zero degrees isotherm becomes horizontal outside the fog region";
  final parameter Modelica.SIunits.MassFraction x_start=x_min
    "Initial absolute humidity in kg water/kg dry air";

  Modelica.SIunits.MassFraction x(start=x_start)
    "Absolute humidity in kg water/kg dry air";
  Modelica.SIunits.SpecificEnthalpy[n_T] hx_T "h_1+x for const T";
  Modelica.SIunits.SpecificEnthalpy[n_h] hx_h(start=h_const) "Const h_1+x";
  Modelica.SIunits.SpecificEnthalpy[n_phi] hx_phi "h_1+x for const phi";
    Modelica.SIunits.SpecificEnthalpy[n_ex] hx_ex "h_1+x for const phi";
  Modelica.SIunits.SpecificEnthalpy[n_T] y_T "Chart enthalpy for const T";
  Modelica.SIunits.SpecificEnthalpy[n_h] y_h "Chart enthalpy for const h";
  Modelica.SIunits.SpecificEnthalpy[n_phi] y_phi "Chart enthalpy for const phi";
    Modelica.SIunits.SpecificEnthalpy[n_ex] y_ex "Chart enthalpy for const phi";
  Medium.BaseProperties[n_T] medium_T "Medium properties for const T";
  Medium.BaseProperties[n_phi] medium_phi "Medium properties for const phi";
  Medium.BaseProperties[n_ex] medium_ex "Medium properties for const phi";
protected
  Modelica.SIunits.Pressure[n_phi] ps_phi
    "Saturation pressure for constant-phi-lines";
  Modelica.SIunits.Temperature[n_phi] T_phi(each start=290);
  Boolean[n_T] fog(start=fill(false, n_T))
    "Triggers events at intersection of isotherms with phi=1";
  Modelica.SIunits.Pressure[n_T] pd "Steam partial pressure along isotherms";
initial equation
  x = x_min;
equation

  der(x) = (x_max - x_min)/t;

  for i in 1:n_T loop
    medium_T[i].T = T_const[i];
    medium_T[i].p = p_const;
    medium_T[i].Xi = {x/(1 + x)};
    hx_T[i] = medium_T[i].h*(medium_T[i].Xi[1]/(1-medium_T[i].Xi[1]) + 1);
    y_T[i] = hx_T[i] - diagSlope*x;

    //trigger events
    pd[i] = medium_T[i].Xi[1]*medium_T[i].MM/Modelica.Media.Air.MoistAir.MMX[
      1]*p_const;
    fog[i] = pd[i] >= Medium.saturationPressure(T_const[i]);
  end for;
  for i in 1:n_h loop
    der(hx_h[i]) = 0.0;
    y_h[i] = hx_h[i] - diagSlope*x;
  end for;
  for i in 1:n_phi loop
    medium_phi[i].p = p_const;
    ps_phi[i] = p_const*x/phi_const[i]/(Medium.MMX[1]/Medium.MMX[2] + x);
    T_phi[i] = if x < 5e-6 then 200 else Exergy.XModelica.Media.Air.MoistAir.saturationTemperature(
      ps_phi[i]);
    medium_phi[i].T = T_phi[i];
    medium_phi[i].Xi = {x/(1 + x)};
    hx_phi[i] = medium_phi[i].h*(medium_T[i].Xi[1]/(1-medium_T[i].Xi[1]) + 1);
    y_phi[i] = hx_phi[i] - diagSlope*x;
  end for;

    for i in 1:n_ex loop
    medium_ex[i].p = p_const;
    Exergy.XModelica.Media.Air.MoistAir.specificExergy(medium_ex[i].state,refEnv)=ex_const[i];
    medium_ex[i].Xi = {x/(1 + x)};
    hx_ex[i] = medium_ex[i].h*(medium_T[i].Xi[1]/(1-medium_T[i].Xi[1]) + 1);
    y_ex[i] = hx_ex[i] - diagSlope*x;
  end for;

  annotation (Documentation(info="<html>
<p>This model produces psychrometric data from the moist air model in this library to be plotted in charts. The two most common chart varieties are the Mollier Diagram and the Psychrometric Chart. The first is widely used in some European countries while the second is more common in the Anglo-American world. Specific enthalpy is plotted over absolute humidity in the Mollier Diagram, it is the other way round in the Psychrometric Chart.<br>
It must be noted that the relationship of both axis variables is not right-angled, the absolute humidity follows a slope which equals the enthalpy of vaporization at 0 &deg;C. For better reading and in order to reduce the fog region the humidity axis is rotated to obtain a right-angled plot. Both charts usually contain additional information as isochores or auxiliary scales for e.g., heat ratios. Those information are omitted in this model and the charts below. Other important features of psychrometric chart data are that all mass specific variables (like absolute humidity, specific enthalpy etc.) are expressed in terms of kg dry air and that their baseline of 0 enthalpy is found at 0 &deg;C and zero humidity.</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Media/Air/Mollier.png\"><br>
<img src=\"modelica://Modelica/Resources/Images/Media/Air/PsycroChart.png\">
</p>

<p>
<b>Legend:</b> blue - constant specific enthalpy, red - constant temperature, black - constant relative humidity</p>

<p>The model provides data for lines of constant specific enthalpy, temperature and relative humidity in a Mollier Diagram or Psychrometric Chart as they were used for the figures above. For limitations and ranges of validity please refer to the <a href=\"modelica://Modelica.Media.Air.MoistAir\">MoistAir package description</a>. Absolute humidity <b>x</b> is increased with time in this model. The specific enthalpies adjusted for plotting are then obtained from:</p>
<ul>
<li><b>y_h</b>: constant specific enthalpy</li>
<li><b>y_T</b>: constant temperature</li>
<li><b>y_phi</b>: constant relative humidity</li>
</ul>
</html>"));
end PsychrometricData1;
