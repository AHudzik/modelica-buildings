within Buildings.Utilities.Psychrometrics;
block Phi_TpWX
  extends BaseClasses.HumidityRatioVaporPressure;
  Modelica.SIunits.Pressure p_w;
  Modelica.Blocks.Interfaces.RealInput X_w(min=0, max=1, nominal=0.01)
                                                                      annotation (Placement(transformation(extent={{-140,-80},{-100,-40}})));

constant Real k = 0.621964713077499 "Ratio of molar masses";
Modelica.SIunits.AbsolutePressure psat "Saturation pressure";

 Modelica.Blocks.Interfaces.RealInput T(final unit="K",
                                           displayUnit="degC",
                                           min = 0) "Temperature"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));

  Modelica.Blocks.Interfaces.RealOutput phi(min = 0, max=1)
    "Relative humidity (0...1)"
                               annotation (Placement(transformation(extent={{100,-10},{120,10}})));

algorithm
  psat:=Buildings.HeatTransfer.Conduction.Functions.p_sat(T);
  p_w:=Buildings.Utilities.Psychrometrics.Functions.pW_X(X_w,p_in_internal);
  phi:= p_w/psat;
end Phi_TpWX;
