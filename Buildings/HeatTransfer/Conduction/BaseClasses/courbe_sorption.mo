within Buildings.HeatTransfer.Conduction.BaseClasses;
model courbe_sorption
extends Modelica.Blocks.Interfaces.SISO;

  replaceable parameter Data.BaseClasses.HygroThermalMaterial material
    "Material from Data.Solids, Data.SolidsPCM or Data.Resistances" annotation (
    Evaluate=true,
    choicesAllMatching=true,
    Placement(transformation(extent={{60,60},{80,80}})));

  Real numerator;
  Real denominator;
  Real b;

public
  Modelica.Blocks.Tables.CombiTable1D tab_sorption(
  final tableOnFile=false,
  final table=[material.sorp_tab_layer],
  final columns=2:2,
  final smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
equation
  if (material.Kunzel) then

        b              =  ( 0.8 * (material.w_80 - material.w_f)) / (material.w_80 - 0.8 * material.w_f);
        numerator      =  ( b -1)  * material.w_f * u;
        denominator    =  ( b - u);
        y              =  numerator / denominator;

  else

  connect(u, tab_sorption.u[1]) annotation (Line(
      points={{-120,0},{-12,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(tab_sorption.y[1], y) annotation (Line(
      points={{11,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  end if;

  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics), Icon(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
          Line(
          points={{-100,-100},{4,-96},{14,-94},{36,-84},{46,-76},{70,-32},{94,
              22},{98,74},{100,90}},
          color={0,0,0},
          smooth=Smooth.None)}));
end courbe_sorption;
