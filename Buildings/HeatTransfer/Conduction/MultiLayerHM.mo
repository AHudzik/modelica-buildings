within Buildings.HeatTransfer.Conduction;
model MultiLayerHM
  "Model for heat and moisture conductance through a solid with multiple material layers"

  Modelica.SIunits.Temperature T[sum(nSta)](each nominal = 300)
    "Temperature at the states";
  Modelica.SIunits.HeatFlowRate Q_flow[sum(nSta)+nLay]
    "Heat flow rate from state i to i+1";
  Modelica.SIunits.MassFlowRate m_flow[sum(nSta)+nLay]
    "Mass flow rate from state i to i+1";
  Real phi[sum(nSta)];
  Modelica.SIunits.MassConcentration w[sum(nSta)];

  parameter Integer switch_w "switch for the water content"
    annotation (Dialog(tab="HM"));
  parameter Integer switch_lamb "switch for the heat conductivity"
    annotation (Dialog(tab="HM"));
  parameter Integer switch_dw
    "switch for the liquid transport coefficient for suction "
    annotation (Dialog(tab="HM"));

parameter Boolean activatesuction;

  extends Buildings.HeatTransfer.Conduction.BaseClasses.PartialConstructionHM;

protected
  SingleLayerHM[nLay] lay(each final A=A, material=layers.material, switch_w=switch_w, switch_lamb=switch_lamb, switch_dw=switch_dw, activatesuction=activatesuction)
    "Material layer"
    annotation (Placement(transformation(extent={{-20,-12},{0,12}})));

public
  Interfaces.HeatMassPort_a heatMassPort_a
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
  Interfaces.HeatMassPort_b heatMassPort_b
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));
equation
 // This section assigns the temperatures and heat flow rates of the layer models to
  // an array that makes plotting the results easier.
  for i in 1:nLay loop
    for j in 1:nSta[i] loop
      T[sum(nSta[k] for k in 1:(i-1)) +j] = lay[i].T[j];
      phi[sum(nSta[k] for k in 1:(i-1)) +j] = lay[i].phi[j];
      w[sum(nSta[k] for k in 1:(i-1)) +j] = lay[i].w[j];
    end for;
    for j in 1:nSta[i]+1 loop
      Q_flow[sum(nSta[k] for k in 1:i-1)+(i-1)+j] = lay[i].Q_flow[j];
      m_flow[sum(nSta[k] for k in 1:i-1)+(i-1)+j] = lay[i].m_flow[j];
    end for;
  end for;
   connect(heatMassPort_a, lay[1].heatMassPort_a) annotation (Line(
      points={{-100,5.55112e-16},{-60,5.55112e-16},{-60,6.10623e-16},{-20,
          6.10623e-16}},
      color={191,0,0},
      smooth=Smooth.None));
  for i in 1:nLay-1 loop
  connect(lay[i].heatMassPort_b, lay[i+1].heatMassPort_a) annotation (Line(
      points={{5.55112e-16,6.10623e-16},{20,6.10623e-16},{20,-20},{-40,-20},{
            -40,6.10623e-16},{-20,6.10623e-16}},
      color={191,0,0},
      smooth=Smooth.None));
  end for;
  connect(lay[nLay].heatMassPort_b, heatMassPort_b) annotation (Line(
      points={{5.55112e-16,6.10623e-16},{49,6.10623e-16},{49,5.55112e-16},{100,
          5.55112e-16}},
      color={191,0,0},
      smooth=Smooth.None));

  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics), Icon(graphics={
        Rectangle(
          extent={{-96,4},{96,-4}},
          lineColor={0,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.HorizontalCylinder),
        Rectangle(
          extent={{-70,52},{-58,-50}},
          lineColor={0,0,0},
          fillColor={192,192,192},
          fillPattern=FillPattern.Backward),
        Polygon(
          points={{-26,16},{-24,16},{-22,12},{-20,6},{-20,2},{-22,-4},{-28,-8},
              {-32,-12},{-40,-14},{-44,-10},{-50,-4},{-52,6},{-50,12},{-48,16},
              {-46,18},{-44,20},{-40,22},{-36,22},{-30,20},{-26,16}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-44,-8},{-36,-10},{-30,-8},{-24,-6},{-28,-8},{-32,-12},{-40,
              -14},{-46,-12},{-50,-4},{-52,6},{-50,12},{-48,16},{-46,18},{-48,8},
              {-48,0},{-44,-8}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-10,52},{2,-50}},
          lineColor={0,0,0},
          fillColor={192,192,192},
          fillPattern=FillPattern.Backward),
        Rectangle(
          extent={{14,52},{26,-50}},
          lineColor={0,0,0},
          fillColor={95,95,95},
          fillPattern=FillPattern.Forward),
        Polygon(
          points={{58,14},{60,14},{62,10},{64,4},{64,0},{62,-6},{56,-10},{52,
              -14},{44,-16},{40,-12},{34,-6},{32,4},{34,10},{36,14},{38,16},{40,
              18},{44,20},{48,20},{54,18},{58,14}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{40,-10},{48,-12},{54,-10},{60,-8},{56,-10},{52,-14},{44,-16},
              {38,-14},{34,-6},{32,4},{34,10},{36,14},{38,16},{36,6},{36,-2},{
              40,-10}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{74,52},{86,-50}},
          lineColor={0,0,0},
          fillColor={95,95,95},
          fillPattern=FillPattern.Forward),
        Line(
          points={{-58,40},{-10,40}},
          color={0,0,255},
          pattern=LinePattern.DashDotDot,
          smooth=Smooth.None,
          arrow={Arrow.Filled,Arrow.Filled},
          thickness=0.5),
        Line(
          points={{26,40},{74,40}},
          color={0,0,255},
          pattern=LinePattern.DashDotDot,
          smooth=Smooth.None,
          arrow={Arrow.Filled,Arrow.Filled},
          thickness=0.5),
        Line(
          points={{-58,-40},{-10,-40}},
          color={0,0,255},
          pattern=LinePattern.DashDotDot,
          smooth=Smooth.None,
          arrow={Arrow.Filled,Arrow.Filled},
          thickness=0.5),
        Line(
          points={{26,-40},{74,-40}},
          color={0,0,255},
          pattern=LinePattern.DashDotDot,
          smooth=Smooth.None,
          arrow={Arrow.Filled,Arrow.Filled},
          thickness=0.5)}),
  Documentation(info="<html>
<p>
This is a model of a heat and moisture conductor with multiple material layers.
The construction has at least one material layer, and each layer has
at least one temperature node. The layers are modeled using an instance of 
<a href=\"Buildings.HeatTransfer.Conduction.SingleLayerHM\">
Buildings.HeatTransfer.Conduction.SingleLayerHM</a>.
</p>
<p>
The construction material is defined by a record of the package
<a href=\"modelica://Buildings.HeatTransfer.Data.OpaqueConstructions\">
Buildings.HeatTransfer.Data.OpaqueConstructions</a>.
This record allows specifying materials that store energy, and material
that are a hygrothermal conductor only with no heat storage.
</p>
<p>
To obtain the surface temperature of the construction, use <code>port_a.T</code> (or <code>port_b.T</code>)
and not the variable <code>T[1]</code> because there is a thermal resistance between the surface
and the temperature state.
</p>
</html>", revisions="<html>
<ul>
<li>
Jun 8 2014, by Antoine Hudzik:<br/>

First implementation.
</li>
</ul>

          </html>"));
end MultiLayerHM;
