within Buildings.HeatTransfer.Interfaces;
connector HeatMassPort_a

  parameter Boolean use_massPort = true;

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort;

   Buildings.HeatTransfer.Interfaces.MassPort massPort if use_massPort;

  annotation (                                 Icon(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
        Polygon(
          points={{0,100},{-72,40},{0,100}},
          pattern=LinePattern.None,
          smooth=Smooth.None,
          fillColor={98,0,98},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{0,100},{0,100}},
          pattern=LinePattern.None,
          smooth=Smooth.None,
          fillColor={98,0,98},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{0,100},{-98,-2},{0,100}},
          pattern=LinePattern.None,
          smooth=Smooth.None,
          fillColor={98,0,98},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{0,100},{0,100}},
          pattern=LinePattern.None,
          smooth=Smooth.None,
          fillColor={98,0,98},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{0,100},{-98,0},{2,-100},{100,0},{0,100}},
          pattern=LinePattern.None,
          smooth=Smooth.None,
          fillColor={98,0,98},
          fillPattern=FillPattern.Solid)}),defaultComponentName = "heatMassPort_a",
    Documentation(info="<HTML>
    <p>This connector is used for 1-dimensional heat and moisture flow between components.
The variables in the connector are:</p>
<pre>
T           Temperature in [K].
Xi_outflow  Mass fraction [kg/kg of dry air].
Q_flow      Heat flow rate in [W].
m_flow      Mass flow rate in [kg/s].
</pre>
<p>According to the Modelica sign convention, a <b>positive</b> heat flow
rate <b>Q_flow</b> and  a <b>positive</b> mass flow
rate <b>m_flow</b> are considered to flow <b>into</b> a component. This
convention has to be used whenever this connector is used in a model
class.</p>
<p>Note, that the two connector classes <b>HeatMassPort_a</b> and
<b>HeatMassPort_b</b> are identical with the only exception of the different
<b>icon layout</b>.</p></html>"));

end HeatMassPort_a;
