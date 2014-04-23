within Buildings.HeatTransfer.Sources;
model PrescribedHumidity

  Modelica.Blocks.Interfaces.RealInput X annotation (Placement(transformation(
          extent={{-140,-20},{-100,20}}, rotation=0)));

  Interfaces.MassPort massPort
    annotation (Placement(transformation(extent={{90,-8},{110,12}})));
equation
  massPort.Xi_outflow = X;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                         graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          pattern=LinePattern.None,
          lineThickness=0.5,
          fillColor={170,170,255},
          fillPattern=FillPattern.Backward),
        Polygon(
          points={{44,-20},{44,20},{82,0},{44,-20}},
          lineColor={0,128,255},
          fillColor={0,128,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{16,-8},{-84,-108}},
          lineColor={0,0,0},
          textString="Kg/Kg (dry air)"),
        Line(
          points={{-64,0},{44,0}},
          color={0,128,255},
          thickness=0.5)}),    Documentation(info="<HTML>
<p>
This model represents a variable humidity ratio boundary condition.
The mass fraction in [Kg/Kg of dry air] is given as input signal <b>X</b>
to the model. The effect is that an instance of this model acts as
an infinite reservoir able to absorb or generate as much water vapour 
as required to keep the humidity at the specified value.
</p>
</HTML>
"), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}}), graphics));
end PrescribedHumidity;
