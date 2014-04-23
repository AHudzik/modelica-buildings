within Buildings.HeatTransfer.Sources;
model FixedHumidity

parameter Modelica.SIunits.MassFraction X "Fixed humidity at port";

  Interfaces.MassPort massPort
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));
equation
  massPort.Xi_outflow = X;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}), graphics={
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
          thickness=0.5)}),
           Documentation(info="<HTML>
<p>
This model defines a fixed humidity ratio at its port in kilogramme per kilogramme of dry air (mass of water vapor present in moist air - to the mass of the dry air),
i.e., it defines a fixed moisture content as a boundary condition.
</p>
</HTML>
"),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
            100,100}}), graphics));
end FixedHumidity;
